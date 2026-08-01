"""
라인 데이터 + DEM → 네트워크 데이터셋용 GPKG

입력: 라인 GPKG (EPSG:5179), DEM GeoTIFF (EPSG:5179)
출력: 16개 비용 필드가 계산된 GPKG (ND.gpkg)

필요 패키지: geopandas, rasterio, shapely, numpy
  pip install geopandas rasterio shapely numpy
"""

import tkinter as tk
from tkinter import filedialog
import geopandas as gpd
import rasterio
import numpy as np
import math
from shapely.geometry import LineString
from shapely.ops import substring
import os
import warnings
import pandas as pd

warnings.filterwarnings("ignore")

# ============================================================
# 파일 선택
# ============================================================
root = tk.Tk()
root.withdraw()

in_lines = filedialog.askopenfilename(
    title="라인 데이터 선택",
    filetypes=[("GeoPackage", "*.gpkg"), ("All", "*.*")]
)
if not in_lines:
    raise SystemExit("라인 데이터 미선택")

in_dem = filedialog.askopenfilename(
    title="DEM 래스터 선택",
    filetypes=[("GeoTIFF", "*.tif *.tiff"), ("All", "*.*")]
)
if not in_dem:
    raise SystemExit("DEM 미선택")

in_water = filedialog.askopenfilename(
    title="수상 폴리곤 선택 (없으면 취소)",
    filetypes=[("GeoPackage", "*.gpkg"), ("All", "*.*")]
)
if in_water:
    print(f"  수상 폴리곤: {in_water} (ID 필드: 바다/하천)")
else:
    print("  수상 폴리곤 없음 → 전 구간 육상 처리")

out_dir = filedialog.askdirectory(title="출력 폴더 선택")
if not out_dir:
    raise SystemExit("출력 폴더 미선택")

# DEM 셀 크기를 분할 간격 기본값으로 사용 (래스터 그래프와 동일 해상도)
with rasterio.open(in_dem) as src:
    if src.crs is None or src.crs.to_epsg() != 5179:
        raise SystemExit(f"DEM CRS 오류: EPSG:5179 필요, 현재 {src.crs}")
    dem_xres, dem_yres = abs(src.res[0]), abs(src.res[1])
if abs(dem_xres - dem_yres) > 0.1:
    print(f"  ⚠ DEM x·y 해상도 불일치: x={dem_xres:.4f}, y={dem_yres:.4f} → x값 사용")
dem_cell_m = dem_xres
print(f"  DEM 셀 크기: {dem_cell_m:.6f}m")

split_input = input(
    f"세그먼트 분할 간격(m, Enter=DEM 셀 크기 {dem_cell_m:.4f}, 0=분할 안 함): "
).strip()
split_interval = float(split_input) if split_input else dem_cell_m

# 수상 구간 속도 (km/h)
V_SEA = 5.0
V_RIVER = 1.2
LAMBDA_DOWN = 0.5


# ============================================================
# 벡터 속도/MET 함수 (NumPy)
# ============================================================
def v_kmh_tobler_vec(g):
    return np.clip(6.0 * np.exp(-3.5 * np.abs(g + 0.05)), 1.0, 7.5)


def v_kmh_kondoseino_vec(g):
    return np.maximum(
        np.where(
            g >= -0.07,
            5.1 * np.exp(-2.25 * np.abs(g + 0.07)),
            5.1 * np.exp(-1.5 * np.abs(g + 0.07)),
        ),
        0.5,
    )


def met_acsm_vec(v_kmh, g):
    v_mpm = v_kmh * 1000.0 / 60.0
    g_eff = np.where(g >= 0, g, LAMBDA_DOWN * g)
    vo2 = 0.1 * v_mpm + 1.8 * v_mpm * g_eff + 3.5
    return np.maximum(vo2 / 3.5, 1.0)


# ============================================================
# Step 1: 라인 분할
# ============================================================
def split_line_at_interval(line, interval):
    if interval <= 0 or line.length <= interval:
        return [line]
    segments = []
    start = 0.0
    total = line.length
    while start < total:
        end = min(start + interval, total)
        seg = substring(line, start, end)
        if seg.length > 0:
            segments.append(seg)
        start = end
    return segments


def split_all_lines(gdf, interval):
    if interval <= 0:
        return gdf.copy()

    attr_cols = [c for c in gdf.columns if c != "geometry"]
    # 속성 행을 한 번만 dict로 변환
    attr_rows = gdf[attr_cols].to_dict("records")

    new_geoms = []
    new_attr_idx = []  # 원본 행 인덱스 → 속성 재사용
    total = len(gdf)

    for i, geom in enumerate(gdf.geometry):
        lines = list(geom.geoms) if geom.geom_type == "MultiLineString" else [geom]
        segs = []
        for ln in lines:
            segs.extend(split_line_at_interval(ln, interval))
        for seg in segs:
            new_geoms.append(seg)
            new_attr_idx.append(i)
        if (i + 1) % 1000 == 0:
            print(f"  분할 [{i + 1}/{total}]")

    # 속성 복원: 인덱스 기반 슬라이싱
    attr_df = pd.DataFrame([attr_rows[j] for j in new_attr_idx])
    attr_df["geometry"] = new_geoms
    result = gpd.GeoDataFrame(attr_df, crs=gdf.crs)
    print(f"  분할 완료: {total}개 → {len(result)}개 세그먼트")
    return result


# ============================================================
# Step 2: DEM Z값 샘플링 (일괄)
# ============================================================
def interpolate_z(gdf, dem_path):
    total = len(gdf)

    # 1) 전체 좌표 수집 + 오프셋 기록
    all_xy = []
    offsets = [0]
    for geom in gdf.geometry:
        coords = list(geom.coords)
        all_xy.extend([(c[0], c[1]) for c in coords])
        offsets.append(len(all_xy))

    print(f"  좌표 수집 완료: {len(all_xy)}개 점")

    # 2) DEM 일괄 샘플링 (1회 호출)
    with rasterio.open(dem_path) as src:
        dem_nodata = src.nodata
        z_all = np.array([v[0] for v in src.sample(all_xy)], dtype=np.float64)

    # nodata 처리
    mask = np.isnan(z_all)
    if dem_nodata is not None:
        mask |= (z_all == dem_nodata)
    z_all[mask] = 0.0

    print(f"  DEM 샘플링 완료: {len(z_all)}개 Z값")

    # 3) geometry 재구축
    new_geoms = []
    for i, geom in enumerate(gdf.geometry):
        coords_2d = list(geom.coords)
        z_slice = z_all[offsets[i]:offsets[i + 1]]
        coords_3d = [
            (c[0], c[1], float(z))
            for c, z in zip(coords_2d, z_slice)
        ]
        new_geoms.append(LineString(coords_3d))
        if (i + 1) % 50000 == 0:
            print(f"  geometry 재구축 [{i + 1}/{total}]")

    gdf = gdf.copy()
    gdf["geometry"] = new_geoms
    print(f"  Z 보간 완료: {total}개 라인")
    return gdf


# ============================================================
# Step 2-b: 수상 폴리곤 공간 조인 → water_id 필드 부여
# ============================================================
def assign_water_id(gdf, water_path):
    water = gpd.read_file(water_path)
    if water.crs != gdf.crs:
        water = water.to_crs(gdf.crs)

    water = water.rename(columns={"ID": "_water_id"})[["_water_id", "geometry"]]

    if "water_id" in gdf.columns:
        gdf = gdf.drop(columns=["water_id"])

    pts = gdf.copy()
    pts["geometry"] = pts.geometry.interpolate(0.5, normalized=True)

    joined = gpd.sjoin(pts, water, how="left", predicate="within")
    joined = joined[~joined.index.duplicated(keep="first")]

    gdf = gdf.copy()
    gdf["water_id"] = joined["_water_id"].reindex(gdf.index).fillna("").values

    cnt_sea = (gdf["water_id"] == "바다").sum()
    cnt_riv = (gdf["water_id"] == "하천").sum()
    cnt_land = len(gdf) - cnt_sea - cnt_riv
    print(f"  수상 판별 완료: 육상 {cnt_land}, 바다 {cnt_sea}, 하천 {cnt_riv}")
    return gdf


# ============================================================
# Step 3: 3D 길이 계산
# ============================================================
def calc_length_3d(geom):
    coords = list(geom.coords)
    total = 0.0
    for j in range(1, len(coords)):
        dx = coords[j][0] - coords[j - 1][0]
        dy = coords[j][1] - coords[j - 1][1]
        dz = coords[j][2] - coords[j - 1][2] if len(coords[j]) > 2 else 0.0
        total += math.sqrt(dx * dx + dy * dy + dz * dz)
    return total


# ============================================================
# Step 4: 전체 필드 계산 (벡터 연산)
# ============================================================
def calculate_fields(gdf):
    total = len(gdf)
    print(f"필드 계산 시작 (총 {total}개 세그먼트)")

    # --- 1) geometry에서 배열 추출 (유일한 루프) ---
    l3d_arr = np.empty(total)
    z_from_arr = np.empty(total)
    z_to_arr = np.empty(total)

    for i, geom in enumerate(gdf.geometry):
        l3d_arr[i] = calc_length_3d(geom)
        coords = list(geom.coords)
        z_from_arr[i] = coords[0][2] if len(coords[0]) > 2 else 0.0
        z_to_arr[i] = coords[-1][2] if len(coords[-1]) > 2 else 0.0
        if (i + 1) % 50000 == 0:
            print(f"  좌표 추출 [{i + 1}/{total}]")

    l2d_arr = gdf.geometry.length.values

    # l3d 보정: 0 이하이면 l2d 사용
    zero_mask = l3d_arr <= 0
    l3d_arr[zero_mask] = l2d_arr[zero_mask]
    l3dkm = l3d_arr / 1000.0

    # --- 2) 경사도 ---
    safe_l2d = np.where(l2d_arr > 0, l2d_arr, 1.0)  # 0 나눗셈 방지
    g_f = (z_to_arr - z_from_arr) / safe_l2d
    g_f[l2d_arr <= 0] = 0.0
    g_b = -g_f

    # --- 3) 속도 (벡터) ---
    sp_tob_f = v_kmh_tobler_vec(g_f)
    sp_tob_b = v_kmh_tobler_vec(g_b)
    sp_ks_f = v_kmh_kondoseino_vec(g_f)
    sp_ks_b = v_kmh_kondoseino_vec(g_b)

    # --- 4) 수상 구간 속도 덮어쓰기 ---
    has_id = "water_id" in gdf.columns
    if has_id:
        id_vals = gdf["water_id"].fillna("").str.strip().values
        sea_mask = id_vals == "바다"
        riv_mask = id_vals == "하천"
        water_mask = sea_mask | riv_mask

        sp_tob_f[sea_mask] = V_SEA
        sp_tob_b[sea_mask] = V_SEA
        sp_ks_f[sea_mask] = V_SEA
        sp_ks_b[sea_mask] = V_SEA
        sp_tob_f[riv_mask] = V_RIVER
        sp_tob_b[riv_mask] = V_RIVER
        sp_ks_f[riv_mask] = V_RIVER
        sp_ks_b[riv_mask] = V_RIVER
    else:
        water_mask = np.zeros(total, dtype=bool)

    # --- 5) 시간 ---
    hr_tob_f = np.where(sp_tob_f > 0, l3dkm / sp_tob_f, 1e6)
    hr_tob_b = np.where(sp_tob_b > 0, l3dkm / sp_tob_b, 1e6)
    hr_ks_f = np.where(sp_ks_f > 0, l3dkm / sp_ks_f, 1e6)
    hr_ks_b = np.where(sp_ks_b > 0, l3dkm / sp_ks_b, 1e6)

    # --- 6) 에너지 소비량 ---
    kcal_ks_f = met_acsm_vec(sp_ks_f, g_f) * hr_ks_f
    kcal_ks_b = met_acsm_vec(sp_ks_b, g_b) * hr_ks_b
    kcal_tob_f = met_acsm_vec(sp_tob_f, g_f) * hr_tob_f
    kcal_tob_b = met_acsm_vec(sp_tob_b, g_b) * hr_tob_b

    # 수상 구간 에너지 = 0
    kcal_ks_f[water_mask] = 0.0
    kcal_ks_b[water_mask] = 0.0
    kcal_tob_f[water_mask] = 0.0
    kcal_tob_b[water_mask] = 0.0

    # --- 7) GeoDataFrame에 일괄 할당 ---
    gdf = gdf.copy()
    gdf["length_3d"] = l3d_arr
    gdf["length_3dkm"] = l3dkm
    gdf["grade_f"] = g_f
    gdf["grade_b"] = g_b
    gdf["sp_tob_f"] = sp_tob_f
    gdf["sp_tob_b"] = sp_tob_b
    gdf["sp_ks_f"] = sp_ks_f
    gdf["sp_ks_b"] = sp_ks_b
    gdf["hr_tob_f"] = hr_tob_f
    gdf["hr_tob_b"] = hr_tob_b
    gdf["hr_ks_f"] = hr_ks_f
    gdf["hr_ks_b"] = hr_ks_b
    gdf["kcal_ks_f"] = kcal_ks_f
    gdf["kcal_ks_b"] = kcal_ks_b
    gdf["kcal_tob_f"] = kcal_tob_f
    gdf["kcal_tob_b"] = kcal_tob_b

    print("필드 계산 완료")
    return gdf


# ============================================================
# QA 검증 (벡터)
# ============================================================
def run_qa(gdf):
    l2d = gdf.geometry.length.values
    l3d = gdf["length_3d"].values
    g_f = gdf["grade_f"].values
    hr_tob_f = gdf["hr_tob_f"].values
    hr_tob_b = gdf["hr_tob_b"].values

    err1 = np.sum((l3d > 0) & (l2d > 0) & (l3d < l2d * 0.999))
    err2 = np.sum((np.abs(g_f) > 0.01) & (g_f > 0) & (hr_tob_f < hr_tob_b * 0.99))
    err = int(err1 + err2)

    print(f"QA 완료: {len(gdf)}개 검사, {err}개 이상치")
    if err > 0:
        print(f"  ⚠ {err}개 세그먼트에서 예상 범위 벗어남 (수상구간 포함 가능)")


# ============================================================
# 메인
# ============================================================
def main():
    print("=" * 60)
    print("라인 데이터 읽기")
    print("=" * 60)
    gdf = gpd.read_file(in_lines)
    print(f"  {len(gdf)}개 피처, CRS: {gdf.crs}")

    if gdf.crs is None or gdf.crs.to_epsg() != 5179:
        raise SystemExit(f"CRS 오류: EPSG:5179 필요, 현재 {gdf.crs}")

    if split_interval > 0:
        print("=" * 60)
        print(f"Step 1: 라인 분할 (간격: {split_interval:.4f}m)")
        print("=" * 60)
        gdf = split_all_lines(gdf, split_interval)
    else:
        print("  분할 없음")

    print("=" * 60)
    print("Step 2: DEM Z값 보간")
    print("=" * 60)
    gdf = interpolate_z(gdf, in_dem)

    if in_water:
        print("=" * 60)
        print("Step 2-b: 수상 폴리곤 공간 조인 → water_id 필드 부여")
        print("=" * 60)
        gdf = assign_water_id(gdf, in_water)

    print("=" * 60)
    print("Step 3: 속성 필드 계산")
    print("=" * 60)
    gdf = calculate_fields(gdf)

    run_qa(gdf)

    print("=" * 60)
    print("저장")
    print("=" * 60)
    out_path = os.path.join(out_dir, "nd.gpkg")
    gdf.to_file(out_path, driver="GPKG")
    print(f"완료: {out_path}")


if __name__ == "__main__":
    main()