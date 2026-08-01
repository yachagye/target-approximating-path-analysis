"""
라인 피처 비용 필드 부여 (위상·스키마 유지)

입력: 라인 GPKG (EPSG:5179), DEM GeoTIFF (EPSG:5179), [수상 폴리곤 GPKG]
출력: 원본 위상·스키마를 유지하고 9개 비용 필드를 추가한 GPKG ({입력명}_cost.gpkg)

추가 필드:
  start_x, start_y, end_x, end_y  (WGS84 경위도, 정방향 기준: 시작=첫 정점, 종료=마지막 정점)
  length_3dkm
  hr_ks_f, hr_ks_b, hr_tob_f, hr_tob_b
  kcal_ks_f, kcal_ks_b, kcal_tob_f, kcal_tob_b

분할 없이 라인 내부 정점 구간(vertex-to-vertex)별로
경사 → 속도 → 시간/에너지를 계산하고 피처 단위로 합산한다.
geometry는 원본 그대로 출력한다(위상·스키마 보존).
속도/MET 함수·상수는 라인_네트워크데이터셋_변환_gpkg.py와 동일하다.

필요 패키지: geopandas, rasterio, shapely, numpy
"""

import os
import warnings
import tkinter as tk
from tkinter import filedialog

import numpy as np
import geopandas as gpd
import rasterio
from pyproj import Transformer

warnings.filterwarnings("ignore")

# 수상 구간 속도 (km/h) — 네트워크 변환 스크립트와 동일
V_SEA = 5.0
V_RIVER = 1.2
LAMBDA_DOWN = 0.5


# ── 속도/MET 함수 (네트워크 변환 스크립트와 동일) ──────────────
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


# ── 정방향 시작·종료 정점 추출 (단일·멀티파트 모두 대응) ───────
def first_coord(geom):
    if geom.geom_type == "MultiLineString":
        return geom.geoms[0].coords[0]
    return geom.coords[0]


def last_coord(geom):
    if geom.geom_type == "MultiLineString":
        return geom.geoms[-1].coords[-1]
    return geom.coords[-1]


# ── 정점·구간 평탄화 ───────────────────────────────────────────
def flatten_segments(gdf):
    """피처별 정점을 모으고, 파트 내부 연속 정점 쌍만 구간으로 평탄화.

    MultiLineString은 파트 경계를 넘는 구간을 만들지 않는다.
    반환: vert_x, vert_y, seg_from, seg_to, seg_feat (모두 ndarray)
    """
    vert_x, vert_y = [], []
    seg_from, seg_to, seg_feat = [], [], []
    vi = 0
    total = len(gdf)
    for fi, geom in enumerate(gdf.geometry):
        parts = list(geom.geoms) if geom.geom_type == "MultiLineString" else [geom]
        for part in parts:
            coords = list(part.coords)
            start = vi
            for c in coords:
                vert_x.append(c[0])
                vert_y.append(c[1])
                vi += 1
            for k in range(len(coords) - 1):
                seg_from.append(start + k)
                seg_to.append(start + k + 1)
                seg_feat.append(fi)
        if (fi + 1) % 1000 == 0:
            print(f"  평탄화 [{fi + 1}/{total}]")
    return (
        np.asarray(vert_x), np.asarray(vert_y),
        np.asarray(seg_from), np.asarray(seg_to), np.asarray(seg_feat),
    )


# ── DEM Z 샘플링 (일괄) ────────────────────────────────────────
def sample_z(vert_x, vert_y, dem_path):
    with rasterio.open(dem_path) as src:
        nodata = src.nodata
        z = np.array(
            [v[0] for v in src.sample(np.column_stack([vert_x, vert_y]))],
            dtype=np.float64,
        )
    mask = np.isnan(z)
    if nodata is not None:
        mask |= (z == nodata)
    z[mask] = 0.0
    return z


# ── 수상 판별 (구간 중점 공간 조인) ────────────────────────────
def assign_water(mid_x, mid_y, water_path, crs):
    water = gpd.read_file(water_path)
    if water.crs != crs:
        water = water.to_crs(crs)
    water = water.rename(columns={"ID": "_wid"})[["_wid", "geometry"]]

    pts = gpd.GeoDataFrame(
        {"_sid": np.arange(len(mid_x))},
        geometry=gpd.points_from_xy(mid_x, mid_y),
        crs=crs,
    )
    joined = gpd.sjoin(pts, water, how="left", predicate="within")
    joined = joined[~joined.index.duplicated(keep="first")]
    return joined["_wid"].reindex(range(len(mid_x))).fillna("").values


# ── main ───────────────────────────────────────────────────────
def main():
    root = tk.Tk()
    root.withdraw()

    in_lines = filedialog.askopenfilename(
        title="라인 데이터 선택",
        filetypes=[("GeoPackage", "*.gpkg"), ("All", "*.*")],
    )
    if not in_lines:
        raise SystemExit("라인 데이터 미선택")

    in_dem = filedialog.askopenfilename(
        title="DEM 래스터 선택",
        filetypes=[("GeoTIFF", "*.tif *.tiff"), ("All", "*.*")],
    )
    if not in_dem:
        raise SystemExit("DEM 미선택")

    in_water = filedialog.askopenfilename(
        title="수상 폴리곤 선택 (없으면 취소)",
        filetypes=[("GeoPackage", "*.gpkg"), ("All", "*.*")],
    )
    if in_water:
        print(f"수상 폴리곤: {in_water} (ID 필드: 바다/하천)")
    else:
        print("수상 폴리곤 없음 → 전 구간 육상 처리")

    out_dir = filedialog.askdirectory(title="출력 폴더 선택")
    if not out_dir:
        raise SystemExit("출력 폴더 미선택")

    # 입력 읽기 + CRS 검증
    gdf = gpd.read_file(in_lines)
    print(f"입력: {len(gdf)}개 피처, CRS: {gdf.crs}")
    if gdf.crs is None or gdf.crs.to_epsg() != 5179:
        raise SystemExit(f"CRS 오류: EPSG:5179 필요, 현재 {gdf.crs}")
    with rasterio.open(in_dem) as src:
        if src.crs is None or src.crs.to_epsg() != 5179:
            raise SystemExit(f"DEM CRS 오류: EPSG:5179 필요, 현재 {src.crs}")

    # 1) 정점·구간 평탄화
    print("정점·구간 평탄화")
    vx, vy, s_from, s_to, s_feat = flatten_segments(gdf)
    print(f"  정점 {len(vx)}개, 구간 {len(s_feat)}개")

    # 2) DEM Z 샘플링
    print("DEM Z 샘플링")
    z = sample_z(vx, vy, in_dem)

    # 3) 구간 기하 (벡터)
    x0, y0, z0 = vx[s_from], vy[s_from], z[s_from]
    x1, y1, z1 = vx[s_to], vy[s_to], z[s_to]
    dx, dy, dz = x1 - x0, y1 - y0, z1 - z0
    l2d = np.sqrt(dx * dx + dy * dy)
    l3d = np.sqrt(dx * dx + dy * dy + dz * dz)
    zero = l3d <= 0
    l3d[zero] = l2d[zero]
    l3dkm = l3d / 1000.0

    # 4) 경사 (forward = 시점→종점, backward = 반대)
    safe_l2d = np.where(l2d > 0, l2d, 1.0)
    g_f = dz / safe_l2d
    g_f[l2d <= 0] = 0.0
    g_b = -g_f

    # 5) 속도
    sp_tob_f = v_kmh_tobler_vec(g_f)
    sp_tob_b = v_kmh_tobler_vec(g_b)
    sp_ks_f = v_kmh_kondoseino_vec(g_f)
    sp_ks_b = v_kmh_kondoseino_vec(g_b)

    # 6) 수상 구간 속도 덮어쓰기
    water_mask = np.zeros(len(s_feat), dtype=bool)
    if in_water:
        print("수상 구간 판별")
        mid_x, mid_y = (x0 + x1) / 2.0, (y0 + y1) / 2.0
        wid = assign_water(mid_x, mid_y, in_water, gdf.crs)
        sea = np.array([str(w).strip() == "바다" for w in wid])
        riv = np.array([str(w).strip() == "하천" for w in wid])
        water_mask = sea | riv
        for sp in (sp_tob_f, sp_tob_b, sp_ks_f, sp_ks_b):
            sp[sea] = V_SEA
            sp[riv] = V_RIVER
        print(f"  육상 {(~water_mask).sum()}, 바다 {sea.sum()}, 하천 {riv.sum()}")

    # 7) 시간 (속도는 항상 양수: Tobler≥1.0, Kondo–Seino≥0.5, 수상≥1.2)
    hr_tob_f = l3dkm / sp_tob_f
    hr_tob_b = l3dkm / sp_tob_b
    hr_ks_f = l3dkm / sp_ks_f
    hr_ks_b = l3dkm / sp_ks_b

    # 8) 에너지 (수상 구간은 0)
    kcal_ks_f = met_acsm_vec(sp_ks_f, g_f) * hr_ks_f
    kcal_ks_b = met_acsm_vec(sp_ks_b, g_b) * hr_ks_b
    kcal_tob_f = met_acsm_vec(sp_tob_f, g_f) * hr_tob_f
    kcal_tob_b = met_acsm_vec(sp_tob_b, g_b) * hr_tob_b
    for arr in (kcal_ks_f, kcal_ks_b, kcal_tob_f, kcal_tob_b):
        arr[water_mask] = 0.0

    # 9) 피처 단위 합산
    n = len(gdf)

    def feat_sum(a):
        return np.bincount(s_feat, weights=a, minlength=n)

    gdf = gdf.copy()
    sx = gdf.geometry.apply(lambda g: first_coord(g)[0]).to_numpy()
    sy = gdf.geometry.apply(lambda g: first_coord(g)[1]).to_numpy()
    ex = gdf.geometry.apply(lambda g: last_coord(g)[0]).to_numpy()
    ey = gdf.geometry.apply(lambda g: last_coord(g)[1]).to_numpy()
    transformer = Transformer.from_crs("EPSG:5179", "EPSG:4326", always_xy=True)
    gdf["start_x"], gdf["start_y"] = transformer.transform(sx, sy)
    gdf["end_x"], gdf["end_y"] = transformer.transform(ex, ey)
    gdf["length_3dkm"] = feat_sum(l3dkm)
    gdf["hr_ks_f"] = feat_sum(hr_ks_f)
    gdf["hr_ks_b"] = feat_sum(hr_ks_b)
    gdf["hr_tob_f"] = feat_sum(hr_tob_f)
    gdf["hr_tob_b"] = feat_sum(hr_tob_b)
    gdf["kcal_ks_f"] = feat_sum(kcal_ks_f)
    gdf["kcal_ks_b"] = feat_sum(kcal_ks_b)
    gdf["kcal_tob_f"] = feat_sum(kcal_tob_f)
    gdf["kcal_tob_b"] = feat_sum(kcal_tob_b)

    # 10) 저장 (geometry 원본 유지)
    base = os.path.splitext(os.path.basename(in_lines))[0]
    out_path = os.path.join(out_dir, f"{base}_cost.gpkg")
    gdf.to_file(out_path, driver="GPKG")
    print(f"완료: {out_path} ({n}개 피처, 좌표 4 + 비용 9 필드 추가)")


if __name__ == "__main__":
    main()