# 목표값 ±margin 편차범위 내 상위 N위 경로 분석 (포인트-포인트 / 경유지 포함)
#
# 입력:
#  1) graph_cache.pkl  (네트워크데이터셋_DiGraph_변환_gpkg_pkl.py 출력물: payload dict)
#  2) 경로 분석.csv
#       route_id, x1, y1, x2, y2, ..., xN, yN, (km_beg, km_end) 또는 (hr_beg, hr_end)
#     - 좌표: EPSG:4326 (경도 x, 위도 y)
#     - 목표값 컬럼은 beg, end 두 개:
#         · beg만 입력(end 빈칸)  → 단일 목표값 T = beg → 편차범위 [T-margin, T+margin]
#         · beg, end 모두 입력     → 구간 [beg, end] → 편차범위 [beg-margin, end+margin]
#         · beg 빈칸               → 목표값 미입력 → 건너뜀
#         · beg > end, 또는 end만 입력 → 오류로 해당 행 건너뜀
#  3) 선택 입력: 장애물 GPKG (line / polygon, 완전 차단)
#
# margin:
#  - 거리 기반(km): ±2.25 km 고정 (5리 × 0.45 km)
#  - 시간 기반(hr): 실행 시 사용자 입력 (기본 0.5 hour)
#
# 알고리즘 (편차범위 = 입력 구간을 margin만큼 확장하여 구간 모드로 처리):
#  - 단일 목표값 T는 [T-margin, T+margin], 구간 [beg,end]는 [beg-margin, end+margin]으로
#    확장한 뒤, 확장 구간에 드는 경로를 metric 오름차순 수집 → 다양성 필터 → 상위 TOP_N개.
#  - margin 범위 내 후보가 없으면 해당 경로는 출력하지 않는다(fallback 없음).
#  - abs_diff는 원본 목표값 기준: 단일은 |비용-T|, 구간은 max(0, beg-비용, 비용-end).
#  [한계] 경유지 2개 이상(3-leg+)은 abs_diff/metric 순 상위 N의 엄밀성이 보장되지 않는
#         근사임(순위 분석과 동일). rank 경계에서 어긋날 수 있음.
#
# 출력:
#  - 경로_point_point_편차범위.gpkg (EPSG:5179)
#     layer: route_target  (rank·margin 필드 포함, 경로당 최대 TOP_N개. route_min_* 미산출)

import os
import csv
import math
import pickle
import heapq
import bisect as _bisect
import time
import tkinter as tk
from tkinter import filedialog

import networkx as nx
import geopandas as gpd
import fiona
from shapely.geometry import LineString, Point
from shapely.ops import substring
from pyproj import Transformer


# 후보 생성 파라미터
K_INIT = 10
K_MAX = 30000

# 편차범위 출력 파라미터
TOP_N_INTERVAL = 30        # 편차범위 최종 출력 수
SIM_THRESHOLD = 0.95       # 다양성 필터: 이 값 이상이면 미세 변형으로 간주하여 제외
MARGIN_KM = 2.25           # 거리 기반 편차범위 ±km (5리 × 0.45 km, 고정)
MARGIN_HR = 0.5            # 시간 기반 편차범위 ±hour 기본값(실행 시 입력 가능)

# 2-leg probe 진단 파라미터
K_PROBE = 300              # 각 leg 기울기 측정용 짧은 열거 수
PROBE_RATIO = 3.0          # 두 leg 기울기 비율이 이 값 이상이면 가파른 쪽을 bisect로

# 노드 좌표 라운딩 (graph_cache.pkl 생성 시 동일 값 사용)
NODE_ROUND_M = 0.01

def _round_xy(x, y, r=NODE_ROUND_M):
    return (round(x / r) * r, round(y / r) * r)


def snap_to_graph(x, y, edge_gdf, G):
    """좌표(EPSG:5179)를 최근접 간선에 투영하고, 필요 시 임시 노드를 삽입한다.

    Returns: (node, snap_dist_m, restore_info_or_None)
    """
    pt = Point(x, y)
    dists = edge_gdf.geometry.distance(pt)

    # 거리순 상위 후보 중 그래프에 존재하는 간선 선택 (이전 분할로 제거된 간선 회피)
    for nearest_idx in dists.nsmallest(10).index:
        row = edge_gdf.loc[nearest_idx]
        u, v = row["u"], row["v"]
        if not G.has_edge(u, v):
            continue

        snap_dist = float(dists.loc[nearest_idx])
        edge_geom = row.geometry

        # 투영점 계산
        d = edge_geom.project(pt)
        proj_pt = edge_geom.interpolate(d)

        # 기존 노드와의 거리 확인 (30m 세그먼트 반경 이내면 기존 노드 사용)
        node_tol = 15.0
        dist_u = math.sqrt((u[0] - proj_pt.x) ** 2 + (u[1] - proj_pt.y) ** 2)
        dist_v = math.sqrt((v[0] - proj_pt.x) ** 2 + (v[1] - proj_pt.y) ** 2)

        if dist_u <= node_tol:
            return u, snap_dist, None
        if dist_v <= node_tol:
            return v, snap_dist, None

        # 간선 내부: 임시 노드 삽입
        temp_node = _round_xy(proj_pt.x, proj_pt.y)
        if G.has_node(temp_node):
            return temp_node, snap_dist, None

        restore = _split_edge(G, u, v, temp_node, d)
        return temp_node, snap_dist, restore

    # fallback: 그래프 노드 직접 탐색
    best_node = None
    best_d2 = float("inf")
    for node in G.nodes():
        d2 = (node[0] - x) ** 2 + (node[1] - y) ** 2
        if d2 < best_d2:
            best_d2 = d2
            best_node = node
    return best_node, math.sqrt(best_d2) if best_node else float("inf"), None


def _split_edge(G, u, v, temp_node, d):
    """간선 (u,v)와 역방향 (v,u)를 temp_node 위치에서 분할한다.

    d: u→v 간선 geometry 위 투영 거리 (m).
    비용은 chain_costs를 기반으로 세그먼트별 선형 보간한다.
    Returns: 복원 정보 dict
    """
    cost_keys = ("length_3dkm", "hour_ks", "hour_tob", "kcal_ks", "kcal_tob")

    fwd = dict(G[u][v])
    has_rev = G.has_edge(v, u)
    rev = dict(G[v][u]) if has_rev else None

    # --- 정방향 (u → v) 분할 ---
    fwd_geom = fwd["geom"]
    chain_nodes = fwd.get("chain_nodes", [u, v])
    chain_costs = fwd.get("chain_costs", [{k: float(fwd[k]) for k in cost_keys}])

    costs_ut, costs_tv, geom_ut, geom_tv = _interpolate_split(
        fwd_geom, chain_nodes, chain_costs, d, cost_keys,
    )

    # --- 역방향 (v → u) 분할 ---
    if has_rev:
        rev_geom = rev["geom"]
        rev_chain = rev.get("chain_nodes", [v, u])
        rev_costs = rev.get("chain_costs", [{k: float(rev[k]) for k in cost_keys}])

        d_rev = rev_geom.project(Point(temp_node))
        costs_vt, costs_tu, geom_vt, geom_tu = _interpolate_split(
            rev_geom, rev_chain, rev_costs, d_rev, cost_keys,
        )

    # --- 그래프 수정 ---
    G.remove_edge(u, v)
    if has_rev:
        G.remove_edge(v, u)

    G.add_node(temp_node)
    G.add_edge(u, temp_node, **costs_ut, geom=geom_ut)
    G.add_edge(temp_node, v, **costs_tv, geom=geom_tv)

    if has_rev:
        G.add_edge(v, temp_node, **costs_vt, geom=geom_vt)
        G.add_edge(temp_node, u, **costs_tu, geom=geom_tu)

    return {
        "temp_node": temp_node,
        "u": u, "v": v,
        "fwd_data": fwd,
        "has_rev": has_rev,
        "rev_data": rev,
    }


def _interpolate_split(edge_geom, chain_nodes, chain_costs, d, cost_keys):
    """간선 geometry의 거리 d에서 비용·geometry를 분할한다.

    Returns: (costs_before, costs_after, geom_before, geom_after)
    """
    # 체인 노드별 누적 거리
    cum = [0.0]
    for i in range(1, len(chain_nodes)):
        cum.append(edge_geom.project(Point(chain_nodes[i])))

    # 분할 세그먼트 인덱스
    seg_i = len(chain_costs) - 1
    for i in range(len(cum) - 1):
        if d <= cum[i + 1] + 0.01:
            seg_i = i
            break

    span = cum[seg_i + 1] - cum[seg_i]
    f = (d - cum[seg_i]) / span if span > 0 else 0.0
    f = max(0.0, min(1.0, f))

    # 비용 보간
    costs_before = {}
    costs_after = {}
    for k in cost_keys:
        before_full = sum(chain_costs[j][k] for j in range(seg_i))
        seg_cost = chain_costs[seg_i][k]
        after_full = sum(chain_costs[j][k] for j in range(seg_i + 1, len(chain_costs)))
        costs_before[k] = before_full + f * seg_cost
        costs_after[k] = (1 - f) * seg_cost + after_full

    # geometry 분할
    geom_before = substring(edge_geom, 0, d, normalized=False)
    geom_after = substring(edge_geom, d, edge_geom.length, normalized=False)

    # substring이 Point를 반환하는 극단 케이스 방어
    split_pt = edge_geom.interpolate(d)
    if geom_before.geom_type == "Point":
        geom_before = LineString([geom_before.coords[0], split_pt.coords[0]])
    if geom_after.geom_type == "Point":
        geom_after = LineString([split_pt.coords[0], geom_after.coords[0]])

    return costs_before, costs_after, geom_before, geom_after


def restore_splits(G, infos):
    """임시 노드와 분할 간선을 제거하고 원본 간선을 복원한다."""
    for info in infos:
        if info is None:
            continue
        temp = info["temp_node"]
        u, v = info["u"], info["v"]

        for succ in list(G.successors(temp)):
            G.remove_edge(temp, succ)
        for pred in list(G.predecessors(temp)):
            G.remove_edge(pred, temp)
        G.remove_node(temp)

        G.add_edge(u, v, **info["fwd_data"])
        if info["has_rev"]:
            G.add_edge(v, u, **info["rev_data"])


def parse_routes_csv(csv_path, mode):
    rows = []

    def _read_with_encoding(encoding_name):
        out = []
        with open(csv_path, "r", encoding=encoding_name, newline="") as f:
            reader = csv.reader(f)
            try:
                header = next(reader)
            except StopIteration:
                raise RuntimeError("CSV가 비어 있습니다.")

            if not header:
                raise RuntimeError("CSV 헤더가 없습니다.")

            beg_col = "km_beg" if mode == "km" else "hr_beg"
            end_col = "km_end" if mode == "km" else "hr_end"
            if header[0].strip() != "route_id":
                raise RuntimeError("첫 컬럼은 route_id 이어야 합니다.")
            if header[-2].strip() != beg_col or header[-1].strip() != end_col:
                raise RuntimeError(f"마지막 두 컬럼은 {beg_col}, {end_col} 이어야 합니다.")

            coord_cols = header[1:-2]
            if len(coord_cols) < 4 or len(coord_cols) % 2 != 0:
                raise RuntimeError("좌표 컬럼 구조가 올바르지 않습니다.")

            for row in reader:
                if not row:
                    continue
                if len(row) < len(header):
                    row = row + [""] * (len(header) - len(row))
                elif len(row) > len(header):
                    row = row[:len(header)]

                route_id = str(row[0]).strip()
                beg_str = str(row[-2]).strip()
                end_str = str(row[-1]).strip()
                if route_id == "":
                    continue

                coords = []
                for i in range(1, len(row) - 2, 2):
                    xs = str(row[i]).strip()
                    ys = str(row[i + 1]).strip()
                    if xs == "" or ys == "":
                        continue
                    coords.append((float(xs), float(ys)))

                if len(coords) < 2:
                    continue

                # 목표값 유효성: end만 입력된 행은 오류로 건너뜀
                if beg_str == "" and end_str != "":
                    print(f"  [skip] {route_id}: end만 입력됨 (beg 누락)")
                    continue

                beg = float(beg_str) if beg_str != "" else None
                end = float(end_str) if end_str != "" else None

                # beg > end 행은 오류로 건너뜀
                if beg is not None and end is not None and beg > end:
                    print(f"  [skip] {route_id}: beg({beg}) > end({end})")
                    continue

                out.append({
                    "route_id": route_id,
                    "coords": coords,
                    "beg": beg,
                    "end": end,
                })
        return out

    try:
        rows = _read_with_encoding("utf-8-sig")
    except UnicodeDecodeError:
        rows = _read_with_encoding("cp949")

    if not rows:
        raise RuntimeError("CSV에서 유효한 경로 데이터를 찾지 못했습니다.")

    seen = set()
    dup = set()
    for row in rows:
        rid = row["route_id"]
        if rid in seen:
            dup.add(rid)
        else:
            seen.add(rid)

    if dup:
        raise RuntimeError(f"CSV에 중복 route_id 존재: {sorted(list(dup))[:10]}")

    return rows


def accumulate_costs(G, path_nodes):
    s = {
        "length_km": 0.0,
        "hour_ks": 0.0,
        "hour_tob": 0.0,
        "kcal_ks": 0.0,
        "kcal_tob": 0.0,
    }
    for a, b in zip(path_nodes[:-1], path_nodes[1:]):
        ed = G[a][b]
        s["length_km"] += float(ed["length_3dkm"])
        s["hour_ks"] += float(ed["hour_ks"])
        s["hour_tob"] += float(ed["hour_tob"])
        s["kcal_ks"] += float(ed["kcal_ks"])
        s["kcal_tob"] += float(ed["kcal_tob"])
    return s


def path_metric_from_costs(costs, weight_attr):
    if weight_attr == "length_3dkm":
        return float(costs["length_km"])
    return float(costs[weight_attr])


def path_to_linestring(G, path_nodes):
    coords = []
    for a, b in zip(path_nodes[:-1], path_nodes[1:]):
        geom = G[a][b].get("geom", None)
        if geom is None or geom.is_empty:
            ax, ay = a
            bx, by = b
            seg = [(ax, ay), (bx, by)]
        else:
            seg = list(geom.coords)

        if not coords:
            coords.extend(seg)
        else:
            if coords[-1] == seg[0]:
                coords.extend(seg[1:])
            else:
                coords.extend(seg)

    return LineString(coords)


def build_edge_gdf(G):
    records = []
    for u, v, data in G.edges(data=True):
        geom = data.get("geom", None)
        if geom is None or geom.is_empty:
            geom = LineString([u, v])

        records.append({
            "u": u,
            "v": v,
            "geometry": geom,
        })

    if not records:
        raise RuntimeError("그래프 엣지가 없습니다.")

    return gpd.GeoDataFrame(records, geometry="geometry", crs="EPSG:5179")


def read_barrier_layer(gpkg_path, barrier_type):
    layers = fiona.listlayers(gpkg_path)
    if not layers:
        raise RuntimeError("장애물 GPKG에서 레이어를 찾지 못했습니다.")

    if len(layers) == 1:
        layer_name = layers[0]
    else:
        print("장애물 GPKG 레이어 목록:")
        for i, lyr in enumerate(layers, start=1):
            print(f"  {i}. {lyr}")
        s = input("사용할 레이어 번호 입력(기본=1): ").strip()
        idx = int(s) if s else 1
        if idx < 1 or idx > len(layers):
            raise RuntimeError("레이어 번호가 유효하지 않습니다.")
        layer_name = layers[idx - 1]

    gdf = gpd.read_file(gpkg_path, layer=layer_name)
    if gdf.empty:
        raise RuntimeError("장애물 레이어가 비어 있습니다.")
    if gdf.crs is None:
        raise RuntimeError("장애물 레이어 CRS가 없습니다.")

    gdf = gdf[gdf.geometry.notna()].copy()
    gdf = gdf[~gdf.geometry.is_empty].copy()
    if gdf.empty:
        raise RuntimeError("유효한 장애물 geometry가 없습니다.")

    gdf = gdf.to_crs(5179)

    geom_types = set(gdf.geom_type.unique().tolist())
    if barrier_type == "line":
        allowed = {"LineString", "MultiLineString"}
        if not geom_types.issubset(allowed):
            raise RuntimeError(f"라인 장애물에는 선 geometry만 허용됩니다: {sorted(geom_types)}")
    elif barrier_type == "polygon":
        allowed = {"Polygon", "MultiPolygon"}
        if not geom_types.issubset(allowed):
            raise RuntimeError(f"폴리곤 장애물에는 면 geometry만 허용됩니다: {sorted(geom_types)}")
    else:
        raise RuntimeError("barrier_type은 line 또는 polygon만 가능합니다.")

    return gdf


def _split_edge_at_chain_idx(G, u, v, chain_idx):
    """축약 간선 (u,v)를 chain_nodes[chain_idx]에서 분할한다.
    chain_costs 합산으로 비용을 정확히 분할하므로 보간 오차가 없다.
    반환 구조는 _split_edge와 동일하여 restore_splits로 복원 가능하다."""
    cost_keys = ("length_3dkm", "hour_ks", "hour_tob", "kcal_ks", "kcal_tob")

    fwd = dict(G[u][v])
    has_rev = G.has_edge(v, u)
    rev = dict(G[v][u]) if has_rev else None

    chain = fwd["chain_nodes"]
    cc = fwd["chain_costs"]
    split_node = chain[chain_idx]

    # 정방향 비용 분할
    costs_us = {k: sum(cc[j][k] for j in range(chain_idx)) for k in cost_keys}
    costs_sv = {k: sum(cc[j][k] for j in range(chain_idx, len(cc))) for k in cost_keys}

    # 정방향 geometry 분할
    fwd_geom = fwd["geom"]
    d = fwd_geom.project(Point(split_node))
    geom_us = substring(fwd_geom, 0, d, normalized=False)
    geom_sv = substring(fwd_geom, d, fwd_geom.length, normalized=False)
    if geom_us.geom_type == "Point":
        geom_us = LineString([geom_us.coords[0], split_node])
    if geom_sv.geom_type == "Point":
        geom_sv = LineString([split_node, geom_sv.coords[-1]])

    # 정방향 분할 간선 속성
    attrs_us = {**costs_us, "geom": geom_us}
    attrs_sv = {**costs_sv, "geom": geom_sv}
    chain_us = chain[:chain_idx + 1]
    cc_us = cc[:chain_idx]
    chain_sv = chain[chain_idx:]
    cc_sv = cc[chain_idx:]
    if len(chain_us) > 2:
        attrs_us["chain_nodes"] = chain_us
        attrs_us["chain_costs"] = cc_us
    if len(chain_sv) > 2:
        attrs_sv["chain_nodes"] = chain_sv
        attrs_sv["chain_costs"] = cc_sv

    # 역방향 처리
    rev_attrs_vs = None
    rev_attrs_su = None
    if has_rev:
        rev_chain = rev["chain_nodes"]
        rev_cc = rev["chain_costs"]
        rev_idx = None
        for ri, rn in enumerate(rev_chain):
            if rn == split_node:
                rev_idx = ri
                break

        if rev_idx is not None:
            costs_vs = {k: sum(rev_cc[j][k] for j in range(rev_idx)) for k in cost_keys}
            costs_su = {k: sum(rev_cc[j][k] for j in range(rev_idx, len(rev_cc))) for k in cost_keys}

            rev_geom = rev["geom"]
            d_rev = rev_geom.project(Point(split_node))
            geom_vs = substring(rev_geom, 0, d_rev, normalized=False)
            geom_su = substring(rev_geom, d_rev, rev_geom.length, normalized=False)
            if geom_vs.geom_type == "Point":
                geom_vs = LineString([geom_vs.coords[0], split_node])
            if geom_su.geom_type == "Point":
                geom_su = LineString([split_node, geom_su.coords[-1]])

            rev_attrs_vs = {**costs_vs, "geom": geom_vs}
            rev_attrs_su = {**costs_su, "geom": geom_su}
            rev_chain_vs = rev_chain[:rev_idx + 1]
            rev_cc_vs = rev_cc[:rev_idx]
            rev_chain_su = rev_chain[rev_idx:]
            rev_cc_su = rev_cc[rev_idx:]
            if len(rev_chain_vs) > 2:
                rev_attrs_vs["chain_nodes"] = rev_chain_vs
                rev_attrs_vs["chain_costs"] = rev_cc_vs
            if len(rev_chain_su) > 2:
                rev_attrs_su["chain_nodes"] = rev_chain_su
                rev_attrs_su["chain_costs"] = rev_cc_su

    # 간선 교체
    split_rev = has_rev and rev_attrs_vs is not None
    G.remove_edge(u, v)
    if split_rev:
        G.remove_edge(v, u)

    G.add_node(split_node)
    G.add_edge(u, split_node, **attrs_us)
    G.add_edge(split_node, v, **attrs_sv)

    if split_rev:
        G.add_edge(v, split_node, **rev_attrs_vs)
        G.add_edge(split_node, u, **rev_attrs_su)

    return {
        "temp_node": split_node,
        "u": u, "v": v,
        "fwd_data": fwd,
        "has_rev": split_rev,
        "rev_data": rev if split_rev else None,
    }


def apply_barrier_edges(G, edges_gdf, barrier_gdf):
    try:
        barrier_union = barrier_gdf.geometry.union_all()
    except Exception:
        barrier_union = barrier_gdf.geometry.unary_union

    if barrier_union is None or barrier_union.is_empty:
        return G, edges_gdf, 0

    minx, miny, maxx, maxy = barrier_union.bounds
    cand_idx = list(edges_gdf.sindex.intersection((minx, miny, maxx, maxy)))
    if not cand_idx:
        return G, edges_gdf, 0

    sub = edges_gdf.iloc[cand_idx]

    # 1단계: 교차 간선 분류
    remove_whole = []       # 통째 제거 대상
    split_targets = []      # (u, v, chain_idx, remove_side) 분할 후 부분 제거 대상
    processed = set()

    def _barrier_inside(node_xy):
        return barrier_union.covers(Point(node_xy[0], node_xy[1]))

    for row in sub.itertuples(index=False):
        geom = row.geometry
        if geom is None or geom.is_empty:
            continue
        if not geom.intersects(barrier_union):
            continue

        u, v = row.u, row.v
        key = (u, v)
        if key in processed:
            continue
        processed.add(key)

        edge_data = G[u][v]
        chain = edge_data.get("chain_nodes")

        if chain is None or len(chain) <= 2:
            # 비축약 간선 (~30m): 통째 제거
            remove_whole.append(key)
            continue

        u_in = _barrier_inside(u)
        v_in = _barrier_inside(v)

        if u_in and v_in:
            remove_whole.append(key)
        elif (not u_in) and v_in:
            # u 외부 → v 내부: 첫 진입 chain_node에서 분할, v쪽 제거
            first_inside_idx = None
            for ci in range(1, len(chain)):
                if _barrier_inside(chain[ci]):
                    first_inside_idx = ci
                    break
            if first_inside_idx is None or first_inside_idx == len(chain) - 1:
                remove_whole.append(key)
            else:
                split_targets.append((u, v, first_inside_idx, "after"))
        elif u_in and (not v_in):
            # u 내부 → v 외부: 첫 이탈 chain_node에서 분할, u쪽 제거
            first_outside_idx = None
            for ci in range(1, len(chain)):
                if not _barrier_inside(chain[ci]):
                    first_outside_idx = ci
                    break
            if first_outside_idx is None or first_outside_idx == len(chain) - 1:
                remove_whole.append(key)
            else:
                split_targets.append((u, v, first_outside_idx, "before"))
        else:
            # u, v 모두 외부이나 중간에 장애물 통과
            remove_whole.append(key)

    if not remove_whole and not split_targets:
        return G, edges_gdf, 0

    G_blocked = G.copy()

    # 2단계: 축약 간선 분할 후 장애물 내부 구간만 제거
    split_remove = []
    for u, v, chain_idx, remove_side in split_targets:
        if not G_blocked.has_edge(u, v):
            continue
        restore = _split_edge_at_chain_idx(G_blocked, u, v, chain_idx)
        split_node = restore["temp_node"]
        if remove_side == "after":
            # split_node → v 구간 제거
            if G_blocked.has_edge(split_node, v):
                split_remove.append((split_node, v))
            if G_blocked.has_edge(v, split_node):
                split_remove.append((v, split_node))
        else:
            # u → split_node 구간 제거
            if G_blocked.has_edge(u, split_node):
                split_remove.append((u, split_node))
            if G_blocked.has_edge(split_node, u):
                split_remove.append((split_node, u))

    # 3단계: 간선 제거
    G_blocked.remove_edges_from(remove_whole)
    G_blocked.remove_edges_from(split_remove)

    # 4단계: 분할로 생성된 간선을 반영하여 edge_gdf 재생성
    blocked_edges_gdf = build_edge_gdf(G_blocked)

    removed_total = len(remove_whole) + len(split_remove)
    print(f"    [barrier] whole_remove={len(remove_whole)}, split_partial={len(split_targets)}, split_remove={len(split_remove)}")
    return G_blocked, blocked_edges_gdf, removed_total


def shortest_path_leg(G, src_node, dst_node, weight_attr):
    nodes = nx.shortest_path(G, source=src_node, target=dst_node, weight=weight_attr)
    costs = accumulate_costs(G, nodes)
    metric = path_metric_from_costs(costs, weight_attr)
    return {
        "nodes": nodes,
        "costs": costs,
        "metric": metric,
    }


def collect_leg_candidates_with_cap(G, src_node, dst_node, weight_attr, metric_cap):
    """T_j 상한 기반 후보 생성. metric_cap 초과 후보 1개 포함 후 종료."""
    candidates = []
    gen = nx.shortest_simple_paths(G, source=src_node, target=dst_node, weight=weight_attr)
    t0 = time.time()

    for path_nodes in gen:
        costs = accumulate_costs(G, path_nodes)
        metric = path_metric_from_costs(costs, weight_attr)

        candidates.append({
            "nodes": path_nodes,
            "costs": costs,
            "metric": metric,
        })

        k = len(candidates)
        if k % 1000 == 0:
            elapsed = time.time() - t0
            print(f"    [cap_collect] k={k}, metric={metric:.4f}, cap={metric_cap:.4f}, elapsed={elapsed:.1f}s")

        if metric > metric_cap:
            break
        if k >= K_MAX:
            break

    elapsed = time.time() - t0
    print(f"    [cap_collect] 완료: k={len(candidates)}, cap={metric_cap:.4f}, elapsed={elapsed:.1f}s")
    return candidates if candidates else None


def merge_leg_paths(leg_paths):
    merged = []
    for leg in leg_paths:
        nodes = leg["nodes"]
        if not merged:
            merged.extend(nodes)
        else:
            merged.extend(nodes[1:])
    return merged


def sum_cost_dicts(cost_list):
    out = {
        "length_km": 0.0,
        "hour_ks": 0.0,
        "hour_tob": 0.0,
        "kcal_ks": 0.0,
        "kcal_tob": 0.0,
    }
    for c in cost_list:
        out["length_km"] += float(c["length_km"])
        out["hour_ks"] += float(c["hour_ks"])
        out["hour_tob"] += float(c["hour_tob"])
        out["kcal_ks"] += float(c["kcal_ks"])
        out["kcal_tob"] += float(c["kcal_tob"])
    return out


def path_edge_set(G, path_nodes):
    """노드열에서 간선 집합과 간선별 길이(length_3dkm)를 반환한다.
    정방향·역방향을 동일 간선으로 보기 위해 frozenset 키를 사용한다."""
    edges = {}
    for a, b in zip(path_nodes[:-1], path_nodes[1:]):
        edges[frozenset((a, b))] = float(G[a][b]["length_3dkm"])
    return edges


def weighted_jaccard(es_a, es_b):
    """길이 가중 Jaccard. es_*: {edge_key: length_3dkm}."""
    keys_a = set(es_a.keys())
    keys_b = set(es_b.keys())
    if not keys_a and not keys_b:
        return 1.0
    inter = keys_a & keys_b
    inter_len = sum(es_a[k] for k in inter)
    union_len = sum(es_a[k] for k in keys_a) + sum(es_b[k] for k in keys_b - keys_a)
    if union_len <= 0:
        return 0.0
    return inter_len / union_len


def _passes_diversity(es, selected_es):
    """이미 선택된 경로들과 모두 Jaccard < SIM_THRESHOLD이면 True."""
    for prev in selected_es:
        if weighted_jaccard(es, prev) >= SIM_THRESHOLD:
            return False
    return True


def _assemble_combo(chosen, candidate_n):
    """leg 후보 리스트를 합쳐 단일 결과 dict로 만든다."""
    merged_nodes = merge_leg_paths(chosen)
    total_costs = sum_cost_dicts([x["costs"] for x in chosen])
    total_metric = sum(float(x["metric"]) for x in chosen)
    return {
        "nodes": merged_nodes,
        "costs": total_costs,
        "metric": float(total_metric),
        "candidate_n": candidate_n,
        "legs_n": len(chosen),
    }


def collect_interval_single_leg(G, src, dst, weight_attr, beg, end, top_n):
    """단일 구간: K-최단 후보를 metric 오름차순 열거.
    [beg, end] 범위 내 후보를 다양성 필터를 거쳐 수집한다.
    metric 오름차순이므로 metric > end 시 즉시 종료.
    구간 내 후보가 없으면 최단 1개(fallback)를 반환한다.

    Returns: (results_list, is_fallback)
    """
    gen = nx.shortest_simple_paths(G, source=src, target=dst, weight=weight_attr)
    selected = []
    selected_es = []
    first_cand = None
    k = 0
    t0 = time.time()

    for path_nodes in gen:
        costs = accumulate_costs(G, path_nodes)
        metric = path_metric_from_costs(costs, weight_attr)
        k += 1

        if first_cand is None:
            first_cand = {"nodes": path_nodes, "costs": costs, "metric": metric}

        if k % 1000 == 0:
            print(f"    [interval] k={k}, metric={metric:.4f}, selected={len(selected)}, elapsed={time.time()-t0:.1f}s")

        if metric > end:
            break

        if metric >= beg:
            es = path_edge_set(G, path_nodes)
            if _passes_diversity(es, selected_es):
                selected.append({"nodes": path_nodes, "costs": costs, "metric": metric})
                selected_es.append(es)
                if len(selected) >= top_n:
                    break

        if k >= K_MAX:
            break

    print(f"    [interval] 완료: k={k}, selected={len(selected)}, elapsed={time.time()-t0:.1f}s")

    if not selected:
        first_cand["candidate_n"] = k
        first_cand["legs_n"] = 1
        return [first_cand], True

    for r in selected:
        r["candidate_n"] = k
        r["legs_n"] = 1
    return selected, False


def collect_interval_two_legs(G, legs, weight_attr, beg, end, top_n, t_j_caps):
    """2-leg 구간: T_j 작은 leg bisect(사전 생성) + 큰 leg lazy iterate.
    [beg, end] 범위 내 조합을 수집 후 metric 오름차순 정렬 → 다양성 필터.
    구간 내 조합이 없으면 최단 조합 1개(fallback)를 반환한다.
    """
    # probe 진단: 각 leg을 K_PROBE개만 짧게 열거하여 metric 증가 기울기를 측정한다.
    # 기울기(후보 100개당 metric 증가량)가 가파를수록 상한까지 후보가 적어 bisect(사전 생성)에 유리하다.
    def _probe_slope(src, dst):
        g = nx.shortest_simple_paths(G, source=src, target=dst, weight=weight_attr)
        first_m = None
        last_m = None
        cnt = 0
        for pn in g:
            m = path_metric_from_costs(accumulate_costs(G, pn), weight_attr)
            if first_m is None:
                first_m = m
            last_m = m
            cnt += 1
            if cnt >= K_PROBE:
                break
        if cnt <= 1:
            return float("inf"), first_m, cnt  # 후보가 1개뿐이면 매우 희소 → bisect 적격
        slope = (last_m - first_m) / cnt * 100.0  # 100개당 증가량
        return slope, first_m, cnt

    slope0, fm0, n0 = _probe_slope(legs[0][0], legs[0][1])
    slope1, fm1, n1 = _probe_slope(legs[1][0], legs[1][1])
    print(f"  [probe] leg0: slope={slope0:.4f}/100, first={fm0:.4f}, n={n0}")
    print(f"  [probe] leg1: slope={slope1:.4f}/100, first={fm1:.4f}, n={n1}")

    # 역할 결정: 기울기가 더 가파른(후보가 빨리 소진되는) leg을 bisect로.
    # 두 기울기 비율이 PROBE_RATIO 이내로 비슷하면 기존 T_j 규칙으로 폴백.
    lo_s, hi_s = (slope0, slope1) if slope0 <= slope1 else (slope1, slope0)
    if lo_s > 0 and hi_s / lo_s >= PROBE_RATIO:
        bi = 0 if slope0 >= slope1 else 1
        li = 1 - bi
        print(f"  [probe] 기울기 기반 역할 결정: bisect_leg={bi}(가파름)")
    else:
        if t_j_caps[0] <= t_j_caps[1]:
            bi, li = 0, 1
        else:
            bi, li = 1, 0
        print(f"  [probe] 기울기 유사 → T_j 규칙 폴백: bisect_leg={bi}")

    print(f"  [interval] 2-leg: bisect_leg={bi}(T_j={t_j_caps[bi]:.4f}), iterate_leg={li}(T_j={t_j_caps[li]:.4f})")

    bisect_cands = collect_leg_candidates_with_cap(
        G, legs[bi][0], legs[bi][1], weight_attr, t_j_caps[bi],
    )
    if bisect_cands is None:
        return None, False

    bm = [float(c["metric"]) for c in bisect_cands]
    N_b = len(bm)
    min_b = bm[0]

    gen = nx.shortest_simple_paths(
        G, source=legs[li][0], target=legs[li][1], weight=weight_attr,
    )

    # [beg, end] 범위 조합을 total 오름차순으로 방출하는 min-heap 병합(k-smallest-pairs).
    # iterate 행(gen, ci 오름차순)을 지연 pull, bisect 배열 bm을 col로 둔다. 각 행은 lo_i
    # (범위 하한 통과 첫 col)에서 seed한다. 행의 in-range 최소 total = ci + bm[lo_i]는 ci에
    # 단조가 아니므로(하한 때문), 활성화는 "다음 행 하한(ci + min_b)이 현재 heap 최소 total
    # 이하인 동안 미리 pull"하는 lower-bound 방식으로 보장한다. 전역 total 오름차순이 보장되며
    # heap·iter_cands 크기는 활성 행 수에 비례한다(조합 전수가 아님).
    min_b = bm[0]
    iter_cands = []          # 지연 pull한 iterate 후보 (인덱스로 참조)
    _buf = []                # 다음 행 미리보기 버퍼 (길이 0 또는 1)
    gen_done = False
    fb_cand = None           # fallback: 최소 total 조합 (iter_cand, j)
    fb_total = None
    merge_heap = []          # (total, fi, j)
    t0 = time.time()

    def _peek_cf():
        """다음 행 ci를 미리본다(소비하지 않음). 소진 시 None."""
        nonlocal gen_done
        if _buf:
            return _buf[0]["metric"]
        if gen_done:
            return None
        try:
            pn = next(gen)
        except StopIteration:
            gen_done = True
            return None
        cst = accumulate_costs(G, pn)
        m = path_metric_from_costs(cst, weight_attr)
        _buf.append({"nodes": pn, "costs": cst, "metric": m})
        return m

    def _take_and_seed():
        """미리본 행을 소비하여 in-range 첫 col을 merge_heap에 seed한다."""
        nonlocal fb_cand, fb_total
        cand = _buf.pop(0)
        cf = cand["metric"]
        fi = len(iter_cands)
        iter_cands.append(cand)
        if fb_cand is None or (cf + min_b) < fb_total:
            fb_total = cf + min_b
            fb_cand = (cand, 0)
        lo = _bisect.bisect_left(bm, beg - cf)
        if lo < N_b and (cf + bm[lo]) <= end:
            heapq.heappush(merge_heap, (cf + bm[lo], fi, lo))

    selected = []
    selected_es = []
    popped = 0

    while True:
        # lower-bound pull-ahead: 다음 행 하한이 현재 heap 최소 이하인 동안 미리 pull
        while True:
            nc = _peek_cf()
            if nc is None:
                break
            if (nc + min_b) > end or nc > t_j_caps[li] or len(iter_cands) >= K_MAX:
                gen_done = True
                _buf.clear()
                break
            if merge_heap and (nc + min_b) > merge_heap[0][0]:
                break
            _take_and_seed()

        if not merge_heap:
            break

        total, fi, j = heapq.heappop(merge_heap)
        if total > end:
            break
        popped += 1

        if popped % 100000 == 0:
            print(
                f"    [iterate] popped={popped}, total={total:.4f}, heap={len(merge_heap)}, pulled={len(iter_cands)}, selected={len(selected)}, elapsed={time.time() - t0:.1f}s")

        cf = iter_cands[fi]["metric"]
        if (j + 1) < N_b and (cf + bm[j + 1]) <= end:
            heapq.heappush(merge_heap, (cf + bm[j + 1], fi, j + 1))

        if total >= beg:
            chosen = [None, None]
            chosen[li] = iter_cands[fi]
            chosen[bi] = bisect_cands[j]
            merged = merge_leg_paths(chosen)
            es = path_edge_set(G, merged)
            if _passes_diversity(es, selected_es):
                selected.append(_assemble_combo(chosen, 0))  # candidate_n은 종료 후 보정
                selected_es.append(es)
                if len(selected) >= top_n:
                    break

    candidate_n = len(iter_cands) + len(bisect_cands)
    print(
        f"    [iterate] 완료: pulled={len(iter_cands)}, popped={popped}, selected={len(selected)}, elapsed={time.time() - t0:.1f}s")

    if not selected:
        if fb_cand is None:
            return None, False
        iter_cand, j = fb_cand
        chosen = [None, None]
        chosen[li] = iter_cand
        chosen[bi] = bisect_cands[j]
        return [_assemble_combo(chosen, candidate_n)], True

    for s in selected:
        s["candidate_n"] = candidate_n

    return selected, False


def collect_interval_three_plus(G, legs, weight_attr, beg, end, top_n, t_j_caps):
    """3-leg 이상 구간: T_j 사전 생성 + 앞쪽 legs min-heap + 마지막 leg bisect.
    [beg, end] 범위 내 조합을 수집 후 metric 오름차순 정렬 → 다양성 필터.
    구간 내 조합이 없으면 최단 조합 1개(fallback)를 반환한다.
    """
    L = len(legs)
    all_leg_candidates = []
    candidate_n = 0
    for idx, (src, dst) in enumerate(legs):
        print(f"    [3+leg] leg {idx+1}/{L} 후보 생성 시작")
        cands = collect_leg_candidates_with_cap(G, src, dst, weight_attr, t_j_caps[idx])
        if cands is None:
            return None, False
        all_leg_candidates.append(cands)
        candidate_n += len(cands)

    front_cands = all_leg_candidates[:-1]
    last_cands = all_leg_candidates[-1]
    F = len(front_cands)
    front_metrics = [[float(c["metric"]) for c in cands] for cands in front_cands]
    front_sizes = [len(m) for m in front_metrics]
    last_metrics = [float(c["metric"]) for c in last_cands]
    N_last = len(last_metrics)
    min_last = last_metrics[0]

    start = tuple([0] * F)
    start_sum = sum(front_metrics[j][0] for j in range(F))
    combo_heap = [(start_sum, start)]  # 앞쪽 조합 min-heap (c_front 오름차순)
    visited = {start}

    # [beg, end] 범위 조합을 total 오름차순으로 방출하는 min-heap 병합(k-smallest-pairs).
    # 앞쪽 조합(combo_heap, c_front 오름차순)을 지연 pull, 마지막 leg(last_metrics)을 col로 둔다.
    # 각 행은 lo_i(범위 하한 통과 첫 col)에서 seed하며, 하한 때문에 행의 in-range 최소 total이
    # c_front에 단조가 아니므로, "다음 행 하한(c_front + min_last)이 현재 heap 최소 total 이하인
    # 동안 미리 pull"하는 lower-bound 방식으로 전역 total 오름차순을 보장한다.
    min_last = last_metrics[0]
    front_pool = []          # [(c_front, state), ...]
    fb = None                # fallback: (state, last_idx=0)
    fb_total = None
    merge_heap = []          # (total, fi, j)
    t0 = time.time()

    def _peek_cf():
        """다음 앞쪽 조합 c_front를 미리본다(소비하지 않음). 소진 시 None."""
        return combo_heap[0][0] if combo_heap else None

    def _take_and_seed():
        """미리본 앞쪽 조합을 소비하여 후속 상태 확장 + in-range 첫 col seed."""
        nonlocal fb, fb_total
        c_front, state = heapq.heappop(combo_heap)
        fi = len(front_pool)
        front_pool.append((c_front, state))
        if fb is None or (c_front + min_last) < fb_total:
            fb_total = c_front + min_last
            fb = (state, 0)
        for jj in range(F):
            new_idx = state[jj] + 1
            if new_idx < front_sizes[jj]:
                ns = list(state)
                ns[jj] = new_idx
                ns = tuple(ns)
                if ns not in visited:
                    visited.add(ns)
                    new_sum = c_front - front_metrics[jj][state[jj]] + front_metrics[jj][new_idx]
                    heapq.heappush(combo_heap, (new_sum, ns))
        lo = _bisect.bisect_left(last_metrics, beg - c_front)
        if lo < N_last and (c_front + last_metrics[lo]) <= end:
            heapq.heappush(merge_heap, (c_front + last_metrics[lo], fi, lo))

    selected = []
    selected_es = []
    popped = 0

    while True:
        while True:
            nc = _peek_cf()
            if nc is None:
                break
            if (nc + min_last) > end:
                break
            if merge_heap and (nc + min_last) > merge_heap[0][0]:
                break
            _take_and_seed()

        if not merge_heap:
            break

        total, fi, j = heapq.heappop(merge_heap)
        if total > end:
            break
        popped += 1

        if popped % 10000 == 0:
            print(
                f"    [heap] popped={popped}, total={total:.4f}, heap={len(merge_heap)}, fronts={len(front_pool)}, selected={len(selected)}, elapsed={time.time() - t0:.1f}s")

        c_front, state = front_pool[fi]
        if (j + 1) < N_last and (c_front + last_metrics[j + 1]) <= end:
            heapq.heappush(merge_heap, (c_front + last_metrics[j + 1], fi, j + 1))

        if total >= beg:
            chosen = [front_cands[k][state[k]] for k in range(F)]
            chosen.append(last_cands[j])
            merged = merge_leg_paths(chosen)
            es = path_edge_set(G, merged)
            if _passes_diversity(es, selected_es):
                selected.append(_assemble_combo(chosen, candidate_n))
                selected_es.append(es)
                if len(selected) >= top_n:
                    break

    print(
        f"    [3+leg] 조합 탐색 완료: popped={popped}, fronts={len(front_pool)}, selected={len(selected)}, {time.time() - t0:.1f}s")

    if not selected:
        if fb is None:
            return None, False
        state, last_idx = fb
        chosen = [front_cands[k][state[k]] for k in range(F)]
        chosen.append(last_cands[last_idx])
        return [_assemble_combo(chosen, candidate_n)], True

    return selected, False


def choose_target_route_interval(G, snapped_nodes, weight_attr, beg, end, top_n=TOP_N_INTERVAL, exclude_fallback=False):
    """구간 [beg, end] 내 경로를 metric 오름차순으로 수집한다.
    다양성 필터(Jaccard >= SIM_THRESHOLD 제외) 통과분에 한해 상위 top_n개를 반환한다.
    구간 내 후보가 없으면 구간에 가장 가까운 1개(최단)를 반환한다.
    exclude_fallback=True(편차범위)이면 구간 내 후보가 없을 때 fallback 없이 None을 반환한다.

    Returns: list of result dict (rank 순서) 또는 None
    """
    legs = list(zip(snapped_nodes[:-1], snapped_nodes[1:]))
    if not legs:
        return None

    L = len(legs)
    print(f"  [interval] legs={L}, weight={weight_attr}, [{beg:.4f}, {end:.4f}], top_n={top_n}")

    if L == 1:
        results, is_fb = collect_interval_single_leg(
            G, legs[0][0], legs[0][1], weight_attr, beg, end, top_n,
        )
    else:
        min_metrics = []
        for src_node, dst_node in legs:
            leg_min = shortest_path_leg(G, src_node, dst_node, weight_attr)
            min_metrics.append(float(leg_min["metric"]))
        s_min = sum(min_metrics)
        t_j_caps = [float(end) - (s_min - min_metrics[j]) for j in range(L)]
        print(f"  [interval] S_min={s_min:.4f}, T_j_caps(end기준)={[f'{c:.4f}' for c in t_j_caps]}")

        if L == 2:
            results, is_fb = collect_interval_two_legs(
                G, legs, weight_attr, beg, end, top_n, t_j_caps,
            )
        else:
            results, is_fb = collect_interval_three_plus(
                G, legs, weight_attr, beg, end, top_n, t_j_caps,
            )

    if not results:
        return None
    if is_fb:
        if exclude_fallback:
            print(f"  [margin] 범위 내 후보 없음 → 제외")
            return None
        print(f"  [interval] 구간 내 후보 없음 → fallback 1개(최단) 출력")
    return results


def build_record(route_id, weight_attr, target_weight_attr, beg, end, mode, result, snap_start, snap_end, G, candidate_n=None, rank=None, margin=None):
    geom = path_to_linestring(G, result["nodes"])
    costs = result["costs"]

    target_metric = path_metric_from_costs(costs, target_weight_attr)

    # abs_diff: 미입력(beg None) → None / 단일(end None) → |metric-beg| / 구간 → 구간거리
    if beg is None:
        abs_diff = None
    elif end is None:
        abs_diff = float(abs(target_metric - beg))
    else:
        abs_diff = float(max(0.0, beg - target_metric, target_metric - end))

    rec = {
        "route_id": route_id,
        "weight_attr": weight_attr,
        "metric_val": float(result["metric"]),
        "abs_diff": abs_diff,
        "length_km": float(costs["length_km"]),
        "hour_ks": float(costs["hour_ks"]),
        "hour_tob": float(costs["hour_tob"]),
        "kcal_ks": float(costs["kcal_ks"]),
        "kcal_tob": float(costs["kcal_tob"]),
        "legs_n": int(result["legs_n"]),
        "snap_m_start": float(snap_start),
        "snap_m_end": float(snap_end),
        "geometry": geom,
    }

    if rank is not None:
        rec["rank"] = int(rank)
    if candidate_n is not None:
        rec["candidate_n"] = int(candidate_n)

    if mode == "km":
        rec["km_beg"] = float(beg) if beg is not None else None
        rec["km_end"] = float(end) if end is not None else None
    else:
        rec["hr_beg"] = float(beg) if beg is not None else None
        rec["hr_end"] = float(end) if end is not None else None

    if margin is not None:
        rec["margin"] = float(margin)

    return rec


def main():
    try:
        root = tk.Tk()
        root.withdraw()

        cache_path = filedialog.askopenfilename(
            title="graph_cache.pkl 선택",
            filetypes=[("Pickle", "*.pkl")],
        )
        if not cache_path:
            raise RuntimeError("graph_cache.pkl을 선택하지 않았습니다.")

        routes_csv = filedialog.askopenfilename(
            title="경로 분석.csv 선택",
            filetypes=[("CSV", "*.csv")],
        )
        if not routes_csv:
            raise RuntimeError("경로 분석.csv를 선택하지 않았습니다.")

        out_dir = filedialog.askdirectory(title="출력 폴더 선택")
        if not out_dir:
            raise RuntimeError("출력 폴더를 선택하지 않았습니다.")

        print("route_target 기준 선택:")
        print("  1) 거리 기반 - length_3dkm / target_km")
        print("  2) 시간 기반 - hour_ks / target_hr")
        print("  3) 시간 기반 - hour_tob / target_hr")
        sel = input("입력(기본=1): ").strip()

        if sel == "2":
            mode = "hr"
            target_weight_attr = "hour_ks"
        elif sel == "3":
            mode = "hr"
            target_weight_attr = "hour_tob"
        else:
            mode = "km"
            target_weight_attr = "length_3dkm"

        # 편차범위 margin: 거리=2.25 km 고정, 시간=사용자 입력
        if mode == "km":
            margin = MARGIN_KM
            print(f"[INFO] 거리 기반 margin = ±{margin} km (고정)")
        else:
            ms = input(f"시간 margin 입력(±hour, 기본={MARGIN_HR}): ").strip()
            margin = float(ms) if ms else MARGIN_HR
            print(f"[INFO] 시간 기반 margin = ±{margin} hour")

        print("장애물 사용 여부:")
        print("  1) 사용 안 함")
        print("  2) line 장애물 사용")
        print("  3) polygon 장애물 사용")
        barrier_sel = input("입력(기본=1): ").strip()

        barrier_type = None
        barrier_path = ""
        if barrier_sel == "2":
            barrier_type = "line"
        elif barrier_sel == "3":
            barrier_type = "polygon"

        if barrier_type is not None:
            barrier_path = filedialog.askopenfilename(
                title=f"{barrier_type} 장애물 GPKG 선택",
                filetypes=[("GeoPackage", "*.gpkg")],
            )
            if not barrier_path:
                raise RuntimeError("장애물 GPKG를 선택하지 않았습니다.")

        out_gpkg = os.path.join(out_dir, "경로_point_point_편차범위.gpkg")

        with open(cache_path, "rb") as f:
            payload = pickle.load(f)

        if isinstance(payload, dict) and "graph" in payload:
            G = payload["graph"]
            crs_epsg = payload.get("crs_epsg", 5179)
        else:
            G = payload
            crs_epsg = 5179

        if crs_epsg != 5179:
            raise RuntimeError("그래프 캐시 CRS가 EPSG:5179가 아닙니다.")

        routes = parse_routes_csv(routes_csv, mode)
        edges_gdf = build_edge_gdf(G)

        if barrier_type is not None:
            barrier_gdf = read_barrier_layer(barrier_path, barrier_type)
            G_work, snap_edge_gdf, removed_n = apply_barrier_edges(G, edges_gdf, barrier_gdf)
            print(f"[INFO] barrier applied: type={barrier_type}, removed_edges={removed_n}")
        else:
            G_work = G
            snap_edge_gdf = edges_gdf
            removed_n = 0
            print("[INFO] barrier not used")

        print(f"[INFO] graph: nodes={G_work.number_of_nodes():,}, edges={G_work.number_of_edges():,}")

        tf = Transformer.from_crs(4326, 5179, always_xy=True)

        feats_target = []

        total = len(routes)
        for i, row in enumerate(routes, start=1):
            route_id = row["route_id"]
            beg = row["beg"]   # None 가능 (미입력)
            end = row["end"]   # None 가능 (단일 모드)
            coords_wgs84 = row["coords"]

            t_route_start = time.time()
            if beg is None:
                target_display = "미입력"
            elif end is None:
                target_display = f"단일 T={beg}"
            else:
                target_display = f"구간 [{beg}, {end}]"
            print(f"\n[{i}/{total}] {route_id} 시작 ({target_display}, coords={len(coords_wgs84)}개)")

            snapped_nodes = []
            snap_dists = []
            restore_infos = []
            for x, y in coords_wgs84:
                x5179, y5179 = tf.transform(float(x), float(y))
                node, snap_dist, restore = snap_to_graph(x5179, y5179, snap_edge_gdf, G_work)
                if node is None:
                    snapped_nodes = []
                    break
                snapped_nodes.append(node)
                snap_dists.append(float(snap_dist))
                if restore is not None:
                    restore_infos.append(restore)

            if len(snapped_nodes) < 2:
                if restore_infos:
                    restore_splits(G_work, restore_infos)
                print(f"[{i}/{total}] {route_id} NO_SNAP")
                continue

            print(f"  snap_m: {[f'{d:.1f}' for d in snap_dists]}")

            target_results = None  # route_target 결과 리스트 (rank 순서)
            if beg is None:
                if restore_infos:
                    restore_splits(G_work, restore_infos)
                print(f"[{i}/{total}] {route_id} SKIP(목표값 미입력)")
                continue
            else:
                # 편차범위: 단일은 [beg-margin, beg+margin], 구간은 [beg-margin, end+margin]
                lo = beg - margin
                hi = (beg if end is None else end) + margin
                try:
                    target_results = choose_target_route_interval(
                        G=G_work,
                        snapped_nodes=snapped_nodes,
                        weight_attr=target_weight_attr,
                        beg=lo,
                        end=hi,
                        exclude_fallback=True,
                    )
                except (nx.NetworkXNoPath, nx.NodeNotFound):
                    target_results = None

                if not target_results:
                    if restore_infos:
                        restore_splits(G_work, restore_infos)
                    print(f"[{i}/{total}] {route_id} NO_PATH(범위 내 후보 없음)")
                    continue

            if target_results is not None:
                for rk, tr in enumerate(target_results, start=1):
                    rec_target = build_record(
                        route_id=route_id,
                        weight_attr=target_weight_attr,
                        target_weight_attr=target_weight_attr,
                        beg=beg,
                        end=end,
                        mode=mode,
                        result=tr,
                        snap_start=snap_dists[0],
                        snap_end=snap_dists[-1],
                        G=G_work,
                        candidate_n=tr["candidate_n"],
                        rank=rk,
                        margin=margin,
                    )
                    feats_target.append(rec_target)

            # 임시 노드 복원
            if restore_infos:
                restore_splits(G_work, restore_infos)

            elapsed_route = time.time() - t_route_start
            print(f"[{i}/{total}] {route_id} done ({elapsed_route:.1f}s)")

        total_saved = len(feats_target)
        if total_saved == 0:
            raise RuntimeError("생성된 경로가 없습니다.")

        if os.path.exists(out_gpkg):
            os.remove(out_gpkg)

        if feats_target:
            gdf_target = gpd.GeoDataFrame(feats_target, geometry="geometry", crs="EPSG:5179")
            gdf_target.to_file(out_gpkg, layer="route_target", driver="GPKG")

        print(f"saved gpkg: {out_gpkg}")

    except Exception as e:
        print(f"[ERROR] {e}")


if __name__ == "__main__":
    main()
