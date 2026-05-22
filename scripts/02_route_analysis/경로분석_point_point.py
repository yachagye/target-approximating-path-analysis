# 목표 거리/시간 기준 최적 경로 분석 (포인트-포인트 / 경유지 포함)
#
# 입력:
#  1) graph_cache.pkl  (네트워크데이터셋_DiGraph_변환_gpkg_pkl.py 출력물: payload dict)
#  2) 경로 분석.csv
#       route_id, x1, y1, x2, y2, ..., xN, yN, target_km 또는 target_hr
#     - 좌표: EPSG:4326 (경도 x, 위도 y)
#     - target 셀이 비어 있으면 해당 경로는 route_min_*만 산출(route_target 생략)
#  3) 선택 입력: 장애물 GPKG
#     - line 또는 polygon만 지원
#     - 장애물은 완전 차단(교차 간선 제거)
#
# route_target 알고리즘 (target 입력 경로 대상):
#  - 단일 구간(경유지 없음): K_INIT/K_MAX 2단계 후보 생성 후 |비용-목표값| 최소 선택
#  - 2-leg(경유지 1개): T_j 작은 leg bisect(사전 생성) + 큰 leg lazy iterate (heap 불필요)
#  - 3-leg 이상: T_j 사전 생성 후 앞쪽 legs min-heap + 마지막 leg bisect
#
# 출력:
#  - 경로_point_point.gpkg (EPSG:5179)
#     layer: route_min_km
#     layer: route_min_hr_ks
#     layer: route_min_hr_tob
#     layer: route_min_kcal_ks
#     layer: route_min_kcal_tob
#     layer: route_target  (target 입력 경로만 등록. route_min_* 레이어의 abs_diff·target_km/hr은 target 미입력 시 NULL)

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

# 노드 좌표 라운딩 (graph_cache.pkl 생성 시 동일 값 사용)
NODE_ROUND_M = 0.01

MIN_WEIGHT_CONFIGS = [
    ("length_3dkm", "route_min_km"),
    ("hour_ks", "route_min_hr_ks"),
    ("hour_tob", "route_min_hr_tob"),
    ("kcal_ks", "route_min_kcal_ks"),
    ("kcal_tob", "route_min_kcal_tob"),
]


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

            target_col = "target_km" if mode == "km" else "target_hr"
            if header[0].strip() != "route_id":
                raise RuntimeError("첫 컬럼은 route_id 이어야 합니다.")
            if header[-1].strip() != target_col:
                raise RuntimeError(f"마지막 컬럼은 {target_col} 이어야 합니다.")

            coord_cols = header[1:-1]
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
                target_str = str(row[-1]).strip()
                if route_id == "":
                    continue

                coords = []
                for i in range(1, len(row) - 1, 2):
                    xs = str(row[i]).strip()
                    ys = str(row[i + 1]).strip()
                    if xs == "" or ys == "":
                        continue
                    coords.append((float(xs), float(ys)))

                if len(coords) < 2:
                    continue

                out.append({
                    "route_id": route_id,
                    "coords": coords,
                    "target": float(target_str) if target_str != "" else None,
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


def k_shortest_leg_candidates(G, src_node, dst_node, weight_attr, target_value):
    """단일 leg용 후보 생성.

    종료 조건:
      (a) metric >= target_value인 후보가 등장하면 그 후보까지 포함하고 종료
      (b) K_INIT까지 (a)가 발생하지 않으면 K_MAX까지 확장
      (c) K_MAX 도달 시 종료
    K-shortest는 비용 오름차순이므로 (a) 이후 |metric - target_value| 단조 증가,
    따라서 (a) 시점 종료에도 |비용 - 목표값| 최소 후보 선택의 정확성이 보존된다.
    """
    candidates = []
    gen = nx.shortest_simple_paths(G, source=src_node, target=dst_node, weight=weight_attr)

    limit = K_INIT
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
            print(f"    [k_shortest] k={k}, metric={metric:.4f}, target={target_value:.4f}, elapsed={elapsed:.1f}s")

        # (a) 목표값 이상 후보 등장 시 즉시 종료 (K_INIT/K_MAX 단계 무관)
        if metric >= target_value:
            break

        # (b) K_INIT 도달 시 K_MAX로 확장, K_MAX 도달 시 종료
        if k >= limit:
            if limit >= K_MAX:
                break
            limit = K_MAX

    elapsed = time.time() - t0
    print(f"    [k_shortest] 완료: k={len(candidates)}, elapsed={elapsed:.1f}s")
    return candidates


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


def choose_target_combo_bisect_last(front_leg_candidates, last_leg_sorted, target_value):
    """3-leg 이상: 앞쪽 legs min-heap 조합 탐색 + 마지막 leg bisect 최적 매칭."""
    L = len(front_leg_candidates)
    last_metrics = [float(c["metric"]) for c in last_leg_sorted]
    N_last = len(last_metrics)

    front_metrics = [[float(c["metric"]) for c in cands] for cands in front_leg_candidates]
    front_sizes = [len(m) for m in front_metrics]

    start = tuple([0] * L)
    start_sum = sum(front_metrics[j][0] for j in range(L))

    heap = [(start_sum, start)]
    visited = {start}

    best_combo = None
    best_diff = float("inf")

    combos_processed = 0
    t0 = time.time()

    while heap:
        c_front, state = heapq.heappop(heap)
        combos_processed += 1

        if combos_processed % 10000 == 0:
            elapsed = time.time() - t0
            print(
                f"    [heap] combos={combos_processed}, heap_size={len(heap)}, c_front={c_front:.4f}, best_diff={best_diff:.4f}, elapsed={elapsed:.1f}s")

        if c_front + last_metrics[0] - target_value >= best_diff:
            break

        remainder = target_value - c_front
        pos = _bisect.bisect_left(last_metrics, remainder)

        for idx in (pos - 1, pos):
            if 0 <= idx < N_last:
                diff = abs(c_front + last_metrics[idx] - target_value)
                if diff < best_diff:
                    best_diff = diff
                    best_combo = (state, idx)

        for j in range(L):
            new_idx = state[j] + 1
            if new_idx < front_sizes[j]:
                new_state = list(state)
                new_state[j] = new_idx
                new_state = tuple(new_state)
                if new_state not in visited:
                    visited.add(new_state)
                    new_sum = c_front - front_metrics[j][state[j]] + front_metrics[j][new_idx]
                    heapq.heappush(heap, (new_sum, new_state))

        if best_diff < 1e-4:
            break

    if best_combo is None:
        return None

    front_state, last_idx = best_combo
    chosen = [front_leg_candidates[j][front_state[j]] for j in range(L)]
    chosen.append(last_leg_sorted[last_idx])
    return chosen


def choose_target_route(G, snapped_nodes, weight_attr, target_value):
    legs = list(zip(snapped_nodes[:-1], snapped_nodes[1:]))
    if not legs:
        return None

    print(f"  [target] legs={len(legs)}, weight={weight_attr}, T={target_value:.4f}")

    # 단일 구간: 기존 K_INIT/K_MAX 방식 유지
    if len(legs) == 1:
        cands = k_shortest_leg_candidates(G, legs[0][0], legs[0][1], weight_attr, target_value)
        if not cands:
            return None
        best = min(cands, key=lambda x: (abs(float(x["metric"]) - target_value), float(x["metric"])))
        return {
            "nodes": best["nodes"],
            "costs": best["costs"],
            "metric": float(best["metric"]),
            "candidate_n": len(cands),
            "legs_n": 1,
        }

    # multi-leg: bisect 최적화
    L = len(legs)

    # 1차 패스: 각 leg 최소 비용
    min_metrics = []
    for src_node, dst_node in legs:
        leg_min = shortest_path_leg(G, src_node, dst_node, weight_attr)
        min_metrics.append(float(leg_min["metric"]))

    s_min = sum(min_metrics)
    t_j_caps = [float(target_value) - (s_min - min_metrics[j]) for j in range(L)]

    print(f"  [target] S_min={s_min:.4f}, min_metrics={[f'{m:.4f}' for m in min_metrics]}")
    print(f"  [target] T_j_caps={[f'{c:.4f}' for c in t_j_caps]}")

    if L == 2:
        # 2-leg: T_j 작은 leg → bisect(사전 생성), 큰 leg → lazy iterate
        if t_j_caps[0] <= t_j_caps[1]:
            bi, li = 0, 1
        else:
            bi, li = 1, 0

        print(f"  [target] 2-leg: bisect_leg={bi}(T_j={t_j_caps[bi]:.4f}), iterate_leg={li}(T_j={t_j_caps[li]:.4f})")

        bisect_cands = collect_leg_candidates_with_cap(
            G, legs[bi][0], legs[bi][1], weight_attr, t_j_caps[bi],
        )
        if bisect_cands is None:
            return None

        bm = [float(c["metric"]) for c in bisect_cands]
        N_b = len(bm)
        min_b = bm[0]

        gen = nx.shortest_simple_paths(
            G, source=legs[li][0], target=legs[li][1], weight=weight_attr,
        )

        best_pair = None
        best_diff = float("inf")
        iter_n = 0
        t0 = time.time()

        for path_nodes in gen:
            costs = accumulate_costs(G, path_nodes)
            ci = path_metric_from_costs(costs, weight_attr)
            iter_n += 1

            if iter_n % 1000 == 0:
                elapsed = time.time() - t0
                print(f"    [iterate] n={iter_n}, ci={ci:.4f}, best_diff={best_diff:.4f}, elapsed={elapsed:.1f}s")

            if ci + min_b - target_value >= best_diff:
                break

            remainder = target_value - ci
            pos = _bisect.bisect_left(bm, remainder)
            for idx in (pos - 1, pos):
                if 0 <= idx < N_b:
                    diff = abs(ci + bm[idx] - target_value)
                    if diff < best_diff:
                        best_diff = diff
                        best_pair = (
                            {"nodes": path_nodes, "costs": costs, "metric": ci},
                            bisect_cands[idx],
                        )

            if ci > t_j_caps[li] or iter_n >= K_MAX or best_diff < 1e-4:
                break

        elapsed = time.time() - t0
        print(f"    [iterate] 완료: n={iter_n}, best_diff={best_diff:.4f}, elapsed={elapsed:.1f}s")

        if best_pair is None:
            return None

        chosen = [None, None]
        chosen[li] = best_pair[0]
        chosen[bi] = best_pair[1]
        candidate_n = iter_n + len(bisect_cands)

    else:
        # 3-leg 이상: T_j 사전 생성 + 앞쪽 legs heap + 마지막 leg bisect
        all_leg_candidates = []
        candidate_n = 0
        for idx, (src_node, dst_node) in enumerate(legs):
            print(f"    [3+leg] leg {idx+1}/{L} 후보 생성 시작")
            cands = collect_leg_candidates_with_cap(G, src_node, dst_node, weight_attr, t_j_caps[idx])
            if cands is None:
                return None
            all_leg_candidates.append(cands)
            candidate_n += len(cands)

        print(f"    [3+leg] 후보 합계: {candidate_n}, heap 조합 탐색 시작")
        t0 = time.time()
        chosen = choose_target_combo_bisect_last(
            all_leg_candidates[:-1], all_leg_candidates[-1], target_value,
        )
        elapsed = time.time() - t0
        print(f"    [3+leg] heap 조합 탐색 완료: {elapsed:.1f}s")

    if chosen is None:
        return None

    merged_nodes = merge_leg_paths(chosen)
    total_costs = sum_cost_dicts([x["costs"] for x in chosen])
    total_metric = sum(float(x["metric"]) for x in chosen)

    return {
        "nodes": merged_nodes,
        "costs": total_costs,
        "metric": float(total_metric),
        "candidate_n": candidate_n,
        "legs_n": L,
    }


def choose_min_route(G, snapped_nodes, weight_attr):
    legs = list(zip(snapped_nodes[:-1], snapped_nodes[1:]))
    if not legs:
        return None

    leg_paths = []
    for src_node, dst_node in legs:
        leg_paths.append(shortest_path_leg(G, src_node, dst_node, weight_attr))

    merged_nodes = merge_leg_paths(leg_paths)
    total_costs = sum_cost_dicts([x["costs"] for x in leg_paths])
    total_metric = sum(float(x["metric"]) for x in leg_paths)

    return {
        "nodes": merged_nodes,
        "costs": total_costs,
        "metric": float(total_metric),
        "legs_n": len(legs),
    }


def build_record(route_id, weight_attr, target_weight_attr, target_value, mode, result, snap_start, snap_end, G, candidate_n=None):
    geom = path_to_linestring(G, result["nodes"])
    costs = result["costs"]

    target_metric = path_metric_from_costs(costs, target_weight_attr)
    rec = {
        "route_id": route_id,
        "weight_attr": weight_attr,
        "metric_val": float(result["metric"]),
        "abs_diff": float(abs(target_metric - target_value)) if target_value is not None else None,
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

    if candidate_n is not None:
        rec["candidate_n"] = int(candidate_n)

    if mode == "km":
        rec["target_km"] = float(target_value) if target_value is not None else None
    else:
        rec["target_hr"] = float(target_value) if target_value is not None else None

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

        out_gpkg = os.path.join(out_dir, "경로_point_point.gpkg")

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

        feats_min_map = {layer_name: [] for _, layer_name in MIN_WEIGHT_CONFIGS}
        feats_target = []

        total = len(routes)
        for i, row in enumerate(routes, start=1):
            route_id = row["route_id"]
            target_value = row["target"]  # None 가능 (target 미입력 시)
            coords_wgs84 = row["coords"]

            t_route_start = time.time()
            target_display = f"{target_value}" if target_value is not None else "미입력"
            print(f"\n[{i}/{total}] {route_id} 시작 (target={target_display}, coords={len(coords_wgs84)}개)")

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

            if target_value is None:
                target_result = None
                print(f"  [target] target 미입력, route_min만 산출")
            else:
                try:
                    target_result = choose_target_route(
                        G=G_work,
                        snapped_nodes=snapped_nodes,
                        weight_attr=target_weight_attr,
                        target_value=target_value,
                    )
                except nx.NetworkXNoPath:
                    target_result = None
                except nx.NodeNotFound:
                    target_result = None

                if target_result is None:
                    if restore_infos:
                        restore_splits(G_work, restore_infos)
                    print(f"[{i}/{total}] {route_id} NO_TARGET_PATH")
                    continue

            for min_weight_attr, layer_name in MIN_WEIGHT_CONFIGS:
                try:
                    min_result = choose_min_route(
                        G=G_work,
                        snapped_nodes=snapped_nodes,
                        weight_attr=min_weight_attr,
                    )
                except nx.NetworkXNoPath:
                    min_result = None
                except nx.NodeNotFound:
                    min_result = None

                if min_result is None:
                    continue

                rec_min = build_record(
                    route_id=route_id,
                    weight_attr=min_weight_attr,
                    target_weight_attr=target_weight_attr,
                    target_value=target_value,
                    mode=mode,
                    result=min_result,
                    snap_start=snap_dists[0],
                    snap_end=snap_dists[-1],
                    G=G_work,
                )
                feats_min_map[layer_name].append(rec_min)

            if target_result is not None:
                rec_target = build_record(
                    route_id=route_id,
                    weight_attr=target_weight_attr,
                    target_weight_attr=target_weight_attr,
                    target_value=target_value,
                    mode=mode,
                    result=target_result,
                    snap_start=snap_dists[0],
                    snap_end=snap_dists[-1],
                    G=G_work,
                    candidate_n=target_result["candidate_n"],
                )
                feats_target.append(rec_target)

            # 임시 노드 복원
            if restore_infos:
                restore_splits(G_work, restore_infos)

            elapsed_route = time.time() - t_route_start
            print(f"[{i}/{total}] {route_id} done ({elapsed_route:.1f}s)")

        total_saved = sum(len(v) for v in feats_min_map.values()) + len(feats_target)
        if total_saved == 0:
            raise RuntimeError("생성된 경로가 없습니다.")

        if os.path.exists(out_gpkg):
            os.remove(out_gpkg)

        wrote_any = False

        for _, layer_name in MIN_WEIGHT_CONFIGS:
            feats = feats_min_map[layer_name]
            if not feats:
                continue

            gdf_min = gpd.GeoDataFrame(feats, geometry="geometry", crs="EPSG:5179")
            if not wrote_any:
                gdf_min.to_file(out_gpkg, layer=layer_name, driver="GPKG")
                wrote_any = True
            else:
                gdf_min.to_file(out_gpkg, layer=layer_name, driver="GPKG", mode="a")

        if feats_target:
            gdf_target = gpd.GeoDataFrame(feats_target, geometry="geometry", crs="EPSG:5179")
            if not wrote_any:
                gdf_target.to_file(out_gpkg, layer="route_target", driver="GPKG")
                wrote_any = True
            else:
                gdf_target.to_file(out_gpkg, layer="route_target", driver="GPKG", mode="a")

        print(f"saved gpkg: {out_gpkg}")

    except Exception as e:
        print(f"[ERROR] {e}")


if __name__ == "__main__":
    main()
