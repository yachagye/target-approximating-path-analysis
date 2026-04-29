# 목표값 ±margin 범위 내 상위 N위 경로 분석 (포인트-포인트 / 경유지 포함)
#
# 입력:
#  1) graph_cache.pkl
#  2) 경로 분석.csv (route_id, x1, y1, ..., xN, yN, target_km 또는 target_hr)
#  3) 선택: 장애물 GPKG (line / polygon)
#
# margin:
#  - 거리 기반(target_km): ±2.25 km 고정 (5리 × 0.45 km)
#  - 시간 기반(target_hr): 실행 시 사용자 입력 (hour 단위)
#
# 출력:
#  - 경로_point_point_편차범위.gpkg (EPSG:5179)
#     layer: route_target  (경로당 최대 TOP_N개, rank 속성 포함)
#     margin 범위 내 후보가 없으면 해당 경로는 출력되지 않음

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
from shapely.geometry import Point, LineString
from shapely.ops import substring
from pyproj import Transformer


TOP_N = 500
K_MAX = 30000
MARGIN_KM = 2.25  # 5리 × 0.45 km

NODE_ROUND_M = 0.01


# ── CSV 파싱 ──────────────────────────────────────────────

def parse_routes_csv(csv_path, mode):
    def _read(enc):
        out = []
        with open(csv_path, "r", encoding=enc, newline="") as f:
            reader = csv.reader(f)
            header = next(reader, None)
            if not header:
                raise RuntimeError("CSV가 비어 있습니다.")

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
                    row += [""] * (len(header) - len(row))
                elif len(row) > len(header):
                    row = row[:len(header)]

                route_id = str(row[0]).strip()
                target_str = str(row[-1]).strip()
                if route_id == "" or target_str == "":
                    continue

                coords = []
                for i in range(1, len(row) - 1, 2):
                    xs, ys = str(row[i]).strip(), str(row[i + 1]).strip()
                    if xs == "" or ys == "":
                        continue
                    coords.append((float(xs), float(ys)))

                if len(coords) < 2:
                    continue

                out.append({"route_id": route_id, "coords": coords, "target": float(target_str)})
        return out

    try:
        rows = _read("utf-8-sig")
    except UnicodeDecodeError:
        rows = _read("cp949")

    if not rows:
        raise RuntimeError("CSV에서 유효한 경로 데이터를 찾지 못했습니다.")

    seen, dup = set(), set()
    for r in rows:
        rid = r["route_id"]
        (dup if rid in seen else seen).add(rid)
    if dup:
        raise RuntimeError(f"CSV에 중복 route_id 존재: {sorted(list(dup))[:10]}")

    return rows


# ── 네트워크 유틸 ─────────────────────────────────────────

def _round_xy(x, y, r=NODE_ROUND_M):
    return (round(x / r) * r, round(y / r) * r)


def snap_to_graph(x, y, edge_gdf, G):
    pt = Point(x, y)
    dists = edge_gdf.geometry.distance(pt)

    for nearest_idx in dists.nsmallest(10).index:
        row = edge_gdf.loc[nearest_idx]
        u, v = row["u"], row["v"]
        if not G.has_edge(u, v):
            continue

        snap_dist = float(dists.loc[nearest_idx])
        edge_geom = row.geometry
        d = edge_geom.project(pt)
        proj_pt = edge_geom.interpolate(d)

        node_tol = 15.0
        dist_u = math.sqrt((u[0] - proj_pt.x) ** 2 + (u[1] - proj_pt.y) ** 2)
        dist_v = math.sqrt((v[0] - proj_pt.x) ** 2 + (v[1] - proj_pt.y) ** 2)

        if dist_u <= node_tol:
            return u, snap_dist, None
        if dist_v <= node_tol:
            return v, snap_dist, None

        temp_node = _round_xy(proj_pt.x, proj_pt.y)
        if G.has_node(temp_node):
            return temp_node, snap_dist, None

        restore = _split_edge(G, u, v, temp_node, d)
        return temp_node, snap_dist, restore

    best_node, best_d2 = None, float("inf")
    for node in G.nodes():
        d2 = (node[0] - x) ** 2 + (node[1] - y) ** 2
        if d2 < best_d2:
            best_d2 = d2
            best_node = node
    return best_node, math.sqrt(best_d2) if best_node else float("inf"), None


def _split_edge(G, u, v, temp_node, d):
    cost_keys = ("length_3dkm", "hour_ks", "hour_tob", "kcal_ks", "kcal_tob")
    fwd = dict(G[u][v])
    has_rev = G.has_edge(v, u)
    rev = dict(G[v][u]) if has_rev else None

    fwd_geom = fwd["geom"]
    chain_nodes = fwd.get("chain_nodes", [u, v])
    chain_costs = fwd.get("chain_costs", [{k: float(fwd[k]) for k in cost_keys}])
    costs_ut, costs_tv, geom_ut, geom_tv = _interpolate_split(fwd_geom, chain_nodes, chain_costs, d, cost_keys)

    if has_rev:
        rev_geom = rev["geom"]
        rev_chain = rev.get("chain_nodes", [v, u])
        rev_costs = rev.get("chain_costs", [{k: float(rev[k]) for k in cost_keys}])
        d_rev = rev_geom.project(Point(temp_node))
        costs_vt, costs_tu, geom_vt, geom_tu = _interpolate_split(rev_geom, rev_chain, rev_costs, d_rev, cost_keys)

    G.remove_edge(u, v)
    if has_rev:
        G.remove_edge(v, u)
    G.add_node(temp_node)
    G.add_edge(u, temp_node, **costs_ut, geom=geom_ut)
    G.add_edge(temp_node, v, **costs_tv, geom=geom_tv)
    if has_rev:
        G.add_edge(v, temp_node, **costs_vt, geom=geom_vt)
        G.add_edge(temp_node, u, **costs_tu, geom=geom_tu)

    return {"temp_node": temp_node, "u": u, "v": v, "fwd_data": fwd, "has_rev": has_rev, "rev_data": rev}


def _interpolate_split(edge_geom, chain_nodes, chain_costs, d, cost_keys):
    cum = [0.0]
    for i in range(1, len(chain_nodes)):
        cum.append(edge_geom.project(Point(chain_nodes[i])))

    seg_i = len(chain_costs) - 1
    for i in range(len(cum) - 1):
        if d <= cum[i + 1] + 0.01:
            seg_i = i
            break

    span = cum[seg_i + 1] - cum[seg_i]
    f = max(0.0, min(1.0, (d - cum[seg_i]) / span if span > 0 else 0.0))

    costs_before, costs_after = {}, {}
    for k in cost_keys:
        before_full = sum(chain_costs[j][k] for j in range(seg_i))
        seg_cost = chain_costs[seg_i][k]
        after_full = sum(chain_costs[j][k] for j in range(seg_i + 1, len(chain_costs)))
        costs_before[k] = before_full + f * seg_cost
        costs_after[k] = (1 - f) * seg_cost + after_full

    geom_before = substring(edge_geom, 0, d, normalized=False)
    geom_after = substring(edge_geom, d, edge_geom.length, normalized=False)
    split_pt = edge_geom.interpolate(d)
    if geom_before.geom_type == "Point":
        geom_before = LineString([geom_before.coords[0], split_pt.coords[0]])
    if geom_after.geom_type == "Point":
        geom_after = LineString([split_pt.coords[0], geom_after.coords[0]])

    return costs_before, costs_after, geom_before, geom_after


def restore_splits(G, infos):
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


def accumulate_costs(G, path_nodes):
    s = {"length_km": 0.0, "hour_ks": 0.0, "hour_tob": 0.0, "kcal_ks": 0.0, "kcal_tob": 0.0}
    for a, b in zip(path_nodes[:-1], path_nodes[1:]):
        ed = G[a][b]
        s["length_km"] += float(ed["length_3dkm"])
        s["hour_ks"] += float(ed["hour_ks"])
        s["hour_tob"] += float(ed["hour_tob"])
        s["kcal_ks"] += float(ed["kcal_ks"])
        s["kcal_tob"] += float(ed["kcal_tob"])
    return s


def path_metric(costs, weight_attr):
    return float(costs["length_km"]) if weight_attr == "length_3dkm" else float(costs[weight_attr])


def sum_cost_dicts(cost_list):
    out = {"length_km": 0.0, "hour_ks": 0.0, "hour_tob": 0.0, "kcal_ks": 0.0, "kcal_tob": 0.0}
    for c in cost_list:
        for k in out:
            out[k] += float(c[k])
    return out


def path_to_linestring(G, path_nodes):
    coords = []
    for a, b in zip(path_nodes[:-1], path_nodes[1:]):
        geom = G[a][b].get("geom", None)
        if geom is None or geom.is_empty:
            seg = [(a[0], a[1]), (b[0], b[1])]
        else:
            seg = list(geom.coords)
        if not coords:
            coords.extend(seg)
        else:
            coords.extend(seg[1:] if coords[-1] == seg[0] else seg)
    return LineString(coords)


def merge_leg_paths(leg_paths):
    merged = []
    for leg in leg_paths:
        nodes = leg["nodes"]
        if not merged:
            merged.extend(nodes)
        else:
            merged.extend(nodes[1:])
    return merged


def build_edge_gdf(G):
    records = []
    for u, v, data in G.edges(data=True):
        geom = data.get("geom", None)
        if geom is None or geom.is_empty:
            geom = LineString([u, v])
        records.append({"u": u, "v": v, "geometry": geom})
    if not records:
        raise RuntimeError("그래프 엣지가 없습니다.")
    return gpd.GeoDataFrame(records, geometry="geometry", crs="EPSG:5179")


# ── 장애물 ────────────────────────────────────────────────

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
        layer_name = layers[idx - 1]

    gdf = gpd.read_file(gpkg_path, layer=layer_name)
    if gdf.empty:
        raise RuntimeError("장애물 레이어가 비어 있습니다.")
    if gdf.crs is None:
        raise RuntimeError("장애물 레이어 CRS가 없습니다.")
    gdf = gdf[gdf.geometry.notna() & ~gdf.geometry.is_empty].copy()
    if gdf.empty:
        raise RuntimeError("유효한 장애물 geometry가 없습니다.")
    gdf = gdf.to_crs(5179)

    geom_types = set(gdf.geom_type.unique().tolist())
    allowed = {"LineString", "MultiLineString"} if barrier_type == "line" else {"Polygon", "MultiPolygon"}
    if not geom_types.issubset(allowed):
        raise RuntimeError(f"장애물 geometry 유형 불일치: {sorted(geom_types)}")
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
        return G, 0

    minx, miny, maxx, maxy = barrier_union.bounds
    cand_idx = list(edges_gdf.sindex.intersection((minx, miny, maxx, maxy)))
    if not cand_idx:
        return G, 0

    sub = edges_gdf.iloc[cand_idx]

    # 1단계: 교차 간선 분류
    remove_whole = []
    split_targets = []
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
            remove_whole.append(key)
            continue

        u_in = _barrier_inside(u)
        v_in = _barrier_inside(v)

        if u_in and v_in:
            remove_whole.append(key)
        elif (not u_in) and v_in:
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
            remove_whole.append(key)

    if not remove_whole and not split_targets:
        return G, 0

    G_blocked = G.copy()

    # 2단계: 축약 간선 분할 후 장애물 내부 구간만 제거
    split_remove = []
    for u, v, chain_idx, remove_side in split_targets:
        if not G_blocked.has_edge(u, v):
            continue
        restore = _split_edge_at_chain_idx(G_blocked, u, v, chain_idx)
        split_node = restore["temp_node"]
        if remove_side == "after":
            if G_blocked.has_edge(split_node, v):
                split_remove.append((split_node, v))
            if G_blocked.has_edge(v, split_node):
                split_remove.append((v, split_node))
        else:
            if G_blocked.has_edge(u, split_node):
                split_remove.append((u, split_node))
            if G_blocked.has_edge(split_node, u):
                split_remove.append((split_node, u))

    # 3단계: 간선 제거
    G_blocked.remove_edges_from(remove_whole)
    G_blocked.remove_edges_from(split_remove)

    removed_total = len(remove_whole) + len(split_remove)
    print(f"    [barrier] whole_remove={len(remove_whole)}, split_partial={len(split_targets)}, split_remove={len(split_remove)}")
    return G_blocked, removed_total


# ── Margin 범위 내 Top-N 경로 탐색 ───────────────────────

def collect_leg_candidates_with_cap(G, src_node, dst_node, weight_attr, metric_cap):
    """T_j 상한 기반 후보 생성. metric_cap 초과 후보 1개 포함 후 종료."""
    candidates = []
    gen = nx.shortest_simple_paths(G, source=src_node, target=dst_node, weight=weight_attr)
    t0 = time.time()

    for path_nodes in gen:
        costs = accumulate_costs(G, path_nodes)
        metric = path_metric(costs, weight_attr)

        candidates.append({"nodes": path_nodes, "costs": costs, "metric": metric})

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


def choose_target_top_n_margin(G, snapped_nodes, weight_attr, target_value, margin, top_n=TOP_N):
    """±margin 범위 내 타겟 근접 상위 top_n개 경로를 반환한다.

    Returns:
        list of dict: [{nodes, costs, metric, candidate_n, legs_n}, ...]
                      |metric - target| 오름차순 정렬, margin 범위 내만 포함
    """
    legs = list(zip(snapped_nodes[:-1], snapped_nodes[1:]))
    if not legs:
        return []

    print(f"  [target] legs={len(legs)}, weight={weight_attr}, T={target_value:.4f}, margin=±{margin:.4f}, top_n={top_n}")

    if len(legs) == 1:
        return _top_n_single_leg(G, legs[0][0], legs[0][1], weight_attr, target_value, margin, top_n)
    else:
        return _top_n_multi_leg(G, legs, weight_attr, target_value, margin, top_n)


def _top_n_single_leg(G, src, dst, weight_attr, target_value, margin, top_n):
    """단일 구간: Yen's K-최단 후보를 T+margin까지 수집, 범위 내 필터, 정렬.

    metric이 오름차순이므로 T+margin 초과 시 이후 후보는 모두 범위 밖 → 즉시 종료.
    """
    gen = nx.shortest_simple_paths(G, source=src, target=dst, weight=weight_attr)
    in_range = []
    total_generated = 0

    t_lo = target_value - margin
    t_hi = target_value + margin

    t0 = time.time()

    for path_nodes in gen:
        costs = accumulate_costs(G, path_nodes)
        metric = path_metric(costs, weight_attr)
        total_generated += 1

        if total_generated % 1000 == 0:
            elapsed = time.time() - t0
            print(f"    [single_leg] k={total_generated}, metric={metric:.4f}, in_range={len(in_range)}, t_hi={t_hi:.4f}, elapsed={elapsed:.1f}s")

        if metric > t_hi:
            break

        if metric >= t_lo:
            in_range.append({"nodes": path_nodes, "costs": costs, "metric": metric})

        if total_generated >= K_MAX:
            break

    elapsed = time.time() - t0
    print(f"    [single_leg] 완료: k={total_generated}, in_range={len(in_range)}, elapsed={elapsed:.1f}s")

    if not in_range:
        return []

    in_range.sort(key=lambda x: (abs(x["metric"] - target_value), x["metric"]))
    results = []
    for c in in_range[:top_n]:
        results.append({
            "nodes": c["nodes"],
            "costs": c["costs"],
            "metric": c["metric"],
            "candidate_n": len(in_range),
            "legs_n": 1,
        })
    return results


def _top_n_multi_leg(G, legs, weight_attr, target_value, margin, top_n):
    """경유지 있는 경우: T_j 사전 생성 + bisect 최적화 + lazy iterate, margin 범위 필터.

    T_j 상한은 (T + margin) 기준으로 완화.
    2-leg: T_j 작은 leg bisect + 큰 leg lazy iterate, margin 범위 내 top_n 수집.
    3-leg 이상: T_j 사전 생성 + 앞쪽 legs min-heap + 마지막 leg bisect, margin 범위 내 top_n 수집.
    """
    L = len(legs)
    t_lo = target_value - margin
    t_hi = target_value + margin

    # 1차 패스: 각 leg 최소 비용
    min_metrics = []
    for src, dst in legs:
        leg_min_path = nx.shortest_path(G, source=src, target=dst, weight=weight_attr)
        leg_min_costs = accumulate_costs(G, leg_min_path)
        min_metrics.append(path_metric(leg_min_costs, weight_attr))

    s_min = sum(min_metrics)
    # T_j 상한: (T + margin) 기준
    t_j_caps = [t_hi - (s_min - min_metrics[j]) for j in range(L)]

    print(f"  [target] S_min={s_min:.4f}, min_metrics={[f'{m:.4f}' for m in min_metrics]}")
    print(f"  [target] T_j_caps (T+margin 기준)={[f'{c:.4f}' for c in t_j_caps]}")

    if L == 2:
        return _top_n_two_legs_margin(G, legs, weight_attr, target_value, margin, top_n, min_metrics, t_j_caps)
    else:
        return _top_n_three_plus_legs_margin(G, legs, weight_attr, target_value, margin, top_n, t_j_caps)


def _top_n_two_legs_margin(G, legs, weight_attr, target_value, margin, top_n, min_metrics, t_j_caps):
    """2-leg margin top_n: T_j 작은 leg bisect + 큰 leg lazy iterate."""
    t_lo = target_value - margin
    t_hi = target_value + margin

    if t_j_caps[0] <= t_j_caps[1]:
        bi, li = 0, 1
    else:
        bi, li = 1, 0

    print(f"  [target] 2-leg: bisect_leg={bi}(T_j={t_j_caps[bi]:.4f}), iterate_leg={li}(T_j={t_j_caps[li]:.4f})")

    bisect_cands = collect_leg_candidates_with_cap(
        G, legs[bi][0], legs[bi][1], weight_attr, t_j_caps[bi],
    )
    if bisect_cands is None:
        return []

    bm = [float(c["metric"]) for c in bisect_cands]
    N_b = len(bm)
    min_b = bm[0]

    gen = nx.shortest_simple_paths(
        G, source=legs[li][0], target=legs[li][1], weight=weight_attr,
    )

    top_heap = []
    seq = 0
    iter_n = 0
    t0 = time.time()

    for path_nodes in gen:
        costs = accumulate_costs(G, path_nodes)
        ci = path_metric(costs, weight_attr)
        iter_n += 1

        if iter_n % 1000 == 0:
            elapsed = time.time() - t0
            worst = -top_heap[0][0] if len(top_heap) >= top_n else float("inf")
            print(f"    [iterate] n={iter_n}, ci={ci:.4f}, worst_diff={worst:.4f}, elapsed={elapsed:.1f}s")

        # 즉시 종료: 최소 합산이 범위 상한 초과
        if ci + min_b > t_hi:
            break

        # bisect 위치 기준 양방향 확장 — margin 범위 내만 수집
        remainder = target_value - ci
        pos = _bisect.bisect_left(bm, remainder)

        worst_diff = -top_heap[0][0] if len(top_heap) >= top_n else float("inf")

        # 오른쪽 확장
        for j in range(pos, N_b):
            total = ci + bm[j]
            if total > t_hi:
                break
            if total >= t_lo:
                diff = abs(total - target_value)
                if len(top_heap) < top_n:
                    heapq.heappush(top_heap, (-diff, seq, ci, j, path_nodes, costs))
                    seq += 1
                    worst_diff = -top_heap[0][0] if len(top_heap) >= top_n else float("inf")
                elif diff < -top_heap[0][0]:
                    heapq.heapreplace(top_heap, (-diff, seq, ci, j, path_nodes, costs))
                    seq += 1
                    worst_diff = -top_heap[0][0]

        # 왼쪽 확장
        for j in range(pos - 1, -1, -1):
            total = ci + bm[j]
            if total < t_lo:
                break
            if total <= t_hi:
                diff = abs(total - target_value)
                if len(top_heap) < top_n:
                    heapq.heappush(top_heap, (-diff, seq, ci, j, path_nodes, costs))
                    seq += 1
                    worst_diff = -top_heap[0][0] if len(top_heap) >= top_n else float("inf")
                elif diff < -top_heap[0][0]:
                    heapq.heapreplace(top_heap, (-diff, seq, ci, j, path_nodes, costs))
                    seq += 1
                    worst_diff = -top_heap[0][0]

        if ci > t_j_caps[li] or iter_n >= K_MAX:
            break

    elapsed = time.time() - t0
    print(f"    [iterate] 완료: n={iter_n}, top={len(top_heap)}, elapsed={elapsed:.1f}s")

    if not top_heap:
        return []

    candidate_n = iter_n + len(bisect_cands)
    ranked = sorted(top_heap, key=lambda x: (-x[0], x[2]))

    results = []
    for neg_diff, _, ci, j, pn_iter, costs_iter in ranked:
        iter_cand = {"nodes": pn_iter, "costs": costs_iter, "metric": ci}
        chosen = [None, None]
        chosen[li] = iter_cand
        chosen[bi] = bisect_cands[j]

        merged_nodes = merge_leg_paths(chosen)
        total_costs = sum_cost_dicts([x["costs"] for x in chosen])
        total_metric = sum(float(x["metric"]) for x in chosen)
        results.append({
            "nodes": merged_nodes,
            "costs": total_costs,
            "metric": float(total_metric),
            "candidate_n": candidate_n,
            "legs_n": 2,
        })

    return results


def _top_n_three_plus_legs_margin(G, legs, weight_attr, target_value, margin, top_n, t_j_caps):
    """3-leg 이상 margin top_n: T_j 사전 생성 + 앞쪽 legs min-heap + 마지막 leg bisect."""
    L = len(legs)
    t_lo = target_value - margin
    t_hi = target_value + margin

    all_leg_candidates = []
    candidate_n = 0
    for idx, (src, dst) in enumerate(legs):
        print(f"    [3+leg] leg {idx+1}/{L} 후보 생성 시작")
        cands = collect_leg_candidates_with_cap(G, src, dst, weight_attr, t_j_caps[idx])
        if cands is None:
            return []
        all_leg_candidates.append(cands)
        candidate_n += len(cands)

    print(f"    [3+leg] 후보 합계: {candidate_n}, heap 조합 탐색 시작")
    t_heap0 = time.time()

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

    combo_heap = [(start_sum, start)]
    visited = {start}

    top_heap = []
    seq = 0

    combos_processed = 0

    while combo_heap:
        c_front, state = heapq.heappop(combo_heap)
        combos_processed += 1

        if combos_processed % 10000 == 0:
            elapsed = time.time() - t_heap0
            worst_diff = -top_heap[0][0] if len(top_heap) >= top_n else float("inf")
            print(
                f"    [heap] combos={combos_processed}, heap_size={len(combo_heap)}, top={len(top_heap)}/{top_n}, c_front={c_front:.4f}, worst_diff={worst_diff:.4f}, elapsed={elapsed:.1f}s")

        # 즉시 종료: 최소 합산이 범위 상한 초과
        if c_front + min_last > t_hi:
            break

        # bisect 위치 기준 양방향 확장 — margin 범위 내만 수집
        remainder = target_value - c_front
        pos = _bisect.bisect_left(last_metrics, remainder)

        for j in range(pos, N_last):
            total = c_front + last_metrics[j]
            if total > t_hi:
                break
            if total >= t_lo:
                diff = abs(total - target_value)
                if len(top_heap) < top_n:
                    heapq.heappush(top_heap, (-diff, seq, total, state, j))
                    seq += 1
                elif diff < -top_heap[0][0]:
                    heapq.heapreplace(top_heap, (-diff, seq, total, state, j))
                    seq += 1

        for j in range(pos - 1, -1, -1):
            total = c_front + last_metrics[j]
            if total < t_lo:
                break
            if total <= t_hi:
                diff = abs(total - target_value)
                if len(top_heap) < top_n:
                    heapq.heappush(top_heap, (-diff, seq, total, state, j))
                    seq += 1
                elif diff < -top_heap[0][0]:
                    heapq.heapreplace(top_heap, (-diff, seq, total, state, j))
                    seq += 1

        # 앞쪽 legs heap 확장
        for j in range(F):
            new_idx = state[j] + 1
            if new_idx < front_sizes[j]:
                new_state = list(state)
                new_state[j] = new_idx
                new_state = tuple(new_state)
                if new_state not in visited:
                    visited.add(new_state)
                    new_sum = c_front - front_metrics[j][state[j]] + front_metrics[j][new_idx]
                    heapq.heappush(combo_heap, (new_sum, new_state))

    elapsed_heap = time.time() - t_heap0
    print(f"    [3+leg] heap 조합 탐색 완료: top={len(top_heap)}, {elapsed_heap:.1f}s")

    if not top_heap:
        return []

    ranked = sorted(top_heap, key=lambda x: (-x[0], x[2]))

    results = []
    for neg_diff, _, combo_sum, state, last_idx in ranked:
        chosen = [front_cands[j][state[j]] for j in range(F)]
        chosen.append(last_cands[last_idx])

        merged_nodes = merge_leg_paths(chosen)
        total_costs = sum_cost_dicts([x["costs"] for x in chosen])
        total_metric = sum(float(x["metric"]) for x in chosen)
        results.append({
            "nodes": merged_nodes,
            "costs": total_costs,
            "metric": float(total_metric),
            "candidate_n": candidate_n,
            "legs_n": len(chosen),
        })

    return results


# ── 레코드 생성 ──────────────────────────────────────────

def build_record(route_id, rank, weight_attr, target_value, margin, mode, result, snap_start, snap_end, G):
    geom = path_to_linestring(G, result["nodes"])
    costs = result["costs"]

    rec = {
        "route_id": route_id,
        "rank": int(rank),
        "weight_attr": weight_attr,
        "metric_val": float(result["metric"]),
        "abs_diff": float(abs(result["metric"] - target_value)),
        "margin": float(margin),
        "length_km": float(costs["length_km"]),
        "hour_ks": float(costs["hour_ks"]),
        "hour_tob": float(costs["hour_tob"]),
        "kcal_ks": float(costs["kcal_ks"]),
        "kcal_tob": float(costs["kcal_tob"]),
        "candidate_n": int(result["candidate_n"]),
        "legs_n": int(result["legs_n"]),
        "snap_m_start": float(snap_start),
        "snap_m_end": float(snap_end),
        "geometry": geom,
    }

    if mode == "km":
        rec["target_km"] = float(target_value)
    else:
        rec["target_hr"] = float(target_value)

    return rec


# ── 메인 ─────────────────────────────────────────────────

def main():
    try:
        root = tk.Tk()
        root.withdraw()

        cache_path = filedialog.askopenfilename(title="graph_cache.pkl 선택", filetypes=[("Pickle", "*.pkl")])
        if not cache_path:
            raise RuntimeError("graph_cache.pkl을 선택하지 않았습니다.")

        routes_csv = filedialog.askopenfilename(title="경로 분석.csv 선택", filetypes=[("CSV", "*.csv")])
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
            mode, target_weight_attr = "hr", "hour_ks"
        elif sel == "3":
            mode, target_weight_attr = "hr", "hour_tob"
        else:
            mode, target_weight_attr = "km", "length_3dkm"

        # margin 설정
        if mode == "km":
            margin = MARGIN_KM
            print(f"[INFO] margin = ±{margin} km (5리 × 0.45 km 고정)")
        else:
            margin_str = input("margin 입력 (±hour, 예: 0.5): ").strip()
            margin = float(margin_str)
            print(f"[INFO] margin = ±{margin} hour")

        print("장애물 사용 여부:")
        print("  1) 사용 안 함")
        print("  2) line 장애물 사용")
        print("  3) polygon 장애물 사용")
        barrier_sel = input("입력(기본=1): ").strip()

        barrier_type = None
        if barrier_sel == "2":
            barrier_type = "line"
        elif barrier_sel == "3":
            barrier_type = "polygon"

        if barrier_type is not None:
            barrier_path = filedialog.askopenfilename(
                title=f"{barrier_type} 장애물 GPKG 선택", filetypes=[("GeoPackage", "*.gpkg")])
            if not barrier_path:
                raise RuntimeError("장애물 GPKG를 선택하지 않았습니다.")

        out_gpkg = os.path.join(out_dir, "경로_point_point_편차범위.gpkg")

        # 그래프 로드
        with open(cache_path, "rb") as f:
            payload = pickle.load(f)
        if isinstance(payload, dict) and "graph" in payload:
            G = payload["graph"]
        else:
            G = payload

        routes = parse_routes_csv(routes_csv, mode)
        edges_gdf = build_edge_gdf(G)

        if barrier_type is not None:
            G_work, removed_n = apply_barrier_edges(G, edges_gdf, read_barrier_layer(barrier_path, barrier_type))
            snap_edge_gdf = build_edge_gdf(G_work)
            print(f"[INFO] barrier applied: type={barrier_type}, removed_edges={removed_n}")
        else:
            G_work = G
            snap_edge_gdf = edges_gdf
            print("[INFO] barrier not used")

        print(f"[INFO] graph: nodes={G_work.number_of_nodes():,}, edges={G_work.number_of_edges():,}")

        tf = Transformer.from_crs(4326, 5179, always_xy=True)

        feats = []
        no_in_range = 0
        total = len(routes)

        for i, row in enumerate(routes, start=1):
            route_id = row["route_id"]
            target_value = float(row["target"])

            t_route_start = time.time()
            print(f"\n[{i}/{total}] {route_id} 시작 (target={target_value}, margin=±{margin})")

            snapped_nodes = []
            snap_dists = []
            restore_infos = []
            for x, y in row["coords"]:
                x5179, y5179 = tf.transform(float(x), float(y))
                node, sd, restore = snap_to_graph(x5179, y5179, snap_edge_gdf, G_work)
                if node is None:
                    snapped_nodes = []
                    break
                snapped_nodes.append(node)
                snap_dists.append(float(sd))
                if restore is not None:
                    restore_infos.append(restore)

            if len(snapped_nodes) < 2:
                if restore_infos:
                    restore_splits(G_work, restore_infos)
                print(f"[{i}/{total}] {route_id} NO_SNAP")
                continue

            print(f"  snap_m: {[f'{d:.1f}' for d in snap_dists]}")

            try:
                top_results = choose_target_top_n_margin(
                    G_work, snapped_nodes, target_weight_attr, target_value, margin)
            except (nx.NetworkXNoPath, nx.NodeNotFound):
                top_results = []

            if not top_results:
                no_in_range += 1
                if restore_infos:
                    restore_splits(G_work, restore_infos)
                print(f"[{i}/{total}] {route_id} NO_IN_RANGE")
                continue

            for rank, result in enumerate(top_results, start=1):
                rec = build_record(
                    route_id=route_id,
                    rank=rank,
                    weight_attr=target_weight_attr,
                    target_value=target_value,
                    margin=margin,
                    mode=mode,
                    result=result,
                    snap_start=snap_dists[0],
                    snap_end=snap_dists[-1],
                    G=G_work,
                )
                feats.append(rec)

            # 임시 노드 복원
            if restore_infos:
                restore_splits(G_work, restore_infos)

            elapsed_route = time.time() - t_route_start
            print(f"[{i}/{total}] {route_id} top={len(top_results)} ({elapsed_route:.1f}s)")

        if not feats:
            raise RuntimeError("margin 범위 내 경로가 없습니다.")

        if no_in_range > 0:
            print(f"[WARN] margin 범위 내 경로 없음: {no_in_range}건")

        if os.path.exists(out_gpkg):
            os.remove(out_gpkg)

        gdf = gpd.GeoDataFrame(feats, geometry="geometry", crs="EPSG:5179")
        gdf.to_file(out_gpkg, layer="route_target", driver="GPKG")
        print(f"saved gpkg: {out_gpkg}  (features={len(feats)})")

    except Exception as e:
        print(f"[ERROR] {e}")


if __name__ == "__main__":
    main()