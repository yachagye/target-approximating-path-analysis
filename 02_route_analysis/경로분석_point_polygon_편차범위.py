# 목표값 ±margin 편차범위 내 상위 N위 경로 분석 (포인트-폴리곤 / 경유지 포함)
#
# 입력:
#  1) graph_cache.pkl
#  2) 출발 포인트.csv
#       거리 모드: route_id, x1, y1, ..., xN, yN, km_beg, km_end
#       시간 모드: route_id, x1, y1, ..., xN, yN, hr_beg, hr_end
#     - 목표값 분기(행별):
#         · beg만 입력(end 빈칸) → 단일 목표값 (T = beg) → 편차범위 [T-margin, T+margin]
#         · beg·end 모두 입력 → 구간 [beg, end] → 편차범위 [beg-margin, end+margin]
#         · beg 빈칸 → 행 건너뜀(목표값 필수)
#         · beg > end → 오류로 행 건너뜀(콘솔 보고)
#  3) 목적 폴리곤.gpkg (route_id 필드 필요)
#  4) 선택: 장애물 GPKG (line / polygon)
#
# margin:
#  - 거리 기반(km): ±5리 × 환산계수(km/리, 실행 시 입력)
#  - 시간 기반(hr): 실행 시 사용자 입력 (기본 0.5 hour)
#
# 알고리즘 (편차범위 = 입력 구간을 margin만큼 확장하여 구간 모드로 처리):
#  - 단일 목표값 T는 [T-margin, T+margin], 구간 [beg,end]는 [beg-margin, end+margin]으로
#    확장한 뒤, 확장 구간에 드는 경로를 metric 오름차순으로 수집한다.
#    · 좌표 1개: 구간에 드는 경계 지점(dest_node) 전부 (각 1경로, 다양성 필터 없음)
#    · 경유지 있음: 앞쪽 leg 조합 + 마지막 polygon leg 범위 수집,
#                  전체 경로 노드열 길이가중 Jaccard 다양성 필터(앞쪽 미세 변형 제거)
#  - margin 범위 내 후보가 없으면 해당 경로는 출력하지 않는다(fallback 없음).
#  - abs_diff는 원본 목표값 기준: 단일은 |비용-T|, 구간은 구간거리 max(0, beg-비용, 비용-end).
#  [한계] 경유지 2개 이상(3-leg+)은 abs_diff/metric 순 상위 N의 엄밀성이 보장되지 않는
#         근사임(순위 분석과 동일). rank 경계에서 어긋날 수 있음.
#  ※ 경유지 포함 경로(2-leg 이상)는 leg 병합 시 우회 가능한 구간 왕복(가짜 고리)을
#    가진 후보를 수집에서 제외한다. 진입로가 유일한 막다른 경유지의 불가피한 왕복은 유지한다.
#
# 출력:
#  - 경로_point_polygon_편차범위.gpkg (EPSG:5179)
#     layer: route_target  (경로당 최대 TOP_N개, rank·margin 속성 포함)
#                          단일 모드는 km_end/hr_end가 NULL (target_val에는 beg 기록)

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


TOP_N = 30
K_MAX = 30000
SIM_THRESHOLD = 0.95       # 다양성 필터: 이 값 이상이면 미세 변형으로 간주하여 제외
MARGIN_RI = 5              # 사료 거리 불확실성 ±리 (고정)
RI_TO_KM_DEFAULT = 0.46    # 리→km 환산계수 기본값(km/리, 실행 시 입력)
MARGIN_HR = 0.5            # 시간 기반 편차범위 ±hour 기본값(실행 시 입력 가능)

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

            beg_col = "km_beg" if mode == "km" else "hr_beg"
            end_col = "km_end" if mode == "km" else "hr_end"
            if header[0].strip() != "route_id":
                raise RuntimeError("첫 컬럼은 route_id 이어야 합니다.")
            if header[-2].strip() != beg_col or header[-1].strip() != end_col:
                raise RuntimeError(f"마지막 두 컬럼은 {beg_col}, {end_col} 이어야 합니다.")

            coord_cols = header[1:-2]
            if len(coord_cols) < 2 or len(coord_cols) % 2 != 0:
                raise RuntimeError("좌표 컬럼 구조가 올바르지 않습니다.")

            for row in reader:
                if not row:
                    continue
                if len(row) < len(header):
                    row += [""] * (len(header) - len(row))
                elif len(row) > len(header):
                    row = row[:len(header)]

                route_id = str(row[0]).strip()
                beg_str = str(row[-2]).strip()
                end_str = str(row[-1]).strip()
                if route_id == "":
                    continue
                if beg_str == "":
                    if end_str != "":
                        print(f"  [skip] {route_id}: end만 입력됨 (beg 누락)")
                    else:
                        print(f"  [skip] {route_id}: 목표값 미입력")
                    continue

                coords = []
                for i in range(1, len(row) - 2, 2):
                    xs, ys = str(row[i]).strip(), str(row[i + 1]).strip()
                    if xs == "" or ys == "":
                        continue
                    coords.append((float(xs), float(ys)))

                if len(coords) < 1:
                    continue

                beg = float(beg_str)
                end = float(end_str) if end_str != "" else None
                if end is not None and beg > end:
                    print(f"  [skip] {route_id}: beg({beg}) > end({end})")
                    continue

                out.append({"route_id": route_id, "coords": coords, "beg": beg, "end": end})
        return out

    try:
        rows = _read("utf-8-sig")
    except UnicodeDecodeError:
        rows = _read("cp949")

    if not rows:
        raise RuntimeError("CSV에서 유효한 시작점 데이터를 찾지 못했습니다.")

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

    costs_us = {k: sum(cc[j][k] for j in range(chain_idx)) for k in cost_keys}
    costs_sv = {k: sum(cc[j][k] for j in range(chain_idx, len(cc))) for k in cost_keys}

    fwd_geom = fwd["geom"]
    d = fwd_geom.project(Point(split_node))
    geom_us = substring(fwd_geom, 0, d, normalized=False)
    geom_sv = substring(fwd_geom, d, fwd_geom.length, normalized=False)
    if geom_us.geom_type == "Point":
        geom_us = LineString([geom_us.coords[0], split_node])
    if geom_sv.geom_type == "Point":
        geom_sv = LineString([split_node, geom_sv.coords[-1]])

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


def compute_bridge_set(G):
    """multi-leg 왕복 판정용 위상 정보. 무방향 그래프의 bridge(절단 간선) 집합과
    articulation point(절단점) 집합을 함께 반환한다(경로 분석 직전 경로당 1회 계산).
    반환 객체는 has_artificial_loop에 그대로 전달한다."""
    UG = nx.Graph()
    for u, v in G.edges():
        if not UG.has_edge(u, v):
            UG.add_edge(u, v)
    bridges = set(frozenset(e) for e in nx.bridges(UG))
    arts = set(nx.articulation_points(UG))
    return {"bridges": bridges, "arts": arts}


def has_artificial_loop(nodes, topo):
    """병합 경로의 '가짜 고리'(우회 가능한 구간 왕복) 여부.
    위상 기준으로 불가피한 왕복(허용)과 가짜 고리(배제)를 구별한다.
    - 왕복(2회 통과)된 간선 중 non-bridge가 하나라도 있으면 가짜 → 배제(True).
    - 왕복 간선이 전부 bridge이거나 없으면, 노드 재방문을 검사하되
      bridge 왕복이 설명하는 노드(왕복된 bridge 간선의 양 끝점)와 절단점(articulation)
      에서의 재방문은 불가피한 통과로 보아 허용하고, 그 외 비절단점 재방문만 가짜 → 배제."""
    bridges = topo["bridges"]
    arts = topo["arts"]
    seen = {}
    repeated = set()
    for i in range(len(nodes) - 1):
        e = frozenset((nodes[i], nodes[i + 1]))
        seen[e] = seen.get(e, 0) + 1
        if seen[e] >= 2:
            repeated.add(e)
    # 1) non-bridge 왕복 간선 → 가짜 고리
    if any(e not in bridges for e in repeated):
        return True
    # 2) 노드 재방문 검사 (이 시점 repeated는 전부 bridge)
    explained = set()
    for e in repeated:
        explained |= set(e)
    cnt = {}
    for n in nodes:
        cnt[n] = cnt.get(n, 0) + 1
    for n, c in cnt.items():
        if c >= 2 and n not in explained and n not in arts:
            return True
    return False



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


# ── 폴리곤 / 장애물 ──────────────────────────────────────

def read_polygon_layer(gpkg_path):
    layers = fiona.listlayers(gpkg_path)
    if not layers:
        raise RuntimeError("GPKG에서 레이어를 찾지 못했습니다.")
    if len(layers) == 1:
        layer_name = layers[0]
    else:
        print("GPKG 레이어 목록:")
        for i, lyr in enumerate(layers, start=1):
            print(f"  {i}. {lyr}")
        s = input("사용할 레이어 번호 입력(기본=1): ").strip()
        idx = int(s) if s else 1
        if idx < 1 or idx > len(layers):
            raise RuntimeError("레이어 번호가 유효하지 않습니다.")
        layer_name = layers[idx - 1]

    gdf = gpd.read_file(gpkg_path, layer=layer_name)
    if gdf.empty:
        raise RuntimeError("폴리곤 레이어가 비어 있습니다.")
    if gdf.crs is None:
        raise RuntimeError("폴리곤 레이어 CRS가 없습니다.")
    if "route_id" not in gdf.columns:
        raise RuntimeError("폴리곤 레이어에 route_id 필드가 없습니다.")

    gdf = gdf.to_crs(5179)
    gdf["route_id"] = gdf["route_id"].astype(str)

    dup = gdf["route_id"].duplicated(keep=False)
    if dup.any():
        dups = sorted(gdf.loc[dup, "route_id"].unique().tolist())
        raise RuntimeError(f"폴리곤 route_id 중복: {dups[:10]}")
    return gdf


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
    gdf = gdf[gdf.geometry.notna() & ~gdf.geometry.is_empty].copy()
    if gdf.empty:
        raise RuntimeError("유효한 장애물 geometry가 없습니다.")
    gdf = gdf.to_crs(5179)

    geom_types = set(gdf.geom_type.unique().tolist())
    if barrier_type == "line":
        allowed = {"LineString", "MultiLineString"}
    else:
        allowed = {"Polygon", "MultiPolygon"}
    if not geom_types.issubset(allowed):
        raise RuntimeError(f"장애물 geometry 유형 불일치: {sorted(geom_types)}")
    return gdf


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
            # u, v 모두 영역 안(barrier 밖)이나 chain 중간이 장애물을 통과
            # 첫 진입점(F→T)과 첫 이탈점(T→F)을 찾아 가운데 장애물 구간만 제거
            idx_in = None
            for ci in range(1, len(chain) - 1):
                if _barrier_inside(chain[ci]):
                    idx_in = ci
                    break
            idx_out = None
            if idx_in is not None:
                for ci in range(idx_in + 1, len(chain)):
                    if not _barrier_inside(chain[ci]):
                        idx_out = ci
                        break
            multi_pass = False
            if idx_out is not None:
                for ci in range(idx_out + 1, len(chain)):
                    if _barrier_inside(chain[ci]):
                        multi_pass = True
                        break
            if idx_in is None or idx_out is None or multi_pass:
                # 진입/이탈점 탐색 실패 또는 다회 통과: 통째 제거로 fallback
                remove_whole.append(key)
            elif idx_out == len(chain) - 1:
                # 이탈점이 v 자신: 진입점 1회 분할 후 뒤쪽(장애물) 구간 제거
                split_targets.append((u, v, idx_in, "after"))
            else:
                # 진입·이탈 두 점 분할 후 가운데 구간만 제거
                split_targets.append((u, v, (idx_in, idx_out), "middle"))

    if not remove_whole and not split_targets:
        return G, edges_gdf, 0

    G_blocked = G.copy()

    # 2단계: 축약 간선 분할 후 장애물 내부 구간만 제거
    split_remove = []
    for u, v, chain_idx, remove_side in split_targets:
        if not G_blocked.has_edge(u, v):
            continue
        if remove_side == "middle":
            # (F,F) 통과: 진입점·이탈점 두 점 분할 후 가운데 장애물 구간만 제거
            idx_in, idx_out = chain_idx
            r1 = _split_edge_at_chain_idx(G_blocked, u, v, idx_in)
            p_in = r1["temp_node"]
            if not G_blocked.has_edge(p_in, v):
                continue
            r2 = _split_edge_at_chain_idx(G_blocked, p_in, v, idx_out - idx_in)
            p_out = r2["temp_node"]
            if G_blocked.has_edge(p_in, p_out):
                split_remove.append((p_in, p_out))
            if G_blocked.has_edge(p_out, p_in):
                split_remove.append((p_out, p_in))
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


# ── 폴리곤 경계 목적지 노드 ──────────────────────────────

def _inside_or_boundary(poly, node_xy):
    return poly.covers(Point(node_xy[0], node_xy[1]))


def find_boundary_destination_nodes(poly, edges_gdf, G):
    """폴리곤 경계를 관통하는 간선에서 목적지 후보 노드를 탐색한다.
    축약 간선의 경우 chain_nodes를 순회하여 경계 진입 chain node에서
    간선을 분할하고, 분할된 노드를 목적지 후보로 채택한다.
    반환: (dest_nodes, boundary_restores)"""
    boundary = poly.boundary
    if boundary.is_empty:
        return [], []
    cand_idx = list(edges_gdf.sindex.intersection(boundary.bounds))
    if not cand_idx:
        return [], []

    sub = edges_gdf.iloc[cand_idx]
    dest_nodes = []
    seen = set()
    boundary_restores = []

    for row in sub.itertuples(index=False):
        geom = row.geometry
        if geom is None or geom.is_empty:
            continue
        if not geom.intersects(boundary):
            continue

        u = row.u
        v = row.v

        u_in = _inside_or_boundary(poly, u)
        v_in = _inside_or_boundary(poly, v)

        if (not u_in) and v_in:
            if not G.has_edge(u, v):
                continue

            edge_data = G[u][v]
            chain = edge_data.get("chain_nodes")

            if chain is None or len(chain) <= 2:
                # 비축약 간선 (~30m): 기존 방식
                if v not in seen:
                    seen.add(v)
                    dest_nodes.append(v)
            else:
                # 축약 간선: 첫 번째 내부 chain node 탐색
                first_inside_idx = None
                for ci in range(1, len(chain)):
                    if _inside_or_boundary(poly, chain[ci]):
                        first_inside_idx = ci
                        break

                if first_inside_idx is None or first_inside_idx == len(chain) - 1:
                    # 첫 내부 노드가 v 자체이거나 탐색 실패
                    if v not in seen:
                        seen.add(v)
                        dest_nodes.append(v)
                else:
                    split_node = chain[first_inside_idx]
                    if split_node not in seen:
                        seen.add(split_node)
                        restore = _split_edge_at_chain_idx(G, u, v, first_inside_idx)
                        boundary_restores.append(restore)
                        dest_nodes.append(split_node)

        elif (not u_in) and (not v_in):
            # (F, F) 관통 케이스: 양 끝 모두 폴리곤 외부이나 chain 중간이 폴리곤 내부를 관통
            # geom.intersects(boundary) True가 이미 확인되었으므로 경계 교차는 보장됨
            # 분기 없는 축약 간선이 폴리곤을 가로지를 때 발생
            if not G.has_edge(u, v):
                continue

            edge_data = G[u][v]
            chain = edge_data.get("chain_nodes")

            # 비축약 간선의 관통은 구조적으로 거의 발생하지 않으며,
            # 그래프 노드 수준에서 진입점을 정의할 수 없으므로 처리 불가
            if chain is None or len(chain) <= 2:
                continue

            # chain 중간에서 첫 번째 외부→내부 전이점 탐색 (양 끝 u, v 제외)
            first_inside_idx = None
            for ci in range(1, len(chain) - 1):
                if _inside_or_boundary(poly, chain[ci]):
                    first_inside_idx = ci
                    break

            if first_inside_idx is not None:
                split_node = chain[first_inside_idx]
                if split_node not in seen:
                    seen.add(split_node)
                    restore = _split_edge_at_chain_idx(G, u, v, first_inside_idx)
                    boundary_restores.append(restore)
                    dest_nodes.append(split_node)

    return dest_nodes, boundary_restores


def collect_clean_dest_candidates(G, src_node, target_dests, all_dest_set, weight_attr):
    """모든 dest_node를 가린 읽기 전용 뷰(core)로 폴리곤 내부를 봉인한 뒤 src에서 단일 Dijkstra.
    각 목적지 후보는 폴리곤에 '처음 진입'하는 깨끗한 최단경로로 복원한다.
    (다른 dest_node 경유 = 폴리곤 선진입 → core 봉인으로 자동 배제되며,
     분할로 생성된 A의 외부측 진입 간선은 비용·geometry 메타데이터를 보유하므로
     accumulate_costs가 정확히 누적한다. A의 내부측 이웃은 core에서 도달 불가(∞)가 되어
     별도 판정 없이 자동 제외된다.)
    반환: {dst_node: {nodes, costs, metric, dst_node}} (깨끗한 진입 경로가 없는 후보는 제외)"""
    core = nx.restricted_view(G, all_dest_set - {src_node}, [])
    dist, paths = nx.single_source_dijkstra(core, src_node, weight=weight_attr)

    out = {}
    for A in target_dests:
        if A == src_node:
            continue
        best = None
        for u in G.predecessors(A):
            if u in dist:
                c = dist[u] + float(G[u][A][weight_attr])
                if best is None or c < best[0]:
                    best = (c, u)
        if best is None:
            continue
        upred = best[1]
        nodes = paths[upred] + [A]
        costs = accumulate_costs(G, nodes)
        out[A] = {
            "nodes": nodes,
            "costs": costs,
            "metric": path_metric(costs, weight_attr),
            "dst_node": A,
        }
    return out


def shortest_path_leg(G, src_node, dst_node, weight_attr):
    nodes = nx.shortest_path(G, source=src_node, target=dst_node, weight=weight_attr)
    costs = accumulate_costs(G, nodes)
    metric = path_metric(costs, weight_attr)
    return {"nodes": nodes, "costs": costs, "metric": metric, "dst_node": dst_node}


# ── Top-N 경로 탐색 ──────────────────────────────────────

def _collect_leg_candidates_with_cap(G, src, dst, weight_attr, metric_cap):
    candidates = []
    gen = nx.shortest_simple_paths(G, source=src, target=dst, weight=weight_attr)
    t0 = time.time()
    for path_nodes in gen:
        costs = accumulate_costs(G, path_nodes)
        metric = path_metric(costs, weight_attr)
        candidates.append({"nodes": path_nodes, "costs": costs, "metric": metric, "dst_node": dst})
        k = len(candidates)
        if k % 1000 == 0:
            elapsed = time.time() - t0
            print(f"    [cap_collect] k={k}, metric={metric:.4f}, cap={metric_cap:.4f}, elapsed={elapsed:.1f}s")
        if metric > metric_cap or k >= K_MAX:
            break
    elapsed = time.time() - t0
    print(f"    [cap_collect] 완료: k={len(candidates)}, cap={metric_cap:.4f}, elapsed={elapsed:.1f}s")
    return candidates if candidates else None


# ── 레코드 생성 ──────────────────────────────────────────

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


def _interval_direct_top_n(G, src_node, dest_nodes, weight_attr, beg, end, top_n):
    """좌표 1개 구간 모드: 구간 [beg,end]에 드는 dest_node를 metric 오름차순으로 전부 수집.
    각 dest_node당 경로가 하나뿐이라 미세 변형이 없으므로 다양성 필터를 적용하지 않는다.
    범위에 드는 dest_node가 없으면 제외한다(fallback 없음)."""
    dest_nodes_set = set(dest_nodes)
    t0 = time.time()
    clean = collect_clean_dest_candidates(G, src_node, dest_nodes, dest_nodes_set, weight_attr)
    all_cands = list(clean.values())
    print(f"    [dest_scan] polygon dest={len(dest_nodes)}, clean={len(all_cands)}, elapsed={time.time()-t0:.1f}s")
    if not all_cands:
        return []

    in_range = [c for c in all_cands if beg <= float(c["metric"]) <= end]

    if not in_range:
        # 편차범위: 범위 내 후보 없으면 제외(fallback 없음)
        print(f"    [margin] 범위 [{beg}, {end}] 내 dest 없음 → 제외")
        return []

    in_range.sort(key=lambda c: float(c["metric"]))
    results = []
    for c in in_range[:top_n]:
        results.append({
            "nodes": c["nodes"], "costs": c["costs"], "metric": float(c["metric"]),
            "dst_node": c["dst_node"], "candidate_n": len(in_range), "legs_n": 1,
        })
    return results


def _interval_via_waypoints_top_n(G, snapped_nodes, dest_nodes, weight_attr, beg, end, top_n):
    """경유지 구간 모드: 앞쪽 leg 최소 고정 + 마지막 polygon leg를 [beg,end] 범위 수집.
    전체 경로 노드열 Jaccard 다양성 필터(앞쪽 point-to-point leg 미세 변형 제거) 후 rank.
    범위에 드는 조합이 없으면 제외한다(fallback 없음)."""
    front_leg_pairs = list(zip(snapped_nodes[:-1], snapped_nodes[1:]))
    dest_nodes_set = set(dest_nodes)
    topo = compute_bridge_set(G)  # 가짜 고리 판정용 (경로당 1회)

    front_min_metrics = []
    for src, dst in front_leg_pairs:
        leg_min = shortest_path_leg(G, src, dst, weight_attr)
        front_min_metrics.append(float(leg_min["metric"]))

    t0 = time.time()
    clean_by_node = collect_clean_dest_candidates(G, snapped_nodes[-1], dest_nodes, dest_nodes_set, weight_attr)
    last_min_by_node = {dst: float(r["metric"]) for dst, r in clean_by_node.items()}
    print(f"    [dest_scan] polygon dest={len(dest_nodes)}, clean={len(last_min_by_node)}, elapsed={time.time()-t0:.1f}s")
    if not last_min_by_node:
        return []

    last_min_metric = min(last_min_by_node.values())
    s_min = sum(front_min_metrics) + last_min_metric

    # t_last 사전 필터: (T+범위 상한) 기준으로 end 사용
    t_last = float(end) - (s_min - last_min_metric)
    filtered_dest = [dst for dst, m in last_min_by_node.items() if m <= t_last]
    if not filtered_dest:
        filtered_dest = [min(last_min_by_node, key=last_min_by_node.get)]

    last_candidates = sorted(
        (clean_by_node[d] for d in filtered_dest),
        key=lambda x: (float(x["metric"]), x["dst_node"]),
    )
    print(f"    [polygon_leg] filtered={len(filtered_dest)}, candidates={len(last_candidates)}")
    if not last_candidates:
        return []

    lm = [float(c["metric"]) for c in last_candidates]
    N_l = len(lm)
    F = len(front_leg_pairs)
    mid = (beg + end) / 2.0

    # [beg, end] 범위 조합을 total 오름차순으로 방출하는 min-heap 병합(k-smallest-pairs).
    # 앞쪽(front)을 c_front 오름차순으로 지연 pull, polygon leg(lm)을 col로 둔다. 각 행은
    # lo(범위 하한 통과 첫 col)에서 seed하며, 하한 때문에 행의 in-range 최소 total이 c_front에
    # 단조가 아니므로, "다음 행 하한(c_front + min_l)이 현재 heap 최소 total 이하인 동안 미리
    # pull"하는 lower-bound 방식으로 전역 total 오름차순을 보장한다.
    # (편차범위: 범위 내 조합이 없으면 제외, fallback 없음.)
    min_l = lm[0]
    front_pool = []          # [(c_front, front_obj), ...]
    leg_cands = None
    front_done = False
    _buf = []                # F==1 미리보기 버퍼 (길이 0 또는 1)
    t0 = time.time()

    if F == 1:
        _front_gen = nx.shortest_simple_paths(
            G, source=front_leg_pairs[0][0], target=front_leg_pairs[0][1], weight=weight_attr)
        _front_cnt = 0

        def _front_peek():
            nonlocal front_done, _front_cnt
            if _buf:
                return _buf[0][0]
            if front_done:
                return None
            try:
                pn = next(_front_gen)
                _front_cnt += 1
            except StopIteration:
                front_done = True
                return None
            cst = accumulate_costs(G, pn)
            cf = path_metric(cst, weight_attr)
            if _front_cnt % 1000 == 0:
                print(f"    [front_gen] front 후보 {_front_cnt}개 생성, cf={cf:.4f}, elapsed={time.time() - t0:.1f}s")
            _buf.append((cf, {"nodes": pn, "costs": cst, "metric": cf}))
            return cf

        def _front_take():
            return _buf.pop(0)
    else:
        leg_cands = []
        for idx, (src, dst) in enumerate(front_leg_pairs):
            t_j = float(end) - (s_min - front_min_metrics[idx])
            cands = _collect_leg_candidates_with_cap(G, src, dst, weight_attr, t_j)
            if not cands:
                return []
            leg_cands.append(cands)
        front_metrics = [[float(c["metric"]) for c in cs] for cs in leg_cands]
        front_sizes = [len(m) for m in front_metrics]
        start = tuple([0] * F)
        start_sum = sum(front_metrics[j][0] for j in range(F))
        combo_heap = [(start_sum, start)]
        visited = {start}

        def _front_peek():
            return combo_heap[0][0] if combo_heap else None

        def _front_take():
            c_front, state = heapq.heappop(combo_heap)
            for k in range(F):
                ni = state[k] + 1
                if ni < front_sizes[k]:
                    ns = list(state)
                    ns[k] = ni
                    ns = tuple(ns)
                    if ns not in visited:
                        visited.add(ns)
                        new_sum = c_front - front_metrics[k][state[k]] + front_metrics[k][ni]
                        heapq.heappush(combo_heap, (new_sum, ns))
            return c_front, state

    merge_heap = []  # (total, fi, j)

    def _take_and_seed():
        cf, obj = _front_take()
        fi = len(front_pool)
        front_pool.append((cf, obj))
        lo = _bisect.bisect_left(lm, beg - cf)
        if lo < N_l and (cf + lm[lo]) <= end:
            heapq.heappush(merge_heap, (cf + lm[lo], fi, lo))

    def assemble(front_obj, j):
        if F == 1:
            chosen = [front_obj, last_candidates[j]]
        else:
            chosen = [leg_cands[k][front_obj[k]] for k in range(F)] + [last_candidates[j]]
        merged = merge_leg_paths(chosen)
        tcosts = sum_cost_dicts([x["costs"] for x in chosen])
        tmetric = sum(float(x["metric"]) for x in chosen)
        return chosen, merged, tcosts, tmetric

    selected = []
    selected_es = []
    popped = 0

    while True:
        # lower-bound pull-ahead: 다음 행 하한이 현재 heap 최소 이하인 동안 미리 pull
        while True:
            nc = _front_peek()
            if nc is None:
                break
            if (nc + min_l) > end:
                break
            if F == 1 and len(front_pool) >= K_MAX:
                front_done = True
                _buf.clear()
                break
            if merge_heap and (nc + min_l) > merge_heap[0][0]:
                break
            _take_and_seed()

        if not merge_heap:
            break

        total, fi, j = heapq.heappop(merge_heap)
        if total > end:
            break
        popped += 1
        if popped % (1000 if F == 1 else 10000) == 0:
            print(
                f"    [iterate] popped={popped}, total={total:.4f}, heap={len(merge_heap)}, fronts={len(front_pool)}, selected={len(selected)}, elapsed={time.time() - t0:.1f}s")

        cf, front_obj = front_pool[fi]
        if (j + 1) < N_l and (cf + lm[j + 1]) <= end:
            heapq.heappush(merge_heap, (cf + lm[j + 1], fi, j + 1))

        if total >= beg:
            chosen, merged, tcosts, tmetric = assemble(front_obj, j)
            if has_artificial_loop(merged, topo):
                continue  # 가짜 고리(non-bridge 왕복) 제외
            es = path_edge_set(G, merged)
            if _passes_diversity(es, selected_es):
                selected.append({"nodes": merged, "costs": tcosts, "metric": tmetric,
                                 "dst_node": last_candidates[j]["dst_node"],
                                 "candidate_n": 0, "legs_n": len(chosen)})
                selected_es.append(es)
                if len(selected) >= top_n:
                    break

    candidate_n = len(front_pool) + len(last_candidates)
    print(
        f"    [iterate] 완료: popped={popped}, fronts={len(front_pool)}, selected={len(selected)}, elapsed={time.time() - t0:.1f}s")

    if not selected:
        print(f"    [margin] 범위 내 조합 없음 → 제외")
        return []

    for s in selected:
        s["candidate_n"] = candidate_n

    return selected


def choose_target_interval_top_n(G, snapped_nodes, dest_nodes, weight_attr, beg, end, top_n=TOP_N):
    """구간 [beg,end] 모드 순위 수집 진입점."""
    if len(snapped_nodes) == 1:
        print(f"  [interval] direct: weight={weight_attr}, [{beg:.4f}, {end:.4f}], top_n={top_n}, dest={len(dest_nodes)}")
        return _interval_direct_top_n(G, snapped_nodes[0], dest_nodes, weight_attr, beg, end, top_n)
    else:
        print(f"  [interval] legs={len(snapped_nodes)}(front={len(snapped_nodes)-1}+polygon=1), [{beg:.4f}, {end:.4f}], top_n={top_n}")
        return _interval_via_waypoints_top_n(G, snapped_nodes, dest_nodes, weight_attr, beg, end, top_n)


def build_record(route_id, rank, weight_attr, beg, end, mode, result, snap_start, snap_last, coord_n, G, margin=None):
    geom = path_to_linestring(G, result["nodes"])
    costs = result["costs"]
    dst_x, dst_y = result["dst_node"]

    metric = float(result["metric"])
    # abs_diff: 단일(end None) → |metric-beg| / 구간 → 구간거리
    if end is None:
        abs_diff = float(abs(metric - beg))
    else:
        abs_diff = float(max(0.0, beg - metric, metric - end))

    rec = {
        "route_id": route_id,
        "rank": int(rank),
        "weight_attr": weight_attr,
        "target_val": float(beg),
        "metric_val": metric,
        "abs_diff": abs_diff,
        "length_km": float(costs["length_km"]),
        "hour_ks": float(costs["hour_ks"]),
        "hour_tob": float(costs["hour_tob"]),
        "kcal_ks": float(costs["kcal_ks"]),
        "kcal_tob": float(costs["kcal_tob"]),
        "coord_n": int(coord_n),
        "legs_n": int(result["legs_n"]),
        "candidate_n": int(result["candidate_n"]),
        "snap_m_start": float(snap_start),
        "snap_m_last": float(snap_last),
        "dst_x": float(dst_x),
        "dst_y": float(dst_y),
        "geometry": geom,
    }

    if mode == "km":
        rec["km_beg"] = float(beg)
        rec["km_end"] = float(end) if end is not None else None
    else:
        rec["hr_beg"] = float(beg)
        rec["hr_end"] = float(end) if end is not None else None

    if margin is not None:
        rec["margin"] = float(margin)

    return rec


# ── 메인 ─────────────────────────────────────────────────

def main():
    try:
        root = tk.Tk()
        root.withdraw()

        cache_path = filedialog.askopenfilename(title="graph_cache.pkl 선택", filetypes=[("Pickle", "*.pkl")])
        if not cache_path:
            raise RuntimeError("graph_cache.pkl을 선택하지 않았습니다.")

        routes_csv = filedialog.askopenfilename(
            title="출발 포인트.csv 선택", filetypes=[("CSV", "*.csv")])
        if not routes_csv:
            raise RuntimeError("출발 포인트.csv를 선택하지 않았습니다.")

        polygon_gpkg = filedialog.askopenfilename(
            title="목적 폴리곤 GPKG 선택 (route_id 필드 필요)", filetypes=[("GeoPackage", "*.gpkg")])
        if not polygon_gpkg:
            raise RuntimeError("목적 폴리곤 GPKG를 선택하지 않았습니다.")

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

        # 편차범위 margin: 거리=±5리 × 환산계수(실행 시 입력), 시간=사용자 입력
        if mode == "km":
            cs = input(f"리→km 환산계수 입력(km/리, 기본={RI_TO_KM_DEFAULT}): ").strip()
            ri_to_km = float(cs) if cs else RI_TO_KM_DEFAULT
            margin = MARGIN_RI * ri_to_km
            print(f"[INFO] 거리 기반 margin = ±{margin:.2f} km (={MARGIN_RI}리 × {ri_to_km} km/리)")
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
        if barrier_sel == "2":
            barrier_type = "line"
        elif barrier_sel == "3":
            barrier_type = "polygon"

        if barrier_type is not None:
            barrier_path = filedialog.askopenfilename(
                title=f"{barrier_type} 장애물 GPKG 선택", filetypes=[("GeoPackage", "*.gpkg")])
            if not barrier_path:
                raise RuntimeError("장애물 GPKG를 선택하지 않았습니다.")

        out_gpkg = os.path.join(out_dir, "경로_point_polygon_편차범위.gpkg")

        # 그래프 로드
        with open(cache_path, "rb") as f:
            payload = pickle.load(f)
        if isinstance(payload, dict) and "graph" in payload:
            G = payload["graph"]
        else:
            G = payload

        routes = parse_routes_csv(routes_csv, mode)
        poly_gdf = read_polygon_layer(polygon_gpkg)
        poly_map = {str(row.route_id): row.geometry for row in poly_gdf.itertuples(index=False)}

        missing_poly = [r["route_id"] for r in routes if r["route_id"] not in poly_map]
        if missing_poly:
            raise RuntimeError(f"CSV route_id에 대응하는 폴리곤이 없습니다: {missing_poly[:10]}")

        edges_gdf = build_edge_gdf(G)

        if barrier_type is not None:
            barrier_gdf = read_barrier_layer(barrier_path, barrier_type)
            G_work, snap_edge_gdf, removed_n = apply_barrier_edges(G, edges_gdf, barrier_gdf)
            print(f"[INFO] barrier applied: type={barrier_type}, removed_edges={removed_n}")
        else:
            G_work = G
            snap_edge_gdf = edges_gdf
            print("[INFO] barrier not used")

        print(f"[INFO] graph: nodes={G_work.number_of_nodes():,}, edges={G_work.number_of_edges():,}")

        tf = Transformer.from_crs(4326, 5179, always_xy=True)

        feats = []
        total = len(routes)

        for i, row in enumerate(routes, start=1):
            route_id = row["route_id"]
            beg = row["beg"]
            end = row["end"]   # None 가능 (단일 모드)
            poly = poly_map[route_id]

            t_route_start = time.time()
            tdisp = f"단일 T={beg}" if end is None else f"구간 [{beg}, {end}]"
            print(f"\n[{i}/{total}] {route_id} 시작 ({tdisp}, coords={len(row['coords'])}개)")

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

            if len(snapped_nodes) < 1:
                if restore_infos:
                    restore_splits(G_work, restore_infos)
                print(f"[{i}/{total}] {route_id} NO_SNAP")
                continue

            print(f"  snap_m: {[f'{d:.1f}' for d in snap_dists]}")

            dest_nodes, boundary_restores = find_boundary_destination_nodes(poly, snap_edge_gdf, G_work)
            if not dest_nodes:
                if boundary_restores:
                    restore_splits(G_work, boundary_restores)
                if restore_infos:
                    restore_splits(G_work, restore_infos)
                print(f"[{i}/{total}] {route_id} NO_BOUNDARY_DEST")
                continue

            print(f"  dest_nodes: {len(dest_nodes)}")

            try:
                # 편차범위: 단일은 [beg-margin, beg+margin], 구간은 [beg-margin, end+margin]
                lo = beg - margin
                hi = (beg if end is None else end) + margin
                top_results = choose_target_interval_top_n(
                    G_work, snapped_nodes, dest_nodes, target_weight_attr, lo, hi)
            except (nx.NetworkXNoPath, nx.NodeNotFound):
                top_results = []

            if not top_results:
                if boundary_restores:
                    restore_splits(G_work, boundary_restores)
                if restore_infos:
                    restore_splits(G_work, restore_infos)
                print(f"[{i}/{total}] {route_id} NO_PATH")
                continue

            for rank, result in enumerate(top_results, start=1):
                rec = build_record(
                    route_id=route_id,
                    rank=rank,
                    weight_attr=target_weight_attr,
                    beg=beg,
                    end=end,
                    mode=mode,
                    result=result,
                    snap_start=snap_dists[0],
                    snap_last=snap_dists[-1],
                    coord_n=len(snapped_nodes),
                    G=G_work,
                    margin=margin,
                )
                feats.append(rec)

            # 임시 노드 복원
            if boundary_restores:
                restore_splits(G_work, boundary_restores)
            if restore_infos:
                restore_splits(G_work, restore_infos)

            elapsed_route = time.time() - t_route_start
            print(f"[{i}/{total}] {route_id} top={len(top_results)} ({elapsed_route:.1f}s)")

        if not feats:
            raise RuntimeError("생성된 경로가 없습니다.")

        if os.path.exists(out_gpkg):
            os.remove(out_gpkg)

        gdf = gpd.GeoDataFrame(feats, geometry="geometry", crs="EPSG:5179")
        gdf.to_file(out_gpkg, layer="route_target", driver="GPKG")
        print(f"saved gpkg: {out_gpkg}  (features={len(feats)})")

    except Exception as e:
        print(f"[ERROR] {e}")


if __name__ == "__main__":
    main()