# 순위 경로 rank 간 Jaccard 유사도 분석
#
# 입력:
#  1) graph_cache.pkl
#  2) 경로_top.gpkg (route_target 레이어, route_id + rank 필드)
#
# 비교 구조:
#  - 각 route_id 내에서 rank 1(최근접 경로)을 기준선으로 고정
#  - rank 2 ~ rank N과의 길이 가중 Jaccard 유사도를 산출
#  - Jaccard >= SIMILARITY_THRESHOLD(0.99)인 rank는 미세 변형으로 간주하여 제외
#  - 새 후보가 이미 선택된 모든 경로(rank 1 포함)와 임계값 미만이어야 채택
#  - 필터 통과 rank에 filtered_rank(1~)를 부여, 상위 FILTER_TOP_N개 출력
#
# 출력:
#  - 유사도_순위.csv (long format)
#     route_id, filtered_rank, original_rank, jaccard
#  - 유사도_순위_필터.gpkg
#     필터 통과 피처 + rank 1 기준선 (filtered_rank=0)

import os
import csv
import pickle
import tkinter as tk
from tkinter import filedialog
import geopandas as gpd
import fiona
from shapely.geometry import LineString, MultiLineString

NODE_ROUND_M = 0.01


def _round_xy(x, y, r=NODE_ROUND_M):
    return (round(x / r) * r, round(y / r) * r)


def load_graph(cache_path):
    with open(cache_path, "rb") as f:
        payload = pickle.load(f)
    if isinstance(payload, dict) and "graph" in payload:
        return payload["graph"]
    return payload


def expand_chain_nodes(G):
    """축약 그래프의 chain_nodes를 전개하여 전체 노드 집합과 간선별 길이 룩업을 생성."""
    all_node_set = set(G.nodes())
    edge_length_lookup = {}

    for u, v, data in G.edges(data=True):
        chain = data.get("chain_nodes", [u, v])
        costs = data.get("chain_costs")

        if costs is None:
            costs = [{"length_3dkm": float(data["length_3dkm"])}]

        for node in chain:
            all_node_set.add(node)

        for i in range(len(chain) - 1):
            seg_u, seg_v = chain[i], chain[i + 1]
            edge_length_lookup[(seg_u, seg_v)] = float(costs[i]["length_3dkm"])

    return all_node_set, edge_length_lookup


def geometry_to_edge_set(geom, node_set):
    if geom is None or geom.is_empty:
        return set()

    if isinstance(geom, LineString):
        parts = [geom]
    elif isinstance(geom, MultiLineString):
        parts = list(geom.geoms)
    else:
        return set()

    edge_set = set()
    for part in parts:
        coords = list(part.coords)
        nodes = []
        for coord in coords:
            rounded = _round_xy(float(coord[0]), float(coord[1]))
            if rounded in node_set:
                if not nodes or nodes[-1] != rounded:
                    nodes.append(rounded)

        if len(nodes) >= 2:
            edge_set.update(zip(nodes[:-1], nodes[1:]))

    return edge_set


def weighted_jaccard(edges_a, edges_b, edge_length_lookup):
    if not edges_a and not edges_b:
        return 0.0

    intersection = edges_a & edges_b
    union = edges_a | edges_b

    def sum_length(edge_set):
        total = 0.0
        for u, v in edge_set:
            if (u, v) in edge_length_lookup:
                total += edge_length_lookup[(u, v)]
            elif (v, u) in edge_length_lookup:
                total += edge_length_lookup[(v, u)]
        return total

    union_km = sum_length(union)
    if union_km <= 0:
        return 0.0

    return sum_length(intersection) / union_km


def read_ranked_routes(gpkg_path, node_set):
    """route_target 레이어에서 {route_id: {rank: edge_set}} 및 GeoDataFrame 반환."""
    layers = fiona.listlayers(gpkg_path)
    if "route_target" not in layers:
        raise RuntimeError("GPKG에 route_target 레이어가 없습니다.")

    gdf = gpd.read_file(gpkg_path, layer="route_target")
    if gdf.empty:
        raise RuntimeError("route_target 레이어가 비어 있습니다.")

    for col in ("route_id", "rank"):
        if col not in gdf.columns:
            raise RuntimeError(f"route_target 레이어에 {col} 필드가 없습니다.")

    result = {}
    for _, row in gdf.iterrows():
        rid = str(row["route_id"]).strip()
        rank = int(row["rank"])
        geom = row.geometry
        if rid == "" or geom is None or geom.is_empty:
            continue

        edges = geometry_to_edge_set(geom, node_set)
        if rid not in result:
            result[rid] = {}
        result[rid][rank] = edges

    return result, gdf


def main():
    try:
        root = tk.Tk()
        root.withdraw()

        cache_path = filedialog.askopenfilename(
            title="graph_cache.pkl 선택", filetypes=[("Pickle", "*.pkl")])
        if not cache_path:
            raise RuntimeError("graph_cache.pkl을 선택하지 않았습니다.")

        gpkg_path = filedialog.askopenfilename(
            title="경로_순위.gpkg 선택", filetypes=[("GeoPackage", "*.gpkg")])
        if not gpkg_path:
            raise RuntimeError("GPKG를 선택하지 않았습니다.")

        out_dir = filedialog.askdirectory(title="출력 폴더 선택")
        if not out_dir:
            raise RuntimeError("출력 폴더를 선택하지 않았습니다.")

        print("[1/3] 그래프 로드...")
        G = load_graph(cache_path)
        all_node_set, edge_length_lookup = expand_chain_nodes(G)
        print(f"  graph nodes={G.number_of_nodes()}, edges={G.number_of_edges()}")
        print(f"  expanded nodes={len(all_node_set)}, segment edges={len(edge_length_lookup)}")

        print("[2/3] ranked 경로 읽기...")
        ranked, gdf = read_ranked_routes(gpkg_path, all_node_set)
        route_ids = sorted(ranked.keys())
        print(f"  route_id: {len(route_ids)}개")

        # 각 route_id별 실제 비교 rank 수 확인
        max_compare_n = 0
        for rid in route_ids:
            ranks = ranked[rid]
            if 1 not in ranks:
                print(f"  [WARN] {rid}: rank 1 없음 — 건너뜀")
                continue
            n = len([r for r in ranks if r != 1])
            if n > max_compare_n:
                max_compare_n = n

        if max_compare_n < 1:
            raise RuntimeError("비교 대상 rank가 없습니다 (모든 경로가 rank 1만 보유).")

        SIMILARITY_THRESHOLD = 0.99
        FILTER_TOP_N = 30

        print(f"[3/3] Jaccard 유사도 산출 (최대 비교 rank {max_compare_n}개)...")
        print(f"  필터: 기선택 경로 대비 Jaccard >= {SIMILARITY_THRESHOLD} 제외, 상위 {FILTER_TOP_N}개 출력")
        rows = []
        filter_map = {}  # (route_id, original_rank) -> filtered_rank
        single_rank_count = 0
        no_valid_count = 0
        for rid in route_ids:
            ranks = ranked[rid]
            if 1 not in ranks:
                continue

            if len(ranks) == 1:
                single_rank_count += 1
                continue

            baseline = ranks[1]
            existing = sorted(r for r in ranks if r != 1)

            filtered = []
            selected_edges = [baseline]  # rank 1 포함
            for r in existing:
                candidate = ranks[r]
                # 이미 선택된 모든 경로와 비교
                is_diverse = True
                for selected in selected_edges:
                    if weighted_jaccard(selected, candidate, edge_length_lookup) >= SIMILARITY_THRESHOLD:
                        is_diverse = False
                        break
                if is_diverse:
                    j_vs_rank1 = weighted_jaccard(baseline, candidate, edge_length_lookup)
                    filtered.append((r, j_vs_rank1))
                    selected_edges.append(candidate)
                    if len(filtered) >= FILTER_TOP_N:
                        break

            if not filtered:
                no_valid_count += 1
                continue

            for frank, (orig_rank, jaccard) in enumerate(filtered, start=1):
                rows.append([rid, frank, orig_rank, f"{jaccard:.4f}"])
                filter_map[(rid, orig_rank)] = frank

        # CSV 출력
        header = ["route_id", "filtered_rank", "original_rank", "jaccard"]
        out_csv = os.path.join(out_dir, "유사도_순위.csv")
        with open(out_csv, "w", encoding="utf-8-sig", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(header)
            writer.writerows(rows)

        print(f"saved: {out_csv}")
        print(f"  유효 행 {len(rows)}개 (경로 {len(set(r[0] for r in rows))}개)")
        if single_rank_count > 0:
            print(f"  [INFO] rank 1만 존재: {single_rank_count}건")
        if no_valid_count > 0:
            print(f"  [INFO] 임계값 미만 경로 없음 (전체 >= {SIMILARITY_THRESHOLD}): {no_valid_count}건")

        # 필터 GPKG 출력
        if filter_map:
            # rank 1 기준선 포함 (filtered_rank=0)
            for rid in set(k[0] for k in filter_map):
                filter_map[(rid, 1)] = 0

            gdf["_rid"] = gdf["route_id"].astype(str).str.strip()
            gdf["_rank"] = gdf["rank"].astype(int)
            gdf["_key"] = list(zip(gdf["_rid"], gdf["_rank"]))
            mask = gdf["_key"].isin(filter_map)
            gdf_out = gdf[mask].copy()
            gdf_out["filtered_rank"] = gdf_out["_key"].map(filter_map)
            gdf_out = gdf_out.drop(columns=["_rid", "_rank", "_key"])
            gdf_out = gdf_out.sort_values(["route_id", "filtered_rank"])
            out_gpkg = os.path.join(out_dir, "유사도_순위_필터.gpkg")
            gdf_out.to_file(out_gpkg, layer="route_target", driver="GPKG")
            print(f"saved: {out_gpkg} ({len(gdf_out)}개 피처)")

    except Exception as e:
        print(f"[ERROR] {e}")


if __name__ == "__main__":
    main()
