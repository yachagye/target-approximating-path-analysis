# 경로 간 Jaccard 유사도 분석
#
# 입력:
#  1) graph_cache.pkl  (네트워크데이터셋_gpkg_pkl.py 출력물)
#  2) 기본 경로.gpkg   (장애물 미적용 분석 결과)
#  3) 제한 경로.gpkg   (장애물 적용 분석 결과, 선택)
#
# 출력:
#  - 유사도_최소_제한.csv
#     기본 GPKG의 route_target을 기준으로,
#     기본 route_min_* 5개
#     + 제한 route_target 1개 + 제한 route_min_* 5개
#     와의 길이 가중 Jaccard 유사도
#
# 전제:
#  - 모든 비교 레이어는 route_id 필드를 가져야 함
#  - route_id 중복은 경로분석 단계 입력에서 이미 차단되었다고 가정

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
    """축약 그래프의 chain_nodes를 전개하여 전체 노드 집합과 간선별 길이 룩업을 생성.

    축약 간선의 chain_nodes에 포함된 중간 노드를 all_node_set에 추가하고,
    원본 세그먼트 단위(연속 chain node 쌍)의 length_3dkm을 edge_length_lookup에 저장한다.
    비축약 그래프에서도 동일하게 작동한다(chain_nodes 미존재 시 [u, v]로 처리).
    """
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
    """길이 가중 Jaccard 유사도."""
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


def read_routes_as_edge_sets(gpkg_path, layer_name, node_set):
    """레이어에서 route_id별 간선 집합 딕셔너리 반환."""
    gdf = gpd.read_file(gpkg_path, layer=layer_name)

    if gdf.empty:
        return {}

    if "route_id" not in gdf.columns:
        raise RuntimeError(f"{layer_name} 레이어에 route_id 필드가 없습니다.")

    result = {}
    for _, row in gdf.iterrows():
        rid = str(row["route_id"]).strip()
        if rid == "":
            continue

        geom = row.geometry
        if geom is None or geom.is_empty:
            continue

        edges = geometry_to_edge_set(geom, node_set)
        result[rid] = edges

    return result


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

        base_gpkg = filedialog.askopenfilename(
            title="기본 경로.gpkg 선택 (장애물 미적용)",
            filetypes=[("GeoPackage", "*.gpkg")],
        )
        if not base_gpkg:
            raise RuntimeError("기본 경로.gpkg를 선택하지 않았습니다.")

        barrier_gpkg = filedialog.askopenfilename(
            title="제한 경로.gpkg 선택 (장애물 적용, 없으면 취소)",
            filetypes=[("GeoPackage", "*.gpkg")],
        )

        out_dir = filedialog.askdirectory(title="출력 폴더 선택")
        if not out_dir:
            raise RuntimeError("출력 폴더를 선택하지 않았습니다.")

        print("[1/4] 그래프 로드...")
        G = load_graph(cache_path)
        all_node_set, edge_length_lookup = expand_chain_nodes(G)
        print(f"  graph nodes={G.number_of_nodes()}, edges={G.number_of_edges()}")
        print(f"  expanded nodes={len(all_node_set)}, segment edges={len(edge_length_lookup)}")

        print("[2/4] 기본 GPKG route_target 읽기...")
        layers = fiona.listlayers(base_gpkg)
        if "route_target" not in layers:
            raise RuntimeError("기본 GPKG에 route_target 레이어가 없습니다.")

        target_edges = read_routes_as_edge_sets(base_gpkg, "route_target", all_node_set)
        route_ids = sorted(target_edges.keys())
        print(f"  route_target: {len(route_ids)}개 경로")

        print("[3/4] 비교 대상 수집...")
        comparisons = []

        base_min_layers = sorted([l for l in layers if l.startswith("route_min_")])
        for layer_name in base_min_layers:
            edges_dict = read_routes_as_edge_sets(base_gpkg, layer_name, all_node_set)
            comparisons.append((layer_name, edges_dict))
            print(f"  {layer_name}: {len(edges_dict)}개")

        if barrier_gpkg:
            bar_layers = fiona.listlayers(barrier_gpkg)

            if "route_target" in bar_layers:
                edges_dict = read_routes_as_edge_sets(
                    barrier_gpkg, "route_target", all_node_set
                )
                comparisons.append(("route_target_제한", edges_dict))
                print(f"  route_target_제한: {len(edges_dict)}개")

            bar_min_layers = sorted([l for l in bar_layers if l.startswith("route_min_")])
            for layer_name in bar_min_layers:
                edges_dict = read_routes_as_edge_sets(
                    barrier_gpkg, layer_name, all_node_set
                )
                label = layer_name + "_제한"
                comparisons.append((label, edges_dict))
                print(f"  {label}: {len(edges_dict)}개")
        else:
            print("  제한 GPKG 미사용")

        print("[4/4] Jaccard 유사도 산출...")
        comp_labels = [c[0] for c in comparisons]
        header = ["route_id"] + comp_labels

        rows = []
        for rid in route_ids:
            row = [rid]
            t_edges = target_edges[rid]

            for _, edges_dict in comparisons:
                if rid in edges_dict:
                    j = weighted_jaccard(t_edges, edges_dict[rid], edge_length_lookup)
                    row.append(f"{j:.4f}")
                else:
                    row.append("")

            rows.append(row)

        out_csv = os.path.join(out_dir, "유사도_최소_제한.csv")
        with open(out_csv, "w", encoding="utf-8-sig", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(header)
            writer.writerows(rows)

        print(f"saved: {out_csv}")
        print(f"  경로 {len(rows)}개 × 비교 대상 {len(comparisons)}개")

    except Exception as e:
        print(f"[ERROR] {e}")


if __name__ == "__main__":
    main()
