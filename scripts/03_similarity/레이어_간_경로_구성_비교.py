# 순수 레이어 간 경로 구성 비교 (길이 가중 Jaccard)
#
# 입력:
#  1) graph_cache.pkl
#  2) 타겟 경로.gpkg   (레이어 1개)
#  3) 비교 경로.gpkg   (레이어 1개)
#
# 출력:
#  - 레이어 비교_summary.csv
#  - 레이어 비교_match.csv
#
# 비교 방식:
#  - 속성/ID 매칭 없음
#  - geometry -> edge set 변환 후 길이 가중 Jaccard 계산
#  - 타겟 레이어 각 피처에 대해 비교 레이어 내 최유사 피처 1개 선택
#  - 레이어 전체 union edge set 기준 구성 차이 산출

import os
import csv
import pickle
import tkinter as tk
from tkinter import filedialog

import fiona
import geopandas as gpd
from shapely.geometry import LineString, MultiLineString


NODE_ROUND_M = 0.01
TARGET_CRS = "EPSG:5179"


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


def get_single_layer_name(gpkg_path):
    layers = fiona.listlayers(gpkg_path)
    if len(layers) != 1:
        raise RuntimeError(
            f"{os.path.basename(gpkg_path)} 에 레이어가 1개가 아닙니다. "
            f"(현재 {len(layers)}개: {list(layers)})"
        )
    return layers[0]


def edge_length_km(edge_set, edge_length_lookup):
    total = 0.0
    for u, v in edge_set:
        if (u, v) in edge_length_lookup:
            total += edge_length_lookup[(u, v)]
        elif (v, u) in edge_length_lookup:
            total += edge_length_lookup[(v, u)]
    return total


def weighted_jaccard_and_parts(edges_a, edges_b, edge_length_lookup):
    intersection = edges_a & edges_b
    only_a = edges_a - edges_b
    only_b = edges_b - edges_a
    union = edges_a | edges_b

    shared_km = edge_length_km(intersection, edge_length_lookup)
    only_a_km = edge_length_km(only_a, edge_length_lookup)
    only_b_km = edge_length_km(only_b, edge_length_lookup)
    union_km = edge_length_km(union, edge_length_lookup)

    if union_km <= 0:
        jaccard = 0.0
    else:
        jaccard = shared_km / union_km

    return jaccard, shared_km, only_a_km, only_b_km


def geometry_to_edge_set(geom, node_set):
    if geom is None or geom.is_empty:
        return set()

    edge_set = set()

    if isinstance(geom, LineString):
        parts = [geom]
    elif isinstance(geom, MultiLineString):
        parts = list(geom.geoms)
    else:
        return set()

    for part in parts:
        coords = list(part.coords)
        nodes = []

        for coord in coords:
            x = float(coord[0])
            y = float(coord[1])
            rounded = _round_xy(x, y)

            if rounded in node_set:
                if not nodes or nodes[-1] != rounded:
                    nodes.append(rounded)

        if len(nodes) >= 2:
            edge_set.update(zip(nodes[:-1], nodes[1:]))

    return edge_set


def read_layer_features(gpkg_path, prefix, node_set, edge_length_lookup):
    layer_name = get_single_layer_name(gpkg_path)
    gdf = gpd.read_file(gpkg_path, layer=layer_name)

    if gdf.empty:
        raise RuntimeError(f"{os.path.basename(gpkg_path)} 의 레이어가 비어 있습니다.")

    if gdf.crs is None:
        raise RuntimeError(f"{os.path.basename(gpkg_path)} 의 CRS가 없습니다.")

    if str(gdf.crs).upper() != TARGET_CRS:
        gdf = gdf.to_crs(TARGET_CRS)

    features = []
    has_route_id = "route_id" in gdf.columns

    for idx, row in gdf.iterrows():
        geom = row.geometry
        if geom is None or geom.is_empty:
            continue

        edges = geometry_to_edge_set(geom, node_set)
        feature_id = f"{prefix}_{len(features) + 1:03d}"
        length_km = edge_length_km(edges, edge_length_lookup)
        route_id = str(row["route_id"]) if has_route_id else ""

        features.append(
            {
                "feature_id": feature_id,
                "route_id": route_id,
                "source_index": idx,
                "geometry": geom,
                "edges": edges,
                "length_km": length_km,
            }
        )

    if not features:
        raise RuntimeError(
            f"{os.path.basename(gpkg_path)} 에서 비교 가능한 선형 피처를 읽지 못했습니다."
        )

    return layer_name, features


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

        target_gpkg = filedialog.askopenfilename(
            title="타겟 경로.gpkg 선택",
            filetypes=[("GeoPackage", "*.gpkg")],
        )
        if not target_gpkg:
            raise RuntimeError("타겟 경로.gpkg를 선택하지 않았습니다.")

        compare_gpkg = filedialog.askopenfilename(
            title="비교 경로.gpkg 선택",
            filetypes=[("GeoPackage", "*.gpkg")],
        )
        if not compare_gpkg:
            raise RuntimeError("비교 경로.gpkg를 선택하지 않았습니다.")

        out_dir = filedialog.askdirectory(title="출력 폴더 선택")
        if not out_dir:
            raise RuntimeError("출력 폴더를 선택하지 않았습니다.")

        print("[1/5] 그래프 로드...")
        G = load_graph(cache_path)
        all_node_set, edge_length_lookup = expand_chain_nodes(G)
        print(f"  graph nodes={G.number_of_nodes()}, edges={G.number_of_edges()}")
        print(f"  expanded nodes={len(all_node_set)}, segment edges={len(edge_length_lookup)}")

        print("[2/5] 타겟 레이어 읽기...")
        target_layer_name, target_features = read_layer_features(
            target_gpkg, "T", all_node_set, edge_length_lookup
        )
        print(f"  layer={target_layer_name}, features={len(target_features)}")

        print("[3/5] 비교 레이어 읽기...")
        compare_layer_name, compare_features = read_layer_features(
            compare_gpkg, "C", all_node_set, edge_length_lookup
        )
        print(f"  layer={compare_layer_name}, features={len(compare_features)}")

        print("[4/5] 레이어 전체 구성 비교...")
        target_union = set()
        for feat in target_features:
            target_union.update(feat["edges"])

        compare_union = set()
        for feat in compare_features:
            compare_union.update(feat["edges"])

        layer_jaccard, shared_km, target_only_km, compare_only_km = weighted_jaccard_and_parts(
            target_union, compare_union, edge_length_lookup
        )

        target_feature_n = len(target_features)
        compare_feature_n = len(compare_features)
        feature_n_diff = compare_feature_n - target_feature_n

        print("[5/5] 타겟 피처별 최유사 상대 탐색...")
        best_match_rows = []

        total = len(target_features)
        for i, t_feat in enumerate(target_features, start=1):
            print(f"  [{i}/{total}] {t_feat['feature_id']} 비교 중...")

            best_compare_id = ""
            best_compare_route_id = ""
            best_jaccard = -1.0
            best_shared_km = 0.0
            best_target_only_km = 0.0
            best_compare_only_km = 0.0

            for c_feat in compare_features:
                jaccard, s_km, t_only_km, c_only_km = weighted_jaccard_and_parts(
                    t_feat["edges"], c_feat["edges"], edge_length_lookup
                )

                if jaccard > best_jaccard:
                    best_jaccard = jaccard
                    best_compare_id = c_feat["feature_id"]
                    best_compare_route_id = c_feat["route_id"]
                    best_shared_km = s_km
                    best_target_only_km = t_only_km
                    best_compare_only_km = c_only_km

            best_match_rows.append(
                [
                    t_feat["feature_id"],
                    t_feat["route_id"],
                    best_compare_id,
                    best_compare_route_id,
                    f"{best_jaccard:.6f}",
                    f"{best_shared_km:.6f}",
                    f"{best_target_only_km:.6f}",
                    f"{best_compare_only_km:.6f}",
                ]
            )

        layer_summary_csv = os.path.join(out_dir, "레이어 비교_summary.csv")
        with open(layer_summary_csv, "w", encoding="utf-8-sig", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(
                [
                    "target_layer",
                    "compare_layer",
                    "target_feature_n",
                    "compare_feature_n",
                    "feature_n_diff",
                    "shared_km",
                    "target_only_km",
                    "compare_only_km",
                    "layer_jaccard",
                ]
            )
            writer.writerow(
                [
                    target_layer_name,
                    compare_layer_name,
                    target_feature_n,
                    compare_feature_n,
                    feature_n_diff,
                    f"{shared_km:.6f}",
                    f"{target_only_km:.6f}",
                    f"{compare_only_km:.6f}",
                    f"{layer_jaccard:.6f}",
                ]
            )

        best_match_csv = os.path.join(out_dir, "레이어 비교_match.csv")
        with open(best_match_csv, "w", encoding="utf-8-sig", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(
                [
                    "target_id",
                    "target_route_id",
                    "best_match_id",
                    "best_match_route_id",
                    "jaccard",
                    "shared_km",
                    "target_only_km",
                    "compare_only_km",
                ]
            )
            writer.writerows(best_match_rows)

        print("")
        print(f"saved: {layer_summary_csv}")
        print(f"saved: {best_match_csv}")
        print("")
        print("완료")
        print(f"  target features : {target_feature_n}")
        print(f"  compare features: {compare_feature_n}")
        print(f"  feature diff    : {feature_n_diff}")
        print(f"  layer jaccard   : {layer_jaccard:.6f}")

    except Exception as e:
        print(f"[ERROR] {e}")


if __name__ == "__main__":
    main()
