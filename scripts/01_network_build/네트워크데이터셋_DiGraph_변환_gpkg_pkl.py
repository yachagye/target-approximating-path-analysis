# 네트워크데이터셋 .gpkg → .DiGraph 변환
# 방향 비대칭 비용 반영, pickle 캐시 저장
# 네트워크가 바뀌지 않는 한 한 번만 실행

import os
import pickle
import math
import tkinter as tk
from tkinter import filedialog

import fiona
from shapely.geometry import shape, LineString
import networkx as nx

# 노드 좌표 라운딩(5179 meter). Integrate를 했더라도 미세 부동소수점 오차 정리.
NODE_ROUND_M = 0.01  # 1 cm


def _round_xy(x, y, r=NODE_ROUND_M):
    return (round(x / r) * r, round(y / r) * r)


def _pick_gpkg():
    root = tk.Tk()
    root.withdraw()
    gpkg_path = filedialog.askopenfilename(
        title="네트워크데이터셋 GPKG 선택 (EPSG:5179, LineString)",
        filetypes=[("GeoPackage", "*.gpkg")],
    )
    if not gpkg_path:
        raise RuntimeError("GPKG 파일을 선택하지 않았습니다.")
    return gpkg_path


def _choose_layer(gpkg_path):
    layers = fiona.listlayers(gpkg_path)
    if not layers:
        raise RuntimeError("GPKG에서 레이어를 찾지 못했습니다.")
    if len(layers) == 1:
        return layers[0]
    print("GPKG 레이어 목록:")
    for i, lyr in enumerate(layers, start=1):
        print(f"  {i}. {lyr}")
    s = input("사용할 레이어 번호 입력(기본=1): ").strip()
    if not s:
        return layers[0]
    idx = int(s)
    if idx < 1 or idx > len(layers):
        raise RuntimeError("레이어 번호가 유효하지 않습니다.")
    return layers[idx - 1]


def build_graph_from_gpkg(gpkg_path, layer_name):
    # 필수 필드
    required = [
        "length_3dkm",
        "hr_ks_f", "hr_ks_b",
        "hr_tob_f", "hr_tob_b",
        "kcal_ks_f", "kcal_ks_b",
        "kcal_tob_f", "kcal_tob_b",
    ]

    G = nx.DiGraph()
    kept = 0
    dropped = 0

    with fiona.open(gpkg_path, layer=layer_name) as src:
        # CRS 확인
        crs = src.crs
        if not crs:
            raise RuntimeError("입력 GPKG 레이어 CRS가 없습니다. EPSG:5179이어야 합니다.")
        # fiona CRS가 epsg를 직접 제공하지 않을 수 있어, 최소한 proj 문자열/epsg 포함 여부 검사
        crs_str = str(crs).lower()
        if ("5179" not in crs_str) and ("epsg:5179" not in crs_str):
            # 완전 엄격 EPSG 판정은 geopandas/pyproj 경유가 필요하나,
            # 여기서는 5179가 아니면 중단(사용자 데이터 확인 유도)
            raise RuntimeError(f"입력 GPKG CRS가 EPSG:5179로 확인되지 않습니다: {crs}")

        # 필드 확인
        props = src.schema.get("properties", {})
        missing = [c for c in required if c not in props]
        if missing:
            raise RuntimeError(f"필수 필드 누락: {missing}")

        for i, feat in enumerate(src, start=1):
            try:
                geom = feat.get("geometry", None)
                if geom is None:
                    dropped += 1
                    continue

                g = shape(geom)
                if g.is_empty:
                    dropped += 1
                    continue

                # LineString / MultiLineString 처리
                if g.geom_type == "LineString":
                    line = g
                elif g.geom_type == "MultiLineString":
                    # 세그먼트 데이터라면 일반적으로 1개 구성. 첫 라인만 사용.
                    line = list(g.geoms)[0]
                else:
                    dropped += 1
                    continue

                coords = list(line.coords)
                if len(coords) < 2:
                    dropped += 1
                    continue

                x1, y1 = coords[0][0], coords[0][1]
                x2, y2 = coords[-1][0], coords[-1][1]
                u = _round_xy(x1, y1)
                v = _round_xy(x2, y2)

                p = feat.get("properties", {})

                length_km = float(p["length_3dkm"])
                if (not math.isfinite(length_km)) or length_km <= 0:
                    dropped += 1
                    continue

                # 정방향(From->To): *_f, 역방향(To->From): *_b
                attrs_f = {
                    "length_3dkm": length_km,
                    "hour_ks": float(p["hr_ks_f"]),
                    "hour_tob": float(p["hr_tob_f"]),
                    "kcal_ks": float(p["kcal_ks_f"]),
                    "kcal_tob": float(p["kcal_tob_f"]),
                    "geom": LineString(coords),
                }
                attrs_b = {
                    "length_3dkm": length_km,
                    "hour_ks": float(p["hr_ks_b"]),
                    "hour_tob": float(p["hr_tob_b"]),
                    "kcal_ks": float(p["kcal_ks_b"]),
                    "kcal_tob": float(p["kcal_tob_b"]),
                    "geom": LineString(coords[::-1]),
                }

                # DiGraph는 (u,v) 중복 시 1개만 유지 → 더 짧은 length 유지
                def _add_or_replace(a, b, attrs):
                    if G.has_edge(a, b):
                        if attrs["length_3dkm"] < G[a][b]["length_3dkm"]:
                            G[a][b].update(attrs)
                    else:
                        G.add_edge(a, b, **attrs)

                _add_or_replace(u, v, attrs_f)
                _add_or_replace(v, u, attrs_b)

                kept += 1

                if i % 200000 == 0:
                    print(f"[{i}] processed... kept={kept}, dropped={dropped}, nodes={G.number_of_nodes()}, edges={G.number_of_edges()}")

            except Exception:
                dropped += 1
                continue

    print(f"done. kept={kept}, dropped={dropped}, nodes={G.number_of_nodes()}, edges={G.number_of_edges()}")
    return G


def contract_degree2(G):
    """degree-2 체인을 축약하여 그래프를 압축한다.

    DiGraph에서 in_degree=2, out_degree=2이고 predecessors==successors인 노드를
    제거하고, 인접 간선을 병합한다. 비용은 합산, geometry는 연결한다.
    축약 간선에는 chain_nodes(원본 노드 순서열)와 chain_costs(구간별 비용)를
    저장하여, 분석 시 간선 분할에 사용할 수 있도록 한다.
    """

    # 1. 축약 대상 노드 식별
    contractible = set()
    for node in G.nodes():
        preds = set(G.predecessors(node))
        succs = set(G.successors(node))
        if len(preds) == 2 and len(succs) == 2 and preds == succs:
            contractible.add(node)

    n_total = G.number_of_nodes()
    n_contract = len(contractible)
    print(f"[contract] degree-2 nodes: {n_contract:,} / {n_total:,} ({100*n_contract/n_total:.1f}%)")

    # 2. 축약 그래프 생성 (junction/endpoint 노드만 보존)
    G_c = nx.DiGraph()
    junction_nodes = [n for n in G.nodes() if n not in contractible]
    for n in junction_nodes:
        G_c.add_node(n)

    # 3. 각 junction 노드의 출력 간선에서 체인을 따라가며 축약 간선 생성
    parallel_skipped = 0

    for start in junction_nodes:
        for first_next in G.successors(start):
            # 체인 추적: start → first_next → ... → end(junction)
            chain = [start]
            current = first_next

            while current in contractible:
                chain.append(current)
                prev = chain[-2]
                succs = [s for s in G.successors(current) if s != prev]
                if len(succs) != 1:
                    break
                current = succs[0]

            chain.append(current)
            end = current

            # 비용 합산 및 geometry 병합
            cost_keys = ("length_3dkm", "hour_ks", "hour_tob", "kcal_ks", "kcal_tob")
            total = {k: 0.0 for k in cost_keys}
            seg_costs = []
            geom_coords = []

            for i in range(len(chain) - 1):
                u, v = chain[i], chain[i + 1]
                ed = G[u][v]
                sc = {k: float(ed[k]) for k in cost_keys}
                seg_costs.append(sc)
                for k in cost_keys:
                    total[k] += sc[k]

                seg_geom = ed.get("geom")
                if seg_geom and not seg_geom.is_empty:
                    coords = list(seg_geom.coords)
                else:
                    coords = [u, v]

                if not geom_coords:
                    geom_coords.extend(coords)
                elif geom_coords[-1] == coords[0]:
                    geom_coords.extend(coords[1:])
                else:
                    geom_coords.extend(coords)

            attrs = {
                **total,
                "geom": LineString(geom_coords),
                "chain_nodes": chain,
                "chain_costs": seg_costs,
            }

            if G_c.has_edge(start, end):
                if total["length_3dkm"] < G_c[start][end]["length_3dkm"]:
                    G_c[start][end].update(attrs)
                parallel_skipped += 1
            else:
                G_c.add_edge(start, end, **attrs)

    if parallel_skipped > 0:
        print(f"[contract] WARNING: parallel chains={parallel_skipped} (shorter kept)")

    print(f"[contract] result: nodes={G_c.number_of_nodes():,}, edges={G_c.number_of_edges():,}")
    return G_c


def main():
    try:
        gpkg_path = _pick_gpkg()
        layer_name = _choose_layer(gpkg_path)

        G = build_graph_from_gpkg(gpkg_path, layer_name)
        G_c = contract_degree2(G)

        root = tk.Tk()
        root.withdraw()
        out_dir = filedialog.askdirectory(title="그래프 캐시 저장 폴더 선택")
        if not out_dir:
            raise RuntimeError("저장 폴더를 선택하지 않았습니다.")

        cache_path = os.path.join(out_dir, "graph_cache.pkl")
        payload = {
            "graph": G_c,
            "crs_epsg": 5179,
            "node_round_m": NODE_ROUND_M,
            "source_gpkg": gpkg_path,
            "source_layer": layer_name,
            "contracted": True,
            "original_nodes": G.number_of_nodes(),
            "original_edges": G.number_of_edges(),
        }

        with open(cache_path, "wb") as f:
            pickle.dump(payload, f, protocol=pickle.HIGHEST_PROTOCOL)

        print(f"saved: {cache_path}")

    except Exception as e:
        print(f"[ERROR] {e}")


if __name__ == "__main__":
    main()
