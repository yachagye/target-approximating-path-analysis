# 그래프 무결성 점검 (graph_cache.pkl)
# - WCC 분포 및 자투리 조각의 끝점→본체 거리
# - self_loop, zero_len, isolated 등 변환 부산물
# - (선택) routes.csv 입력 시 경로별 기점·종점 도달성 검사

import csv
import pickle
import tkinter as tk
from tkinter import filedialog

import networkx as nx
from shapely.geometry import Point, LineString
from pyproj import Transformer


def load_graph(path):
    with open(path, "rb") as f:
        payload = pickle.load(f)
    if isinstance(payload, dict) and "graph" in payload:
        return payload["graph"]
    return payload


def build_main_geoms(G, main_set):
    """본체 간선 geometry 목록 (자투리 끝점 거리 측정용)."""
    geoms = []
    for u, v, data in G.subgraph(main_set).edges(data=True):
        geom = data.get("geom")
        if geom is None or geom.is_empty:
            geom = LineString([u, v])
        geoms.append(geom)
    return geoms


def check_components(G):
    comps = sorted(nx.weakly_connected_components(G), key=len, reverse=True)
    main_set = comps[0]
    print(f"[WCC] 개수={len(comps)}  (본체 {len(main_set):,} nodes)")

    if len(comps) > 1:
        main_geoms = build_main_geoms(G, main_set)
        for i, comp in enumerate(comps[1:], start=1):
            sub = G.subgraph(comp)
            glen = 0.0
            for u, v, data in sub.edges(data=True):
                geom = data.get("geom")
                if geom is not None and not geom.is_empty:
                    glen = max(glen, geom.length)

            # 자투리 끝점(degree-1)이 본체에서 떨어진 거리 — 단절 지점은 끝점에서 발생
            endpoints = [n for n in comp if sub.degree(n) <= 1] or list(comp)
            end_dists = []
            for n in endpoints:
                pt = Point(n[0], n[1])
                end_dists.append(min(g.distance(pt) for g in main_geoms))
            d_min, d_max = min(end_dists), max(end_dists)

            # 끝점 중 하나라도 본체에 닿아 있으면(0.1m 이내) 노드 병합 실패에 의한 단절
            if d_min <= 0.1:
                verdict = "★ 단절 (끝점이 본체에 닿았으나 노드 미연결 → 변환 병합 실패)"
            elif d_min <= 5.0:
                verdict = f"△ 근접 단절 의심 (끝점이 본체에서 {d_min:.1f}m)"
            else:
                verdict = "별개 도로 (본체와 떨어짐)"
            print(f"  자투리 WCC {i}: {len(comp)} nodes, 최대간선={glen:.1f}m, "
                  f"끝점→본체 최소={d_min:.1f}m / 최대={d_max:.1f}m → {verdict}")

    return main_set


def check_edge_anomalies(G):
    self_loops = list(nx.selfloop_edges(G))
    zero_len = 0
    for u, v, data in G.edges(data=True):
        geom = data.get("geom")
        if geom is None or geom.is_empty or geom.length < 1e-6:
            zero_len += 1
    isolated = [n for n in G.nodes() if G.degree(n) == 0]

    print(f"[anomaly] self_loop={len(self_loops)}, zero_len_edge={zero_len}, "
          f"isolated_node={len(isolated)}")
    for u, v in self_loops:
        print(f"  self_loop: {u}")


def parse_coords(parts):
    nums = [p.strip() for p in parts[1:]]
    coords = []
    for j in range(0, len(nums) - 1, 2):
        a, b = nums[j], nums[j + 1]
        if a and b:
            try:
                coords.append((float(a), float(b)))
            except ValueError:
                break
    return coords


def check_routes(G, routes_csv, main_set):
    """좌표를 최근접 노드로 매칭하여 기점·종점 도달성을 검사한다.
    (분석의 간선 투영 스냅과 미세 차이는 있으나, 경로 단절 판정에는 충분)"""
    nodes_list = list(G.nodes())
    tf = Transformer.from_crs(4326, 5179, always_xy=True)

    def nearest_node(x, y):
        best, best_d2 = None, float("inf")
        for n in nodes_list:
            d2 = (n[0] - x) ** 2 + (n[1] - y) ** 2
            if d2 < best_d2:
                best_d2, best = d2, n
        return best, best_d2 ** 0.5

    n_total = n_fail = 0
    with open(routes_csv, "r", encoding="utf-8-sig") as f:
        reader = csv.reader(f)
        next(reader, None)
        for parts in reader:
            if not parts or not parts[0].strip():
                continue
            rid = parts[0].strip()
            coords = parse_coords(parts)
            if len(coords) < 2:
                print(f"[{rid}] 좌표 부족({len(coords)})")
                continue

            n_total += 1
            (x1, y1), (x2, y2) = coords[0], coords[-1]
            sx1, sy1 = tf.transform(x1, y1)
            sx2, sy2 = tf.transform(x2, y2)
            n1, d1 = nearest_node(sx1, sy1)
            n2, d2 = nearest_node(sx2, sy2)

            if nx.has_path(G, n1, n2):
                print(f"[{rid}] OK (기점 snap≈{d1:.1f}m, 종점 snap≈{d2:.1f}m)")
            else:
                n_fail += 1
                bad = []
                if n1 not in main_set:
                    bad.append("기점 본체밖")
                if n2 not in main_set:
                    bad.append("종점 본체밖")
                if not bad:
                    bad.append("방향 도달불가")
                print(f"[{rid}] ★ 경로불가 (기점 snap≈{d1:.1f}m, 종점 snap≈{d2:.1f}m) "
                      f"→ {', '.join(bad)}")

    print(f"\n[경로 점검 요약] 총 {n_total}건 중 실패 {n_fail}건")


def main():
    try:
        root = tk.Tk()
        root.withdraw()

        cache_path = filedialog.askopenfilename(
            title="graph_cache.pkl 선택", filetypes=[("Pickle", "*.pkl")])
        if not cache_path:
            raise RuntimeError("graph_cache.pkl을 선택하지 않았습니다.")

        routes_csv = filedialog.askopenfilename(
            title="(선택) routes.csv — 취소 시 그래프만 점검", filetypes=[("CSV", "*.csv")])

        G = load_graph(cache_path)
        print(f"[graph] nodes={G.number_of_nodes():,}, edges={G.number_of_edges():,}")

        main_set = check_components(G)
        check_edge_anomalies(G)

        if routes_csv:
            print()
            check_routes(G, routes_csv, main_set)

    except Exception as e:
        print(f"[ERROR] {e}")


if __name__ == "__main__":
    main()