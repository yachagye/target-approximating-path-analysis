# Target-Approximating Path Analysis for Historical GIS

조선시대 경기 광주의 사례를 중심으로 한 역사 GIS 경로 분석 모델
**Target-approximating path analysis model for Historical GIS — A case study of Gwangju in Gyeonggi during the Joseon Dynasty**

---

## 개요 (Overview)

본 저장소는 사료에 기재된 거리에 가장 부합하는 경로를 네트워크 위에서 직접 탐색하는 GIS 경로 분석 모델의 분석 스크립트를 배포합니다. 기존 최소 비용 경로(Least-Cost Path) 분석의 목적 함수를 비용 최소화 `min f`에서 목표값 근접 `min |f − T|`로 전환하여, 사료에 기록된 거리값 자체가 경로 탐색의 입력이 되도록 설계되었습니다.

분석은 사료 기록의 형식에 따라 다음 두 유형으로 구분됩니다.

- **지점 간 분석 (Point-to-point):** 출발지와 도착지가 모두 특정 지점인 경우
- **지점-경계 분석 (Point-to-boundary):** 도착지가 군현 경계 폴리곤인 경우 (예: 사방경계 기록)

이 외에도 순위 분석, 편차범위 분석, 경로 간·레이어 간 길이 가중 유사도 분석을 포함하여 경로의 추정 신뢰도와 성격을 정량적으로 평가합니다.

This repository implements a GIS-based path analysis model that searches a network for routes most closely matching distances recorded in historical sources. By transforming the objective function from cost minimization (`min f`) to target approximation (`min |f − T|`), recorded values themselves become inputs for the path search.

---

## 인용 (Citation)

> 양정현 (). 역사 GIS를 위한 목표값 근접 경로 분석 모델의 설계와 적용
> — 조선시대 경기 광주의 사례를 중심으로. *(투고 학술지명)*.
>
> Yang, J. (). Design and Application of a Target-Approximating
> Path Analysis Model for Historical GIS: A Case Study of Gwangju in Gyeonggi
> during the Joseon Dynasty. *(Journal name)*.

*※ 게재 후 권·호·페이지·DOI를 갱신할 예정*

---

## 저장소 구성

```
target-approximating-path-analysis/
├── README.md
├── LICENSE
├── requirements.txt
├── .gitignore
├── docs/
│   └── 운영_지침.md           # 모델·알고리즘 상세 명세
└── scripts/
    ├── 01_network_build/      # 네트워크 데이터셋 구축
    │   ├── 라인_네트워크데이터셋_변환_gpkg.py
    │   └── 네트워크데이터셋_DiGraph_변환_gpkg_pkl.py
    ├── 02_route_analysis/     # 경로 분석 (기본·순위·편차범위)
    │   ├── 경로분석_point_point.py
    │   ├── 경로분석_point_polygon.py
    │   ├── 경로분석_point_point_순위.py
    │   ├── 경로분석_point_polygon_순위.py
    │   ├── 경로분석_point_point_편차범위.py
    │   └── 경로분석_point_polygon_편차범위.py
    └── 03_similarity/         # 경로·레이어 유사도 분석
        ├── 경로_간_유사도_분석_최소비용_영역제한.py
        ├── 경로_간_유사도_분석_순위.py
        ├── 경로_간_유사도_분석_편차범위.py
        └── 레이어_간_경로_구성_비교.py
```

---

## 실행 환경

- Python 3.10 이상
- 주요 라이브러리: NetworkX, GeoPandas, Shapely, pyproj, fiona, NumPy, rasterio
- 운영 체제: Windows / macOS / Linux 모두 지원

설치는 다음과 같이 수행합니다.

```bash
pip install -r requirements.txt
```

ArcGIS Pro·QGIS 등 상용·전용 GIS 소프트웨어 없이 Python 오픈소스 라이브러리만으로 전체 파이프라인이 실행됩니다.

---

## 실행 순서

스크립트는 다음 3단계로 운영됩니다. 입력 데이터(GIS 도로망, DEM, 군현 폴리곤, 사료 좌표 CSV 등)는 별도로 배포되며, 배포처는 추후 README에 갱신할 예정입니다.

### 1단계: 네트워크 그래프 변환 (`01_network_build/`)

```bash
# 1-1. 라인 데이터 + DEM → 네트워크 데이터셋(GPKG)
python scripts/01_network_build/라인_네트워크데이터셋_변환_gpkg.py

# 1-2. 네트워크 데이터셋(GPKG) → DiGraph 캐시(.pkl)
python scripts/01_network_build/네트워크데이터셋_DiGraph_변환_gpkg_pkl.py
```

### 2단계: 경로 분석 (`02_route_analysis/`)

분석 목적에 따라 아래 6개 스크립트 중 선택하여 실행합니다.

| 분석 유형 | 지점 간 | 지점-경계 |
|---|---|---|
| 기본 (단일 최적해) | `경로분석_point_point.py` | `경로분석_point_polygon.py` |
| 순위 분석 (상위 N) | `경로분석_point_point_순위.py` | `경로분석_point_polygon_순위.py` |
| 편차범위 분석 (±margin) | `경로분석_point_point_편차범위.py` | `경로분석_point_polygon_편차범위.py` |

각 스크립트는 실행 시점에 분석 유형(거리/시간), 장애물 적용 여부 등을 대화식으로 선택합니다.

### 3단계: 유사도 분석 (`03_similarity/`)

2단계 출력물을 입력으로 하여 경로 간·레이어 간 유사도를 산출합니다.

| 스크립트 | 입력 | 출력 |
|---|---|---|
| `경로_간_유사도_분석_최소비용_영역제한.py` | 기본 분석 GPKG (+ 영역 제한 GPKG) | 유사도 CSV |
| `경로_간_유사도_분석_순위.py` | 순위 분석 GPKG | 유사도 CSV + 필터 GPKG |
| `경로_간_유사도_분석_편차범위.py` | 편차범위 분석 GPKG | 유사도 CSV + 필터 GPKG |
| `레이어_간_경로_구성_비교.py` | 두 개의 경로 GPKG | 요약 CSV + 매칭 CSV |

상세한 알고리즘 명세와 입출력 구조는 [`docs/운영_지침.md`](docs/운영_지침.md)를 참조하십시오.

---

## 라이선스 (License)

본 코드는 MIT License로 배포됩니다. 자세한 내용은 [`LICENSE`](LICENSE) 파일을 참조하십시오.

This code is released under the MIT License. See [`LICENSE`](LICENSE) for details.

---

## 연락처 (Contact)

양정현 (Yang, Jung hyun)
국립순천대학교 학술연구교수
yachagye@naver.com
