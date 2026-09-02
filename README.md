# Target-Approximating Path Analysis for Historical GIS

역사 GIS를 위한 목표값 근접 경로 분석 모델
**Target-approximating path analysis model for Historical GIS**

---

## 개요 (Overview)

본 저장소는 사료에 기재된 거리에 가장 부합하는 경로를 네트워크 위에서 직접 탐색하는 GIS 경로 분석 모델의 분석 스크립트를 배포합니다. 기존 최소 비용 경로(Least-Cost Path) 분석의 목적 함수를 비용 최소화 `min f`에서 목표값 근접 `min |f − T|`로 전환하여, 사료에 기록된 거리값 자체가 경로 탐색의 입력이 되도록 설계되었습니다.

분석은 사료 기록의 형식에 따라 다음 두 유형으로 구분됩니다.

- **지점 간 분석 (Point-to-point):** 출발지와 도착지가 모두 특정 지점인 경우
- **지점-경계 분석 (Point-to-boundary):** 도착지가 군현 경계 폴리곤인 경우 (예: 사방경계 기록)

목표값은 단일 값 `T` 외에 구간 `[beg, end]`로도 지정할 수 있어, 사료의 거리 기록이 범위로 주어지는 경우에 대응합니다.

이에 더하여 순위 분석, 편차범위 분석, 경로 간·레이어 간 길이 가중 유사도 분석을 통해 경로의 추정 신뢰도와 성격을 정량적으로 평가할 수 있습니다.

This repository implements a GIS-based path analysis model that searches a network for routes most closely matching distances recorded in historical sources. By transforming the objective function from cost minimization (`min f`) to target approximation (`min |f − T|`), recorded values themselves become inputs for the path search. Targets may be given as a single value `T` or as an interval `[beg, end]`.

---

## 인용 (Citation)

> 양정현 2026. 역사 GIS를 위한 목표값 근접 경로 분석 모델의 설계와 적용
> — 조선시대 경기 광주의 사례를 중심으로. *조선시대사학보 117*.
>
> Yang, J. 2026. Design and Application of a Target-Approximating
> Path Analysis Model for Historical GIS: A Case Study of Gwangju in Gyeonggi
> during the Joseon Dynasty. *THE CHOSON DYNASTY HISTORY ASSOCIATION 117*.

> 양정현 2026. 조선시대 포천과 영평의 도로 관계 기록 분석과 재현.
> *문화역사지리 38(2)*.
>
> Yang, J. 2026. Analysis and Reconstruction of Road Relationship Records
> of Pocheon and Yeongpyeong in the Joseon Dynasty.
> *Journal of Cultural and Historical Geography 38(2)*.

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
    │   ├── 라인_비용필드_계산.py
    │   ├── 네트워크데이터셋_DiGraph_변환_gpkg_pkl.py
    │   └── 그래프_무결성_점검.py
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

스크립트는 다음 3단계로 운영됩니다. 분석에 사용하는 GIS 데이터(군현 단위 도로 네트워크, 경로 분석 결과 등)는 디지털 역사지리 위키 역지사지의 [조선 군현도로 GIS](https://www.hisgeo.info/wiki/조선_군현도로_GIS) 문서를 통해 군현 단위로 배포되며, 데이터 구성과 인용 표기는 해당 문서와 각 군현 하위 문서를 참조하십시오.

### 1단계: 네트워크 그래프 변환 (`01_network_build/`)

```bash
# 1-1. 라인 데이터 + DEM → 네트워크 데이터셋(GPKG)
python scripts/01_network_build/라인_네트워크데이터셋_변환_gpkg.py

# 1-2. 네트워크 데이터셋(GPKG) → DiGraph 캐시(.pkl)
python scripts/01_network_build/네트워크데이터셋_DiGraph_변환_gpkg_pkl.py

# (보조) 환산계수 표본 노선 등 참조용 라인 레이어에 비용 필드 부여
python scripts/01_network_build/라인_비용필드_계산.py

# (선택) DiGraph 캐시의 무결성 점검
python scripts/01_network_build/그래프_무결성_점검.py
```

- `라인_비용필드_계산.py`는 경로 분석 네트워크 구축(1-1)과는 용도가 다른 보조 스크립트로, 라인을 분할하지 않고 원본 위상과 스키마를 유지한 채 정점 구간별 경사 계산으로 9개 비용 필드(3차원 거리·시간·에너지)를 피처 단위로 부여합니다. 리 수가 기록된 복원 노선 레이어에 3D 표면거리를 부여하여 里 환산계수를 산출하는 데 사용합니다.
- `그래프_무결성_점검.py`는 DiGraph 캐시의 약연결요소(WCC) 분포, 자투리 조각의 본체 이격 거리, self loop·영길이 간선 등 변환 부산물을 점검하며, 경로 CSV를 입력하면 경로별 기점·종점의 도달성 검사를 함께 수행합니다.

### 2단계: 경로 분석 (`02_route_analysis/`)

분석 목적에 따라 아래 6개 스크립트 중 선택하여 실행합니다.

| 분석 유형 | 지점 간 | 지점-경계 |
|---|---|---|
| 기본 (단일 최적해) | `경로분석_point_point.py` | `경로분석_point_polygon.py` |
| 순위 분석 (상위 N) | `경로분석_point_point_순위.py` | `경로분석_point_polygon_순위.py` |
| 편차범위 분석 (±margin) | `경로분석_point_point_편차범위.py` | `경로분석_point_polygon_편차범위.py` |

각 스크립트는 실행 시 분석 유형(거리/시간)과 장애물 적용 여부 등을 사용자가 선택합니다. 입력 CSV의 목표값 컬럼(`km_beg`/`km_end` 또는 `hr_beg`/`hr_end`)은 행별로 다음과 같이 분기합니다.

- `beg`만 입력 → 단일 목표값 모드 (`T = beg`, `|비용 − T|` 최소 경로)
- `beg`·`end` 모두 입력 → 구간 모드 `[beg, end]`
- `beg` 빈칸 → 목표값 미입력 (최소 비용 경로 `route_min_*`만 산출)

### 3단계: 유사도 분석 (`03_similarity/`)

2단계 출력물을 입력으로 하여 경로 간·레이어 간 유사도를 산출합니다.

| 스크립트 | 입력 | 출력 |
|---|---|---|
| `경로_간_유사도_분석_최소비용_영역제한.py` | 기본 분석 GPKG (+ 영역 제한 GPKG) | 유사도 CSV |
| `경로_간_유사도_분석_순위.py` | 순위 분석 GPKG | 유사도 CSV |
| `경로_간_유사도_분석_편차범위.py` | 편차범위 분석 GPKG | 유사도 CSV |
| `레이어_간_경로_구성_비교.py` | 두 개의 경로 GPKG | 요약 CSV + 매칭 CSV |

분석 모델의 구체적 절차와 출력 구조는 [docs/운영_지침.md](docs/운영_지침.md)를 참조하십시오.

---

## 변경 이력 (Changelog)

### v2.0 (2026)

- **경유지 경로의 구간 왕복 보정**: 경유지 포함 경로의 구간 병합 시 발생하는 왕복을 네트워크 위상(절단 간선·절단점)으로 판정하여, 우회 가능한 왕복이 포함된 후보는 배제하고 진입로가 유일하여 불가피한 왕복은 유지
- **里 환산계수의 산출**: 고정 환산계수의 일괄 적용 대신 분석 대상 지역을 통과하는 참조 노선에서 환산계수를 직접 산출하는 방식 도입. 이를 위해 원본 위상·스키마를 유지한 채 라인 레이어에 비용 필드를 부여하는 `라인_비용필드_계산.py` 추가

### v1.0 (2026)

- 최초 공개: 목표값 근접 경로 분석 모델(지점 간·지점-경계), 순위·편차범위 분석, 경로 간·레이어 간 유사도 분석 스크립트

---

## 라이선스 (License)

본 코드는 MIT License로 배포됩니다. 자세한 내용은 [`LICENSE`](LICENSE) 파일을 참조하십시오.

This code is released under the MIT License. See [`LICENSE`](LICENSE) for details.

---

## 연락처 (Contact)

양정현 (Yang, Jung hyun)
국립순천대학교 학술연구교수
yachagye@naver.com
