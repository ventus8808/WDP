import ee
from typing import List


# 常量
DRIVE_FOLDER: str = "WONDER"
M2_TO_KM2: float = 1e-6

# 依据项目其它数据的年份范围（例如 CDC/JRC 1999-2020）
YEARS: List[int] = list(range(1999, 2021))

# NLCD 2019 发布对应年份集
NLCD_YEARS: List[int] = [2001, 2004, 2006, 2008, 2011, 2013, 2016, 2019]


def initialize_ee(project_hint: str = "") -> None:
    """初始化 Earth Engine，支持传入 Cloud Project。

    优先使用传入的 project_hint（例如 'nlcd-469307'），否则尝试默认项目。
    若本机未完成认证，请先运行 `earthengine authenticate`。
    """
    try:
        if project_hint:
            ee.Initialize(project=project_hint)
        else:
            ee.Initialize()
    except Exception as e:
        raise RuntimeError(
            (
                "无法初始化 Earth Engine。"
                "若提示缺少项目，请传入有效的 Cloud Project ID；"
                "若提示认证问题，请先运行 `earthengine authenticate` 完成登录。"
            )
        ) from e


def round4(x):
    """将 ee.Number 保留 4 位小数（服务端表达式）。"""
    return ee.Number(x).multiply(10000).round().divide(10000)


def get_counties(add_area: bool = False) -> ee.FeatureCollection:
    """获取美国县级边界（TIGER/2018/Counties 或 2019）并可选添加面积(km^2)。"""
    # 优先 2019，若不可用回退 2018
    try:
        fc = ee.FeatureCollection('TIGER/2019/Counties')
    except Exception:
        fc = ee.FeatureCollection('TIGER/2018/Counties')

    fc = fc.select(['GEOID'])  # 仅保留 GEOID 与几何

    if not add_area:
        return fc

    def add_area_km2(feature):
        area_km2 = ee.Number(feature.geometry().area(maxError=1)).multiply(M2_TO_KM2)
        return feature.set({'total_area_km2': round4(area_km2)})

    return fc.map(add_area_km2)


