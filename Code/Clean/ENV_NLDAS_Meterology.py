#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NLDAS气象数据县级聚合脚本
将0.125° NLDAS栅格数据聚合到县级，生成年度和季节气象指标
"""

import os
import sys
from pathlib import Path
import glob
import numpy as np
import pandas as pd
import xarray as xr
import geopandas as gpd
from shapely.geometry import Point, Polygon
from scipy.sparse import csr_matrix
import warnings
import yaml
warnings.filterwarnings('ignore')

# ==================== Configuration ====================
START_YEAR = 1999
END_YEAR = 2019
MAX_COUNTIES_FOR_TEST = None
# =================================================

def load_paths():
    """Loads paths from config.yaml."""
    project_root = Path(__file__).resolve().parents[2]
    config_path = project_root / "config.yaml"
    if not config_path.exists():
        sys.exit(f"ERROR: Configuration file not found: {config_path}")

    with config_path.open("r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)

    try:
        nldas_config = cfg["data_sources"]["nldas"]
        nldas_folder = (project_root / nldas_config["original"]).resolve()
        output_folder = (project_root / nldas_config["processed"]).resolve()

        tiger_config = cfg["data_sources"]["tiger"]
        county_shape_file = (project_root / tiger_config["shapefile"]).resolve()

    except KeyError as e:
        sys.exit(f"ERROR: config.yaml is missing a required path: {e}")

    if not nldas_folder.exists():
        sys.exit(f"ERROR: NLDAS input directory not found: {nldas_folder}")
    if not county_shape_file.exists():
        sys.exit(f"ERROR: County shapefile not found: {county_shape_file}")

    output_folder.mkdir(parents=True, exist_ok=True)
    return nldas_folder, county_shape_file, output_folder

def get_nldas_files(nldas_folder, start_year, end_year):
    """获取指定年份范围的NLDAS文件；自动补充起始年上一年12月用于DJF"""
    print(f"查找 {start_year}-{end_year} 年的NLDAS文件...")

    all_files = []
    # 先尝试加入起始年前一年的12月（用于DJF）
    if start_year - 1 >= 0:
        pattern_prev_dec = f"NLDAS_FORA0125_M.A{start_year-1:04d}12.020.nc*"
        files_prev_dec = glob.glob(os.path.join(nldas_folder, pattern_prev_dec))
        if files_prev_dec:
            all_files.append((start_year - 1, 12, files_prev_dec[0]))
        else:
            print("未找到起始年前一年的12月数据，DJF(起始年)将缺少12月。")

    for year in range(start_year, end_year + 1):
        for month in range(1, 13):
            pattern = f"NLDAS_FORA0125_M.A{year:04d}{month:02d}.020.nc*"
            files = glob.glob(os.path.join(nldas_folder, pattern))
            if files:
                all_files.append((year, month, files[0]))
            else:
                print(f"未找到 {year:04d}-{month:02d} 的数据")

    print(f"找到 {len(all_files)} 个文件")
    return all_files

def load_county_boundaries(shape_file):
    """加载县级边界数据"""
    print(f"加载县级边界数据: {shape_file}")

    gdf = gpd.read_file(shape_file)

    # 限制县的数量（如果设置了MAX_COUNTIES_FOR_TEST）
    if MAX_COUNTIES_FOR_TEST is not None:
        print(f"仅处理前 {MAX_COUNTIES_FOR_TEST} 个县进行测试。")
        gdf = gdf.head(MAX_COUNTIES_FOR_TEST)

    print(f"加载了 {len(gdf)} 个县")

    # 转换到等积投影用于面积计算
    gdf_albers = gdf.to_crs('EPSG:5070')

    return gdf, gdf_albers

def create_spatial_weight_matrix(gdf_albers, nldas_sample_file):
    """创建空间权重矩阵"""
    print("构建空间权重矩阵...")

    # 加载NLDAS样本文件获取网格信息
    with xr.open_dataset(nldas_sample_file, engine='netcdf4') as ds:
        lons = ds.lon.values
        lats = ds.lat.values

    print(f"  NLDAS网格: {len(lons)} x {len(lats)} = {len(lons) * len(lats)} 个网格点")

    # 创建网格多边形
    print("  创建网格多边形...")
    grid_polygons = []
    grid_indices = []

    for i, lat in enumerate(lats):
        for j, lon in enumerate(lons):
            grid_poly = Polygon([
                (lon - 0.0625, lat - 0.0625),
                (lon + 0.0625, lat - 0.0625),
                (lon + 0.0625, lat + 0.0625),
                (lon - 0.0625, lat + 0.0625)
            ])
            grid_polygons.append(grid_poly)
            grid_indices.append((i, j))

    print(f"  创建了 {len(grid_polygons)} 个网格多边形")

    # 转换网格到等积投影
    print("  转换网格到等积投影...")
    grid_gdf = gpd.GeoDataFrame(
        geometry=grid_polygons,
        crs='EPSG:4326'
    ).to_crs('EPSG:5070')

    # 计算县-网格交叠面积
    print("  计算县-网格交叠面积...")
    weight_matrix = np.zeros((len(gdf_albers), len(grid_polygons)))

    for county_idx, county_geom in enumerate(gdf_albers.geometry):
        if county_idx % 50 == 0:
            print(f"    [{county_idx+1:4d}/{len(gdf_albers)}] 处理县 {county_idx + 1}/{len(gdf_albers)}")

        # 检查县是否有有效几何
        if county_geom is None or county_geom.is_empty:
            continue

        try:
            # 计算与所有网格的交叠
            for grid_idx, grid_geom in enumerate(grid_gdf.geometry):
                if grid_geom is not None and not grid_geom.is_empty:
                    intersection = county_geom.intersection(grid_geom)
                    if not intersection.is_empty:
                        weight_matrix[county_idx, grid_idx] = intersection.area
        except Exception as e:
            print(f"    处理县 {county_idx} 时出错: {e}")
            continue

        # 归一化权重（按县面积）
        county_area = county_geom.area
        if county_area > 0:
            weight_matrix[county_idx, :] /= county_area

    print(f"  完成 {len(gdf_albers)} 个县的权重计算")

    # 转换为稀疏矩阵
    print("  转换为稀疏矩阵...")
    weight_sparse = csr_matrix(weight_matrix)

    print(f"  权重矩阵形状: {weight_sparse.shape}")
    print(f"  非零元素: {weight_sparse.nnz}")

    return weight_sparse, grid_indices

def load_nldas_monthly_data(file_path):
    """加载单个月NLDAS数据"""
    try:
        with xr.open_dataset(file_path, engine='netcdf4') as ds:
            data = {}

            # 2米气温 (K -> ℃)
            if 'Tair' in ds:
                data['tas'] = ds['Tair'].values[0] - 273.15

            # 风速 (合成U和V分量)
            if 'Wind_E' in ds and 'Wind_N' in ds:
                u = ds['Wind_E'].values[0]
                v = ds['Wind_N'].values[0]
                data['wind'] = np.sqrt(u**2 + v**2)

            # 总降水量 (kg m-2 -> mm/month)
            if 'Rainf' in ds:
                # 根据官方文档，Rainf单位是kg m-2，这是月累计值
                # 1 kg/m² = 1 mm (水的密度)
                data['prcp'] = ds['Rainf'].values[0]  # 直接使用，单位已经是mm/month

            # 比湿 (kg kg-1 -> 相对湿度%)
            if 'Qair' in ds and 'Tair' in ds and 'PSurf' in ds:
                # 从比湿计算相对湿度
                q = ds['Qair'].values[0]  # kg/kg
                t = ds['Tair'].values[0]  # K
                p = ds['PSurf'].values[0]  # Pa

                # 计算相对湿度
                t_c = t - 273.15  # 转换为摄氏度
                # 饱和水汽压 (Magnus公式)
                es = 6.112 * np.exp(17.67 * t_c / (t_c + 243.5)) * 100  # hPa -> Pa
                # 实际水汽压
                e = q * p / (0.622 + q * 0.378)
                # 相对湿度
                data['rh'] = np.clip(e / es * 100, 0, 100)

            # 短波辐射 (W m-2 -> W/m²)
            if 'SWdown' in ds:
                data['swrad'] = ds['SWdown'].values[0]

            # 长波辐射 (W m-2 -> W/m²)
            if 'LWdown' in ds:
                data['lwrad'] = ds['LWdown'].values[0]

            # 地表气压 (Pa -> kPa)
            if 'PSurf' in ds:
                data['psurf'] = ds['PSurf'].values[0] / 1000  # Pa -> kPa

            # 对流有效位能 (J kg-1 -> J/kg)
            if 'CAPE' in ds:
                data['cape'] = ds['CAPE'].values[0]

            # 潜在蒸发 (kg m-2 -> mm/month)
            if 'PotEvap' in ds:
                data['potevap'] = ds['PotEvap'].values[0]  # 直接使用，单位已经是mm/month

            return data

    except Exception as e:
        print(f"加载文件失败 {file_path}: {e}")
        return None

def aggregate_to_counties(weight_matrix, grid_data, grid_indices):
    """将网格数据聚合到县级"""
    grid_vector = np.zeros(len(grid_indices))

    for idx, (i, j) in enumerate(grid_indices):
        if i < grid_data.shape[0] and j < grid_data.shape[1]:
            value = grid_data[i, j]
            if not np.isnan(value):
                grid_vector[idx] = value

    county_values = weight_matrix.dot(grid_vector)

    # 处理空值：如果县的权重和为0，则设为NaN
    weight_sums = weight_matrix.sum(axis=1).A1
    county_values[weight_sums == 0] = np.nan

    return county_values

def process_nldas_data_by_year(nldas_files, weight_matrix, grid_indices, gdf, output_folder):
    """按年份处理NLDAS数据并聚合到县级，每年保存一次"""
    print("开始按年份处理NLDAS数据...")

    # 按年份分组文件，但需要特殊处理起始年的上一年12月数据
    files_by_year = {}
    prev_year_dec_files = []  # 存储起始年前一年的12月数据

    for year, month, file_path in nldas_files:
        # 如果是起始年前一年的12月，单独存储
        if year == START_YEAR - 1 and month == 12:
            prev_year_dec_files.append((year, month, file_path))
        else:
            if year not in files_by_year:
                files_by_year[year] = []
            files_by_year[year].append((year, month, file_path))

    print(f"将处理 {len(files_by_year)} 年的数据")
    if prev_year_dec_files:
        print(f"找到起始年前一年12月数据: {len(prev_year_dec_files)} 个文件")

    all_results = []

    for year in sorted(files_by_year.keys()):
        year_files = files_by_year[year]
        print(f"\n处理 {year} 年数据 ({len(year_files)} 个月)...")

        year_results = []

        # 如果是起始年，需要加入前一年12月的数据用于DJF计算
        if year == START_YEAR and prev_year_dec_files:
            print(f"  加入前一年12月数据用于 {year} 年DJF计算...")
            for prev_year, prev_month, prev_file_path in prev_year_dec_files:
                print(f"    [PREV] 处理 {prev_year:04d}-{prev_month:02d}: {os.path.basename(prev_file_path)[:50]}...")

                monthly_data = load_nldas_monthly_data(prev_file_path)
                if monthly_data is None:
                    print(f"      跳过 {prev_year:04d}-{prev_month:02d} 由于数据加载失败")
                    continue

                county_data = {}
                for var_name, grid_data in monthly_data.items():
                    county_values = aggregate_to_counties(weight_matrix, grid_data, grid_indices)
                    county_data[var_name] = county_values

                for county_idx in range(len(gdf)):
                    result_row = {
                        'GEOID': gdf.iloc[county_idx]['GEOID'],
                        'year': prev_year,  # 保持原始年份，在统计时会处理
                        'month': prev_month
                    }

                    for var_name, values in county_data.items():
                        result_row[f'{var_name}_{prev_month:02d}'] = values[county_idx]

                    year_results.append(result_row)

        # 处理当年的数据
        for file_idx, (y, month, file_path) in enumerate(year_files):
            print(f"  [{file_idx+1:2d}/{len(year_files)}] 处理 {y:04d}-{month:02d}: {os.path.basename(file_path)[:50]}...")

            monthly_data = load_nldas_monthly_data(file_path)
            if monthly_data is None:
                print(f"    跳过 {y:04d}-{month:02d} 由于数据加载失败")
                continue

            county_data = {}
            for var_name, grid_data in monthly_data.items():
                county_values = aggregate_to_counties(weight_matrix, grid_data, grid_indices)
                county_data[var_name] = county_values

            for county_idx in range(len(gdf)):
                result_row = {
                    'GEOID': gdf.iloc[county_idx]['GEOID'],
                    'year': y,
                    'month': month
                }

                for var_name, values in county_data.items():
                    result_row[f'{var_name}_{month:02d}'] = values[county_idx]

                year_results.append(result_row)

        # 计算该年的年度和季节统计
        print(f"  计算 {year} 年统计量...")
        year_df = pd.DataFrame(year_results)
        annual_df = calculate_annual_seasonal_stats(year_df, target_year=year)

        # 只保存起始年及之后的数据
        if year >= START_YEAR:
            # Rename GEOID to COUNTY_FIPS for consistency with other datasets
            annual_df = annual_df.rename(columns={"GEOID": "COUNTY_FIPS"})

            year_output_file = f"NLDAS_{year}.csv"
            year_output_path = os.path.join(output_folder, year_output_file)
            annual_df.to_csv(year_output_path, index=False, encoding='utf-8')

            print(f"  {year} 年数据已保存: {year_output_path}")
            print(f"  {year} 年数据形状: {annual_df.shape}")

            all_results.append(annual_df)
        else:
            print(f"  跳过 {year} 年数据保存（早于起始年）")

    # 合并所有年份的数据
    if all_results:
        final_df = pd.concat(all_results, ignore_index=True)
        print(f"\n所有年份处理完成！")
        print(f"总数据形状: {final_df.shape}")
        return final_df
    else:
        print("处理失败，没有生成数据")
        return pd.DataFrame()

def calculate_annual_seasonal_stats(df, target_year=None):
    """计算年度和季节统计"""
    print("计算年度和季节统计...")

    # 如果指定了目标年份，只处理该年份的数据
    if target_year is not None:
        df = df[df['year'] == target_year].copy()
        print(f"  只处理 {target_year} 年的数据")

    # DJF跨年处理：将上一年12月计入下一年的DJF
    # 思路：对每个GEOID，复制一份12月数据并把year加1，再合并后做分组
    df_adj = df.copy()
    is_dec = df_adj['month'] == 12
    if is_dec.any():
        print("  处理DJF跨年数据...")
        dec_rows = df_adj[is_dec].copy()
        dec_rows.loc[:, 'year'] = dec_rows['year'] + 1
        df_adj = pd.concat([df_adj, dec_rows], ignore_index=True)
        print(f"  已复制 {len(dec_rows)} 条12月数据到下一年")

    # 获取所有唯一的GEOID和year组合
    unique_combinations = df_adj[['GEOID', 'year']].drop_duplicates()
    total_combinations = len(unique_combinations)
    print(f"  开始计算 {total_combinations} 个县-年组合的统计量...")

    annual_stats = []
    processed_count = 0

    for (geoid, year), group in df_adj.groupby(['GEOID', 'year']):
        processed_count += 1
        if processed_count % 50 == 0:
            print(f"    [{processed_count:4d}/{total_combinations}] 处理 {geoid}-{year}")

        stats = {'GEOID': geoid, 'year': year}

        # 年度统计
        for var in ['tas', 'wind', 'prcp', 'rh', 'swrad', 'lwrad', 'psurf', 'cape', 'potevap']:
            monthly_cols = [col for col in group.columns if col.startswith(f'{var}_')]
            if monthly_cols:
                values = group[monthly_cols].values.flatten()
                values = values[~np.isnan(values)]

                if len(values) > 0:
                    if var in ['prcp', 'potevap']:
                        stats[f'{var}_sum_annual'] = np.sum(values)  # 累计值
                    else:
                        stats[f'{var}_mean_annual'] = np.mean(values)  # 平均值

        # 季节统计
        seasons = {
            'DJF': [12, 1, 2],  # 已在上方通过复制12月至下一年实现跨年
            'MAM': [3, 4, 5],
            'JJA': [6, 7, 8],
            'SON': [9, 10, 11]
        }

        for season, months in seasons.items():
            for var in ['tas', 'wind', 'prcp', 'rh', 'swrad', 'lwrad', 'psurf', 'cape', 'potevap']:
                monthly_cols = [col for col in group.columns if col.startswith(f'{var}_')]
                season_cols = [col for col in monthly_cols if any(f'{m:02d}' in col for m in months)]

                if season_cols:
                    values = group[season_cols].values.flatten()
                    values = values[~np.isnan(values)]

                    if len(values) > 0:
                        if var in ['prcp', 'potevap']:
                            stats[f'{var}_sum_{season}'] = np.sum(values)  # 累计值
                        else:
                            stats[f'{var}_mean_{season}'] = np.mean(values)  # 平均值

        annual_stats.append(stats)

    print(f"  统计计算完成，共处理 {processed_count} 个县-年组合")

    out = pd.DataFrame(annual_stats)

    # 如果指定了目标年份，只返回该年份的数据
    if target_year is not None:
        out = out[out['year'] == target_year]
        print(f"  只输出 {target_year} 年的数据，共 {len(out)} 条记录")
    else:
        # 否则确保只输出起始年到结束年范围
        out = out[(out['year'] >= START_YEAR) & (out['year'] <= END_YEAR)]
        print(f"  最终输出 {len(out)} 条记录")

    return out

def main():
    """Main execution function."""
    print("NLDAS Meteorology Data Aggregation to County Level")
    print("="*60)

    # Load paths from config
    nldas_folder, county_shape_file, output_folder = load_paths()
    output_file = f"NLDAS_{START_YEAR}_{END_YEAR}.csv"

    # 1. Get NLDAS file list
    print("Step 1: Get NLDAS file list")
    nldas_files = get_nldas_files(nldas_folder, START_YEAR, END_YEAR)
    if not nldas_files:
        print("No NLDAS files found!")
        return

    # 2. Load county boundaries
    print("\nStep 2: Load county boundaries")
    gdf, gdf_albers = load_county_boundaries(county_shape_file)

    # 3. Build or load spatial weight matrix
    print("\nStep 3: Build or load spatial weight matrix")
    weight_matrix_path = output_folder / "county_grid_weights.npz"
    if weight_matrix_path.exists():
        print("Loading existing weight matrix...")
        weight_data = np.load(weight_matrix_path)
        weight_matrix = csr_matrix((weight_data['data'], weight_data['indices'], weight_data['indptr']), shape=weight_data['shape'])
        grid_indices = [tuple(row) for row in weight_data['grid_indices']]
    else:
        print("Creating new weight matrix...")
        weight_matrix, grid_indices = create_spatial_weight_matrix(gdf_albers, nldas_files[0][2])
        print("Saving weight matrix...")
        np.savez(weight_matrix_path, data=weight_matrix.data, indices=weight_matrix.indices, indptr=weight_matrix.indptr, shape=weight_matrix.shape, grid_indices=np.array(grid_indices))

    # 4. Process NLDAS data by year
    print("\nStep 4: Process NLDAS data by year")
    annual_df = process_nldas_data_by_year(nldas_files, weight_matrix, grid_indices, gdf, output_folder)

    # 5. Save combined results
    if not annual_df.empty:
        print("\nStep 5: Save combined results")

        # Rename GEOID to COUNTY_FIPS for consistency with other datasets
        annual_df = annual_df.rename(columns={"GEOID": "COUNTY_FIPS"})

        output_path = output_folder / output_file
        annual_df.to_csv(output_path, index=False, encoding='utf-8')
        print(f"  Combined data saved to: {output_path}")

    if not annual_df.empty:
        print(f"\n处理完成！")
        print(f"结果数据形状: {annual_df.shape}")
        print(f"年份范围: {annual_df['year'].min()} - {annual_df['year'].max()}")
        print(f"县数量: {annual_df['COUNTY_FIPS'].nunique()}")

        # 显示结果示例
        print(f"\n数据示例:")
        print(annual_df.head())

        # 显示变量列表
        print("\n生成的变量:")
        var_cols = [col for col in annual_df.columns if col not in ['COUNTY_FIPS', 'year']]
        for var in sorted(var_cols):
            print(f"  {var}")

        print(f"\n所有处理完成！共生成 {len(var_cols)} 个气象变量")

        # 显示年度文件列表
        print("\n生成的年度文件:")
        import glob
        year_files = glob.glob(os.path.join(output_folder, "NLDAS_*.csv"))
        for year_file in sorted(year_files):
            file_size = os.path.getsize(year_file) / (1024 * 1024)  # MB
            print(f"  {os.path.basename(year_file)} ({file_size:.1f} MB)")
    else:
        print("\n处理失败，没有生成数据")

if __name__ == "__main__":
    main()
