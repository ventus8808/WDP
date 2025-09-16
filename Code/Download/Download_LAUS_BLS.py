#!/usr/bin/env python3
"""
LAUS BLS 数据下载脚本
下载 BLS (Bureau of Labor Statistics) 的 LAUS (Local Area Unemployment Statistics) 数据
数据保存到 Data/Original/LAUS/ 目录
"""

import requests
import sys
from pathlib import Path

# 导入配置文件
sys.path.append(str(Path(__file__).parent.parent))
from config import get_data_dir, ensure_dir

def download_laus_data(years=None, force_download=False):
    """
    下载 LAUS BLS 数据
    
    Args:
        years: 年份列表，如果为 None 则使用默认年份
        force_download: 是否强制重新下载已存在的文件
    """
    # 默认年份（两位数格式）
    if years is None:
        years = [93, 8, 22]  # 对应1993, 2008, 2022
    
    # 获取目标目录
    dest_dir = get_data_dir("laus")
    ensure_dir(dest_dir)
    
    print(f"下载目录: {dest_dir}")
    
    # 添加请求头来模拟浏览器访问
    headers = {
        'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36',
        'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.7',
        'Accept-Language': 'en-US,en;q=0.9',
        'Accept-Encoding': 'gzip, deflate, br',
        'Referer': 'https://www.bls.gov/lau/',
        'Sec-Fetch-Dest': 'document',
        'Sec-Fetch-Mode': 'navigate',
        'Sec-Fetch-Site': 'same-origin',
        'Sec-Fetch-User': '?1',
        'Upgrade-Insecure-Requests': '1'
    }
    
    for year in years:
        # 格式化年份为两位数
        year_str = f"{year:02d}"
        url = f"https://www.bls.gov/lau/laucnty{year_str}.xlsx"
        
        # 构建完整的文件路径
        filename = f"laucnty{year_str}.xlsx"
        file_path = dest_dir / filename
        
        # 检查文件是否已存在
        if file_path.exists() and not force_download:
            print(f"[SKIP] {year_str} 文件已存在: {file_path}")
            continue
        
        print(f"[GET ] {year_str} <- {url}")
        
        try:
            response = requests.get(url, headers=headers)
            if response.status_code == 200:
                # 保存文件
                with open(file_path, "wb") as f:
                    f.write(response.content)
                
                # 获取文件大小
                file_size_mb = file_path.stat().st_size / (1024 * 1024)
                print(f"[SAVE] {year_str} -> {file_path} ({file_size_mb:.2f} MB)")
            else:
                print(f"[FAIL] {year_str} 下载失败 (状态码: {response.status_code})")
                
        except Exception as e:
            print(f"[ERROR] {year_str} 下载出错: {e}")

def list_existing_files():
    """列出已存在的 LAUS 文件"""
    dest_dir = get_data_dir("laus")
    if not dest_dir.exists():
        print("LAUS 目录不存在")
        return
    
    files = list(dest_dir.glob("laucnty*.xlsx"))
    if not files:
        print("没有找到 LAUS 文件")
        return
    
    print(f"已存在的 LAUS 文件 ({len(files)} 个):")
    # 按年份排序
    files.sort(key=lambda x: int(x.stem.replace("laucnty", "")))
    
    for file in files:
        year = file.stem.replace("laucnty", "")
        size_mb = file.stat().st_size / (1024 * 1024)
        print(f"  {year}: {file.name} ({size_mb:.2f} MB)")

def main():
    """主函数"""
    # 可以通过命令行参数指定年份
    if len(sys.argv) > 1:
        try:
            years = [int(year) for year in sys.argv[1:]]
            print(f"使用命令行指定的年份: {years}")
        except ValueError:
            print("错误: 年份必须是整数")
            sys.exit(1)
    else:
        years = None
        print("使用默认年份: [93, 8, 22] (对应1993, 2008, 2022)")
    
    download_laus_data(years)
    print("\n下载完成!")

if __name__ == "__main__":
    main()
