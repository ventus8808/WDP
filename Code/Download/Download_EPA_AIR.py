from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import time
import zipfile
from pathlib import Path
from urllib import request, error

# 导入配置文件
import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent))
from config import get_data_dir, ensure_dir

# python3 Download_EPA_Air.py --start 1999 --end 2024 --filename daily_42401

# 修改这里的文件名中间部分来下载不同的EPA Air数据
# 例如: "daily_44201" 用于NO2数据, "daily_42401" 用于SO2数据, "daily_PRESS" 用于气压数据
FILE_NAME_PART = "daily_42602"

# 基础URL前缀
BASE_PREFIX = "https://aqs.epa.gov/aqsweb/airdata"


def get_effective_file_part(arg_filename: str | None) -> str:
    return (arg_filename or FILE_NAME_PART).strip()


def ensure_destination_dir(file_part: str) -> Path:
    # 使用配置文件中的相对路径
    dest_dir = get_data_dir("epa_air") / file_part
    return ensure_dir(dest_dir)


def download_file(url: str, dest_path: Path, retries: int = 3, backoff_seconds: float = 2.0) -> None:
    last_exc: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            # 使用临时文件，成功后原子替换，避免中断产生损坏文件
            tmp_path = dest_path.with_suffix(dest_path.suffix + ".part")
            with request.urlopen(url) as resp, tmp_path.open("wb") as out:
                # 流式写入，避免大文件占用内存
                chunk_size = 1024 * 256
                while True:
                    chunk = resp.read(chunk_size)
                    if not chunk:
                        break
                    out.write(chunk)
            tmp_path.replace(dest_path)
            return
        except (error.HTTPError, error.URLError) as exc:
            last_exc = exc
            if attempt < retries:
                time.sleep(backoff_seconds * attempt)
            else:
                raise


def extract_archive(archive_path: Path, dest_dir: Path) -> None:
    suffix = archive_path.suffix.lower()
    if suffix == ".zip":
        with zipfile.ZipFile(archive_path, 'r') as zip_ref:
            zip_ref.extractall(dest_dir)
        return
    if suffix == ".rar":
        # 优先使用 unar，其次 unrar
        unar = shutil.which("unar")
        if unar:
            subprocess.run([unar, "-o", str(dest_dir), str(archive_path)], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            return
        unrar = shutil.which("unrar")
        if unrar:
            subprocess.run([unrar, "x", "-o+", str(archive_path), str(dest_dir)], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            return
        raise RuntimeError("未找到可用的RAR解压工具，请安装 'unar' 或 'unrar'")
    # 其他后缀不支持
    raise RuntimeError(f"不支持的压缩格式: {suffix}")


def download_year(year: int, file_part: str, dest_dir: Path) -> None:
    # 目标文件名（不含后缀）
    basename = f"{file_part}_{year}"
    csv_path = dest_dir / f"{basename}.csv"
    zip_path = dest_dir / f"{basename}.zip"
    rar_path = dest_dir / f"{basename}.rar"

    # 如果CSV已存在则跳过
    if csv_path.exists() and csv_path.stat().st_size > 0:
        print(f"[SKIP] {year} CSV已存在: {csv_path}")
        return
    # 如果压缩包已存在（且非空），也跳过
    if zip_path.exists() and zip_path.stat().st_size > 0:
        print(f"[SKIP] {year} ZIP已存在: {zip_path}")
        return
    if rar_path.exists() and rar_path.stat().st_size > 0:
        print(f"[SKIP] {year} RAR已存在: {rar_path}")
        return

    # 依次尝试 .zip, .rar, .csv
    for ext in (".zip", ".rar", ".csv"):
        url = f"{BASE_PREFIX}/{basename}{ext}"
        out_path = dest_dir / f"{basename}{ext}"
        print(f"[GET ] {year} <- {url}")
        try:
            download_file(url, out_path)
            size_mb = out_path.stat().st_size / (1024 * 1024)
            print(f"[SAVE] {year} -> {out_path} ({size_mb:.2f} MB)")

            # 如果是压缩包则解压并删除
            if ext in (".zip", ".rar"):
                print(f"[EXT] {year} 解压缩中...")
                extract_archive(out_path, dest_dir)
                if csv_path.exists():
                    csv_size_mb = csv_path.stat().st_size / (1024 * 1024)
                    print(f"[DONE] {year} -> {csv_path} ({csv_size_mb:.2f} MB)")
                out_path.unlink(missing_ok=True)
                print(f"[DEL ] {year} 已删除压缩文件")
            return
        except error.HTTPError as http_err:
            if getattr(http_err, 'code', None) == 404:
                # 尝试下一种扩展名
                print(f"[MISS] {year} 远端不存在 {ext}，尝试其它格式")
                continue
            raise
        except Exception as exc:
            print(f"[FAIL] {year} 下载失败: {exc}")
            raise
    # 如果所有格式都失败
    print(f"[FAIL] {year} 未找到可用的下载格式 (.zip/.rar/.csv)")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "下载 AQS daily_*_YYYY.(zip|rar|csv) 到 'Data/Original/EPA Air/<name>/'，"
            "支持 --filename 覆盖文件内默认配置；如为压缩包将自动解压，CSV 直接保留。"
        )
    )
    parser.add_argument("--start", type=int, default=1999, help="起始年份（含）")
    parser.add_argument("--end", type=int, default=2024, help="结束年份（含）")
    parser.add_argument(
        "--filename",
        type=str,
        default=None,
        help=(
            "文件名中间部分，例如 daily_44201、daily_42401、daily_PRESS；"
            "如提供则覆盖脚本内置 FILE_NAME_PART。"
        ),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    if args.start > args.end:
        print("起始年份不能大于结束年份", file=sys.stderr)
        return 2

    file_part = get_effective_file_part(args.filename)
    dest_dir = ensure_destination_dir(file_part)
    print(f"下载目录: {dest_dir}")

    for year in range(args.start, args.end + 1):
        try:
            download_year(year, file_part, dest_dir)
        except Exception as exc:  # noqa: BLE001 - 顶层输出错误信息
            print(f"[FAIL] {year}: {exc}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
