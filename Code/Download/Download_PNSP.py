import sys
import time
from pathlib import Path
import requests
import yaml

def load_paths() -> tuple[Path, str, str]:
    """Loads required paths and URL patterns from the config.yaml file."""
    project_root = Path(__file__).resolve().parents[2]
    config_path = project_root / "config.yaml"
    if not config_path.exists():
        sys.exit(f"ERROR: Configuration file not found at {config_path}")

    with config_path.open("r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)

    try:
        pnsp_config = cfg["data_sources"]["usgs_pnsp"]
        output_rel_path = pnsp_config["original"]
        base_url = cfg["data_sources"]["usgs_pnsp"]["download_urls"]["base"]
        filename_pattern = cfg["data_sources"]["usgs_pnsp"]["download_urls"]["filename_pattern"]
    except KeyError as e:
        sys.exit(f"ERROR: config.yaml is missing a required key: {e}")

    output_dir = (project_root / output_rel_path).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    return output_dir, base_url, filename_pattern

def download_file(url: str, destination: Path):
    """Downloads a file from a URL to a specified destination."""
    try:
        response = requests.get(url, stream=True, timeout=30)
        response.raise_for_status()  # Raise an exception for bad status codes (4xx or 5xx)

        with open(destination, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        print(f"  -> Successfully saved to {destination.name}")
        return True
    except requests.exceptions.HTTPError as e:
        print(f"  -> FAILED. HTTP Error: {e.response.status_code} for URL: {url}")
    except requests.exceptions.RequestException as e:
        print(f"  -> FAILED. Error downloading {url}: {e}")
    return False

def main():
    """Main function to download PNSP data for the specified years."""
    output_dir, base_url, filename_pattern = load_paths()
    start_year = 1992
    end_year = 2012

    print(f"Starting download of PNSP data from {start_year} to {end_year}.")
    print(f"Files will be saved in: {output_dir}\n")

    for year in range(start_year, end_year + 1):
        filename = filename_pattern.format(year=year)
        file_url = f"{base_url}{filename}"
        output_path = output_dir / filename

        print(f"Downloading data for {year}...")
        print(f"  URL: {file_url}")

        if output_path.exists():
            print(f"  -> File already exists. Skipping.")
            continue

        download_file(file_url, output_path)
        time.sleep(1)  # Be polite to the server

    print("\nDownload process completed.")

if __name__ == "__main__":
    main()

