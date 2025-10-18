import pandas as pd
import yaml
from pathlib import Path

def main():
    """
    Analyze climate zone distributions in the Climate_Zone dataset.
    """
    # Find project root and load config
    project_root = Path(__file__).resolve().parents[2]
    config_path = project_root / "config.yaml"
    
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    # Define input path
    data_path = project_root / config['eqi_aamr_outputs']['base_dir'] / "Climate_Zone.csv"
    
    # Read the dataframe
    print("Reading Climate_Zone.csv...")
    df = pd.read_csv(data_path, dtype={'COUNTY_FIPS': str})
    
    # 检查 DOE 8 区的县数分布
    doe_8_dist = df.groupby('doe_zone_number')['COUNTY_FIPS'].nunique()
    print("DOE 8 Zone Distribution:")
    print(doe_8_dist)
    print(f"\nMinimum counties: {doe_8_dist.min()}")

    # 检查 DOE 15 类的县数分布
    doe_15_dist = df.groupby('doe_zone_code')['COUNTY_FIPS'].nunique()
    print("\nDOE Zone Code Distribution:")
    print(doe_15_dist.sort_values())

    # 检查 Köppen 5 大类
    koppen_major_dist = df.groupby('koppen_major')['COUNTY_FIPS'].nunique()
    print("\nKöppen Major Distribution:")
    print(koppen_major_dist)

    # 检查 Köppen 22 类
    koppen_full_dist = df.groupby('koppen_code')['COUNTY_FIPS'].nunique()
    print("\nKöppen Full Code Distribution:")
    print(koppen_full_dist.sort_values())

if __name__ == '__main__':
    main()