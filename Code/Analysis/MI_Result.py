import pandas as pd
import logging

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class MIResultProcessor:
    """
    处理多重插补结果的类
    """
    
    def __init__(self):
        """
        初始化处理器
        """
        logger.info("初始化MI结果处理器")
    
    def load_mice_data(self, file_path: str) -> pd.DataFrame:
        """
        读取MICE插补数据
        
        参数:
            file_path (str): 插补数据文件路径
            
        返回:
            DataFrame: 插补数据
        """
        logger.info(f"读取MICE数据: {file_path}")
        df = pd.read_csv(file_path)
        logger.info(f"数据形状: {df.shape}")
        return df
    
    def convert_to_long_format(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        将数据转换为正确格式（修复COUNTY_FIPS和EQI_Period格式，并按正确顺序排列癌症类型）
        
        返回:
            修复后的数据框
        """
        logger.info("=== 修复数据格式并按正确顺序排列癌症类型 ===")
        
        # 创建修复后的数据副本
        fixed_df = df.copy()
        
        # 修复COUNTY_FIPS格式，确保是5位数字的字符串
        fixed_df['COUNTY_FIPS'] = fixed_df['COUNTY_FIPS'].apply(
            lambda x: f"{int(x):05d}" if pd.notna(x) and str(x).isdigit() else x
        )
        
        # 修复EQI_Period格式，确保是字符串格式
        fixed_df['EQI_Period'] = fixed_df['EQI_Period'].apply(
            lambda x: '0005' if str(x) == '5' else 
                     '0610' if str(x) == '610' else str(x)
        )
        
        # 定义癌症类型的正确顺序
        cancer_type_order = [
            'C00_C97',   # All Cancers
            'C15_C26',   # Digestive System
            'C18_C21',   # Colorectal
            'C25',       # Pancreatic
            'C30_C39',   # Respiratory System
            'C34',       # Lung and Bronchus
            'C50',       # Female Breast
            'C51_C58',   # Female Genital System
            'C60_C63',   # Male Genital System
            'C61',       # Prostate
            'C64_C68',   # Urinary System
            'C76_C80',   # Other and Unspecified Primary Sites
            'C81_C96'    # Hematopoietic and Lymphoid Tissues
        ]
        
        # 按照指定顺序重新排列数据
        fixed_df['Cancer_Type'] = pd.Categorical(fixed_df['Cancer_Type'], categories=cancer_type_order, ordered=True)
        fixed_df = fixed_df.sort_values(['Cancer_Type', 'Time_Period', 'Lag_Years', 'COUNTY_FIPS']).reset_index(drop=True)
        
        logger.info(f"修复后数据形状: {fixed_df.shape}")
        
        return fixed_df
    
    def restore_smoking_rate(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        从原始数据中恢复正确的Smoking_Rate值
        
        参数:
            df (DataFrame): 处理后的数据
            
        返回:
            DataFrame: 恢复了正确Smoking_Rate的数据
        """
        logger.info("恢复正确的Smoking_Rate值")
        
        # 读取原始数据中的Smoking_Rate值
        original_df = pd.read_csv("/Users/maguoli/Repository/WDP/Data/Processed/df_EQI_AAMR/EQI_AAMR_Point.csv")
        
        # 确保COUNTY_FIPS是字符串格式，便于匹配
        original_df['COUNTY_FIPS'] = original_df['COUNTY_FIPS'].astype(str).str.zfill(5)
        df['COUNTY_FIPS'] = df['COUNTY_FIPS'].astype(str).str.zfill(5)
        
        # 创建一个映射字典，以COUNTY_FIPS为键，Smoking_Rate为值
        smoking_rate_map = dict(zip(original_df['COUNTY_FIPS'], original_df['Smoking_Rate']))
        
        # 使用映射字典更新处理后数据的Smoking_Rate列
        df['Smoking_Rate'] = df['COUNTY_FIPS'].map(smoking_rate_map)
        
        logger.info("Smoking_Rate值恢复完成")
        return df
    
    def save_result(self, df: pd.DataFrame, output_path: str):
        """
        保存结果到CSV文件
        
        参数:
            df (DataFrame): 要保存的数据
            output_path (str): 输出文件路径
        """
        logger.info(f"保存结果到: {output_path}")
        df.to_csv(output_path, index=False)
        logger.info("结果保存完成")

def main():
    """
    主函数
    """
    # 创建处理器实例
    processor = MIResultProcessor()
    
    # 定义文件路径
    input_file = "/Users/maguoli/Repository/WDP/Data/Processed/df_EQI_AAMR/EQI_AAMR_Point_MICE.csv"
    output_file = "/Users/maguoli/Repository/WDP/Data/Processed/df_EQI_AAMR/EQI_AAMR_Point_MICE.csv"
    
    # 加载数据
    df = processor.load_mice_data(input_file)
    
    # 转换为正确格式
    fixed_df = processor.convert_to_long_format(df)
    
    # 恢复正确的Smoking_Rate值
    fixed_df = processor.restore_smoking_rate(fixed_df)
    
    # 保存结果
    processor.save_result(fixed_df, output_file)
    
    logger.info("处理完成")

if __name__ == "__main__":
    main()