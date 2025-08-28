import pandas as pd

def extract_pressure_columns(input_file: str, output_file: str):
    """
    从气压数据集中提取指定字段，并保存为新的CSV文件。
    :param input_file: 原始CSV文件路径
    :param output_file: 输出CSV文件路径
    """
    # 需要保留的字段
    columns_to_keep = [
        "staPresMean",
        "staPresMinimum",
        "staPresMaximum",
        "corPres",
        "corPresExpUncert"
    ]

    # 读取CSV
    df = pd.read_csv(input_file)

    # 检查字段是否存在
    missing = [col for col in columns_to_keep if col not in df.columns]
    if missing:
        raise ValueError(f"输入文件缺少字段: {missing}")

    # 只保留需要的列
    df_filtered = df[columns_to_keep]

    # 保存新文件
    df_filtered.to_csv(output_file, index=False)
    print(f"提取完成，结果已保存到 {output_file}")


if __name__ == "__main__":
    # 修改成你的文件路径
    input_file = "//dataset/NEON_barometric_pressure.csv"  # 原始数据
    output_file = "//dataset/barometric_pressure.csv"  # 提取结果
    extract_pressure_columns(input_file, output_file)
