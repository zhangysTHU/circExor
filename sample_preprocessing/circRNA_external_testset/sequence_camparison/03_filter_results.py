import pandas as pd
import os

def filter_csv(original_csv, matched_ids_file, output_csv, id_col='circID'):
    if not os.path.exists(matched_ids_file):
        print(f"Error: {matched_ids_file} not found.")
        return

    with open(matched_ids_file, 'r') as f:
        matched_ids = set(line.strip() for line in f if line.strip())

    df = pd.read_csv(original_csv)
    original_count = len(df)
    
    # 格式化需一致以便比对
    df['normalized_id'] = df[id_col].astype(str).str.replace(" ", "_").str.replace("/", "_")
    
    # 过滤掉存在于 matched_ids 中的行
    filtered_df = df[~df['normalized_id'].isin(matched_ids)].drop(columns=['normalized_id'])
    
    filtered_df.to_csv(output_csv, index=False)
    
    print(f"[Results for {os.path.basename(original_csv)}]")
    print(f"  - 原始记录: {original_count}")
    print(f"  - 剔除记录: {original_count - len(filtered_df)}")
    print(f"  - 剩余记录: {len(filtered_df)}")
    print(f"  - 文件已保存至: {output_csv}\n")

if __name__ == "__main__":
    out_dir = "./outputs"
    
    configs = [
        {
            "csv": "/lulabdata3/huangkeyun/zhangys/RNA_locator/ML_python_scripts/reference_preprocessing/heldout_circRNA_3.20/exoRbase_output.bak.csv",
            "lst": f"{out_dir}/target_exoRbase_matched_ids.lst",
            "out": f"{out_dir}/exoRbase_filtered.csv",
            "id_col": "circID"
        },
        {
            "csv": "/lulabdata3/huangkeyun/zhangys/RNA_locator/ML_python_scripts/reference_preprocessing/heldout_circRNA_3.20/filtered_celline_arrary_df.csv",
            "lst": f"{out_dir}/target_celline_matched_ids.lst",
            "out": f"{out_dir}/celline_filtered.csv",
            "id_col": "CircRNA"
        }
    ]
    
    for cfg in configs:
        filter_csv(cfg["csv"], cfg["lst"], cfg["out"], id_col=cfg["id_col"])