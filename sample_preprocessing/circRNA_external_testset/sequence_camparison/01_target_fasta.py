import pandas as pd

import os

def csv_to_fasta(csv_path, output_fasta, id_col='circID', seq_col='Sequence'):

    df = pd.read_csv(csv_path)

    with open(output_fasta, 'w') as f:

        for idx, row in df.iterrows():

                               

            seq_id = str(row[id_col]).replace(" ", "_").replace("/", "_")

            seq = row[seq_col]

            f.write(f">{seq_id}\n{seq}\n")

    print(f"Created {output_fasta} with {len(df)} sequences.")

if __name__ == "__main__":

    out_dir = "./outputs"

    os.makedirs(out_dir, exist_ok=True)

    

                   

    exo_csv = "circExor/reference_preprocessing/heldout_circRNA_3.20/exoRbase_output.bak.csv"

    csv_to_fasta(exo_csv, f"{out_dir}/target_exoRbase.fasta", id_col='circID', seq_col='Sequence')

    

                  

    cell_csv = "circExor/reference_preprocessing/heldout_circRNA_3.20/filtered_celline_arrary_df.csv"

    csv_to_fasta(cell_csv, f"{out_dir}/target_celline.fasta", id_col='CircRNA', seq_col='Sequence')
