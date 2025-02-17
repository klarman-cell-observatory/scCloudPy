import time

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from sccloud.io import read_input


def convert_to_parquet(data, output_name, nthreads):
    data.obs.index.name = "sccloud.cell.barcode"
    df = data.obs.reset_index()
    whitelist = ['X_pca', 'X_rpca', 'X_tsne', 'X_fitsne', 'X_umap', 'X_fle', 'X_net_tsne', 'X_net_umap', 'X_net_fle']
    whitelist_3d = ['X_diffmap_pca']
    for key in data.obsm.keys():
        if key in whitelist:
            df["{}_1".format(key)] = data.obsm[key][:, 0]
            df["{}_2".format(key)] = data.obsm[key][:, 1]
        elif key in whitelist_3d:
            df["{}_1".format(key)] = data.obsm[key][:, 0]
            df["{}_2".format(key)] = data.obsm[key][:, 1]
            df["{}_3".format(key)] = data.obsm[key][:, 2]

    metadata_table = pa.Table.from_pandas(df, nthreads=nthreads)

    df_expr = pd.DataFrame(data=data.X.toarray(), columns=data.var_names)
    parquet_table = pa.Table.from_pandas(df_expr, nthreads=nthreads)

    for i in range(metadata_table.num_columns - 1, 0, -1):
        column = metadata_table[i]
        if column.name != "__index_level_0__":
            parquet_table = parquet_table.add_column(0, column)

    output_file = output_name + ".parquet"
    pq.write_table(parquet_table, output_file)
    print(output_file + " is written!")


def run_conversion(input_h5ad_file, output_name, nthreads):
    start = time.time()
    data = read_input(input_h5ad_file)
    end = time.time()
    print(
        "Time spent for loading the expression matrix is {:.2f}s.".format(end - start)
    )

    start = time.time()
    convert_to_parquet(data, output_name, nthreads)
    end = time.time()
    print("Time spent on generating the PARQUET file is {:.2f}s.".format(end - start))
