import os
import argparse
import pandas as pd

def get_results(out_dir):

    inj_dirs = []
    for file_name in os.listdir(out_dir):
        if file_name.startswith("inj_"):
            inj_dirs.append(file_name)

    df_combined = pd.read_csv(f'{out_dir}/{inj_dirs[0]}/inj_results.csv', index_col=0)
    df_combined['inj_index'] = inj_dirs[0]

    for inj_name in inj_dirs[1:]:
        df_i = pd.read_csv(f'{out_dir}/{inj_name}/inj_results.csv', index_col=0)
        df_i['inj_index'] = inj_name

        df_combined = pd.concat([df_combined, df_i])

    df_combined = df_combined.reset_index(drop=True)
    df_combined.to_csv(f'{out_dir}/combined_inj_results.csv')

if __name__=='__main__':
    parser = argparse.ArgumentParser(prog='Result gatherer for offline injection pipeline',
                                     epilog='Feel free to contact me if you have questions - rsenzel@mpifr-bonn.mpg.de')

    parser.add_argument('--out_dir', metavar='dir', required=True, help='output directory')
    args = parser.parse_args()

    get_results(args.out_dir)