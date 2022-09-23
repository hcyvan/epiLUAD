import pandas as pd
import numpy as np
import time
import os
import argparse

parser = argparse.ArgumentParser(description='meth matrix filter and select')
parser.add_argument('-i', '--input', required=True, help='The input bed ratio file with header')
parser.add_argument('-t', '--tss', required=True, help='The tss bed file: chrom, start, end')
parser.add_argument('-u', '--up', default=10000, type=int, help='tss upstream')
parser.add_argument('-d', '--down', default=10000, type=int, help='tss downstream')
parser.add_argument('-o', '--out', help='The out matrix')
parser.add_argument('-v', '--version', action='version', version='0.1')

args = parser.parse_args()
output = args.out
tss_file_txt = args.tss
ratio_txt = args.input
UP = args.up
DOWN = args.down

tss_file = pd.read_csv(tss_file_txt, sep='\t')
ratio = pd.read_csv(ratio_txt, sep='\t')

sample = ratio.columns[3:].to_list()
SAMPLE = len(sample)

top_count = np.zeros((UP + DOWN + 1, SAMPLE))
top_ratio = np.zeros((UP + DOWN + 1, SAMPLE))

ratio_map = dict()
for chrom, data in ratio.groupby("#chrom"):
    ratio_map[chrom] = data

total = len(tss_file)
start = time.time()
idx = 0
for i in range(len(tss_file)):
    if idx % 10000 == 0:
        print("index {}/{}, pass {} min".format(idx, total, round((time.time() - start) / 60), 4))
    idx += 1
    chrom = tss_file.iloc[i, 0]
    tss = tss_file.iloc[i, 1] - 1
    sign = tss_file.iloc[i, 2]
    ratio2 = ratio_map[chrom]
    interval = [tss - UP, tss + DOWN]
    if sign == '-':
        interval = [tss - DOWN, tss + UP]
    data = ratio2.loc[(ratio2['start'] >= interval[0]) & (ratio2['start'] <= interval[1])]
    left = pd.DataFrame({
        'chrom': chrom,
        'start': list(range(interval[0], interval[1] + 1))
    })
    a = pd.merge(left, data, how='left')
    if sign == '-':
        a = a.iloc[::-1]
    a[a == -1] = np.nan
    mm = a.iloc[:, 4:].to_numpy()
    mm_count = mm.copy()
    mm_count[~np.isnan(mm)] = 1
    mm_count[np.isnan(mm)] = 0
    top_count += mm_count
    mm[np.isnan(mm)] = 0
    top_ratio += mm
np.seterr(invalid='ignore')
ratio_matrix = top_ratio / top_count
matrix = pd.DataFrame(np.round_(ratio_matrix, 2), columns=ratio.columns[3:])
matrix['tss'] = list(range(-UP, DOWN + 1))
matrix.set_index('tss', inplace=True)

with open(output, "w") as f:
    matrix.to_csv(f, sep="\t")

end = time.time()
print()
print("use {} min !!".format((end - start) / 60))
