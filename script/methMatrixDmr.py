import gzip
import time
import numpy as np
import argparse
import sys
# TODO: remove class column
parser = argparse.ArgumentParser(description='Extract methylation level of DMR')
parser.add_argument('-i', '--input', required=True, help='The input bgzip bed')
parser.add_argument('-b', '--bed', required=True, help='The input dmr bed')
parser.add_argument('-o', '--out', help='The out file')
parser.add_argument('-v', '--version', action='version', version='0.1')

args = parser.parse_args()

gzip_file = args.input
bed_file = args.bed
output = args.out

if output:
    out = open(output, 'w')
else:
    out = sys.stdout

gzip_f = gzip.open(gzip_file, 'r')
headers = gzip_f.readline().decode().strip().split('\t')
ncol = len(headers) - 3
top = np.empty((0, ncol), np.float64)

bed_f = open(bed_file, 'r')
bed = bed_f.readline()
while bed.startswith('#'):
    bed = bed_f.readline()
bed = bed.strip().split('\t')

headers_new = ['chrom', 'start', 'end', 'class']
headers_new.extend(headers[3:])
out.write('\t'.join(headers_new) + '\n')
out.flush()
start = time.time()
idx = 0
for line in gzip_f:
    if idx % 100000 == 0:
        sys.stderr.write("index {}, pass {} min\n".format(idx, round((time.time() - start) / 60), 4))
    idx += 1
    line = line.decode()
    if line.startswith('#'):
        continue
    line = line.strip().split('\t')
    if line[0] != bed[0] or int(line[1]) < int(bed[1]) or int(line[1]) > int(bed[2]) - 1:
        if top.shape[0] > 0:
            mean_ratio = np.round(np.nanmean(top, axis=0), 3).astype('str')
            dmr_line = bed[0:4]
            dmr_line.extend(mean_ratio)
            out.write('\t'.join(dmr_line) + '\n')
            out.flush()
            top = np.empty((0, ncol), np.float64)
            bed = bed_f.readline().strip().split('\t')
        else:
            if line[0] == bed[0] and int(line[1]) > int(bed[2]) - 1:
                bed = bed_f.readline().strip().split('\t')
    else:
        ratios = np.array([line[3:]]).astype(np.float64)
        ratios[ratios == -1] = np.nan
        top = np.append(top, ratios, axis=0)

end = time.time()
sys.stderr.write('\n')
sys.stderr.write("use {} min !!\n".format((end - start) / 60))
