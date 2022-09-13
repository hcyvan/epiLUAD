import argparse
import sys

parser = argparse.ArgumentParser(description='merge dmc to dmr')
parser.add_argument('-i', '--input', required=True, help='dmc bed file. "chrom, start, end, class" is necessary')
parser.add_argument('-m', '--minimum', default=3, help='minimum dmc in the dmr')
parser.add_argument('-d', '--dist', default=200, type=int, help='maximum distance between two dmc')
parser.add_argument('-o', '--output', help='dmr bed file. header: "chrom, start, end, class, dmc_number, length"')
parser.add_argument('-v', '--version', action='version', version='0.1')
args = parser.parse_args()

dmc_file = args.input
dist = args.dist
minimum = args.minimum
dmr_file = args.output

dmr_chrom = None
dmr_start = None
dmr_end = None
dmr_class = None
dmr_cpg_num = 0

if dmr_file:
    out = open(dmr_file, 'w')
else:
    out = sys.stdout

for line in open(dmc_file, 'r'):
    dmc = line.strip().split('\t')
    if dmr_cpg_num == 0:
        dmr_chrom = dmc[0]
        dmr_start = dmc[1]
        dmr_end = dmc[2]
        dmr_class = dmc[3]
        dmr_cpg_num = 1
    else:
        if dmr_chrom == dmc[0] and int(dmc[2]) - int(dmr_end) <= dist and dmr_class == dmc[3]:
            dmr_end = dmc[2]
            dmr_cpg_num += 1
        else:
            if dmr_cpg_num >= minimum:
                out.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(dmr_chrom, dmr_start, dmr_end, dmr_class, dmr_cpg_num, int(dmr_end) - int(dmr_start)))
            dmr_chrom = dmc[0]
            dmr_start = dmc[1]
            dmr_end = dmc[2]
            dmr_class = dmc[3]
            dmr_cpg_num = 1
