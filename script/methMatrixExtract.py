#!/usr/bin/env python
import sys
import numpy as np
import pandas as pd
from collections import OrderedDict
import gzip
import time
import argparse


def parse_arguments():
    parser = argparse.ArgumentParser(prog='methMatrixDmr2.py', description='extract dmr from matrix')
    parser.add_argument('-v', '--version', action='version', version='0.1')
    parser.add_argument('-q', '--quiet', default=False, type=bool, help='print run details to stderr')
    parser.add_argument('-b', '--g-bed', required=True, type=argparse.FileType('r'),
                        help='The G.bed file generate by mcall or G.bed matrix generate by methMatrixGenerate.py')
    parser.add_argument('-r', '--region-bed', required=True, type=argparse.FileType('r'),
                        help='The region to extract region methy ratio')
    parser.add_argument('-o', '--out-file', required=True, type=argparse.FileType('w'), help='The output file')
    return parser.parse_args()


class Reader:
    def __init__(self, file_name, is_zip=None):
        if is_zip is None:
            if file_name.endswith('.gz'):
                self.__is_zip = True
            else:
                self.__is_zip = False
        else:
            self.__is_zip = is_zip
        if self.__is_zip:
            self._f = gzip.open(file_name, 'r')
        else:
            self._f = open(file_name, 'r')

    def readline(self):
        line = self._f.readline()
        if line:
            if self.__is_zip:
                return line.decode()
            else:
                return line
        return None

    def __iter__(self):
        return self

    def __next__(self):
        line = self._f.readline()
        if line:
            if self.__is_zip:
                return line.decode()
            else:
                return line
        else:
            raise StopIteration

    def close(self):
        self._f.close()


class Meth:
    def __init__(self, meth_str):
        if isinstance(meth_str, str):
            items = meth_str.split('\t')
        else:
            items = meth_str
        if len(items) < 4:
            raise Exception('Meth parse error: {}'.format(meth_str))
        self.chrom = items[0]
        self.start = int(items[1])
        self.end = int(items[2])
        self.levels = [float(x) for x in items[3:]]

    def __repr__(self):
        return "[{} {} {} {}]".format(self.chrom, self.start, self.end, ' '.join([str(i) for i in self.levels]))


class Region:
    def __init__(self, items, samples=1):
        if isinstance(items, str):
            items = items.split('\t')
        self.chrom = items[0]
        self.start = int(items[1])
        self.end = int(items[2])
        self.ratio_array = []
        self.samples = samples

    def __repr__(self):
        return "<Region {}:{}-{}>".format(self.chrom, self.start, self.end)

    def add_meth(self, meth: Meth, min_cpg=1):
        """
        :param meth:
        :param min_cpg:
        :return:
            -1: meth is on the left of this region and not counted
            1: meth is on the right of this region and not counted
            0: meth is counted in this region
        """
        if meth.start < self.start:
            return -1
        elif meth.start >= self.end:
            return 1
        else:
            self.ratio_array.append(meth.levels)
            return 0

    def parse_ratios(self):
        if len(self.ratio_array):
            arr = np.array(self.ratio_array)
            arr[arr == -1] = np.nan
            return np.nanmean(arr, axis=0)
        else:
            return [np.nan] * self.samples


class RegionArray:
    def __init__(self):
        self.arr = []
        self.top = 0

    def append(self, region: Region):
        self.arr.append(region)

    def add_meth(self, meth: Meth, min_cpg=1):
        for i in range(self.top, len(self.arr)):
            region: Region = self.arr[i]
            ret = region.add_meth(meth, min_cpg)
            if ret < 0:
                return
            elif ret > 0:
                self.top = i


class RegionMap:
    def __init__(self, region_bed, samples=1, header=None):
        self.map: OrderedDict[str, RegionArray] = OrderedDict()
        r_bed = pd.read_csv(region_bed, sep='\t', header=None)
        for chr in r_bed[0].unique().tolist():
            if chr.startswith('#'):
                continue
            bed = r_bed[r_bed[0] == chr]
            bed = bed.astype({1: int, 2: int})
            bed = bed.sort_values(by=[1])
            regions = RegionArray()
            for row in bed.values:
                regions.append(Region(row, samples=samples))
            self.map[chr] = regions
            self.header=header

    def print_info(self):
        total = 0
        for k, v in self.map.items():
            region_num = len(v.arr)
            total += region_num
            print('{}\t{} regions'.format(k, region_num))
        print()
        print("total:\t {} regions".format(total))

    def add_meth(self, meth: Meth, min_cpg=1):
        regions = self.map.get(meth.chrom)
        if regions:
            regions.add_meth(meth, min_cpg)

    #
    def print(self, of=sys.stdout, order=None):
        if self.header:
            of.write(self.header + "\n")
        if order:
            for k in order:
                v = self.map.get(k)
                if v:
                    for region in v.arr:
                        level = "\t".join([str(x) for x in region.parse_ratios()])
                        of.write('{}\t{}\t{}\t{}\n'.format(region.chrom, region.start, region.end, level))
        else:
            for k, v in self.map.items():
                for region in v.arr:
                    level = "\t".join([str(x) for x in region.parse_ratios()])
                    of.write('{}\t{}\t{}\t{}\n'.format(region.chrom, region.start, region.end, level))


def extract_region_methy_level(g_bed_file, region_bed, out_file, verbose=False):
    start = time.time()
    g_bed_f = Reader(g_bed_file)
    line = g_bed_f.readline()
    header = None
    if line.startswith("#chrom"):
        header = line.strip()
        line = g_bed_f.readline()
    g_bed = line.strip().split('\t')
    rmap = RegionMap(region_bed, samples=len(g_bed) - 3, header=header)
    rmap.print_info()
    i = 0
    while True:
        if verbose and i % 5000000 == 0:
            sys.stderr.write("iter:{}\tuse:{}s\n".format(i, (int(time.time() - start))))
        i += 1
        meth = Meth(g_bed)
        rmap.add_meth(meth)
        line = g_bed_f.readline()
        if not line:
            break
        g_bed = line.strip().split('\t')
    sys.stderr.write("Finished! use:{}s\twriting to {}\n".format(int(time.time() - start), out_file))
    with open(out_file, 'w') as f:
        rmap.print(of=f)


if __name__ == '__main__':
    args = parse_arguments()
    extract_region_methy_level(g_bed_file=args.g_bed.name,
                               region_bed=args.region_bed.name,
                               out_file=args.out_file.name,
                               verbose=(not args.quiet))
