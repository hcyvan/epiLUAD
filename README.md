# epiLUAD
## How to run?
+ config.default.yaml to config.yaml, and set variables
## script

### dmc2dmr.py
Merge SR-DMc to SR-DMR
```
python ./script/dmc2dmr.py -i ./data/intermediate/srdmc.s2.bed -o ./data/intermediate/srdmr.s2.bed
```
### methMatrixExtract.py
Get the methylation level of SR-DMC and SR-DMR from methylation level matrix

```
python ./script/methMatrixExtract.py -b /path/to/merge_ratio_bed.gz -r ./data/intermediate/wgbs/srdmc.bed -o ./data/intermediate/wgbs/srdmc.ratio.bed
```
### methMatrixGenerate.py

Generate methylation level matrix from *.G.bed file of MOABS.

eg:
```
# generate methylation level matrix. Only keep CpG sites with depth greater than or equal to 10 in all samples  
python methMatrixGenerate.py -i G.bed.list -c hg38.all.cpgSite.bed -o ./merge.ratio.d10.bed -r B -d 10 -e one

```