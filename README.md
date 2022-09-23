# epiLUAD

## script

**dmc2dmr.py**

**methMatrixGenerate.py**

Generate methylation level matrix from *.G.bed file of MOABS.

eg:
```shell
# generate methylation level matrix. Only keep CpG sites with depth greater than or equal to 10 in all samples  
python methMatrixGenerate.py -i G.bed.list -c hg38.all.cpgSite.bed -o ./merge.ratio.d10.bed -r B -d 10 -e one

```