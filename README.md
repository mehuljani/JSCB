# JSCB

Usage:
```
python jscb.py inputfile.gb
```

Input:
JS-CB takes genbank(full) as input file

Output:
JS-CB outputs three tab separated files:
1) JSCB_output.tsv: Has three columns- Genomic Island (GI)-ID, GI-start co-ordinate and GI-end co-ordinate.
2) JSCB_ouput.gi: Contains four columns- Gene number, Cluster id, Cluster size and Gene length.
3) JSCB_output.clus: Contains four columns- Cluster id, Cluster size, Gene number, and Gene length.

Genomic island prediction can be found in JSCB_output.tsv file. JSCB_output.gi/.clus files contains cluster configuration information.

Requirements:
The program requires Biopython. The binary JSCB file is compiled using F95 (Fortran95) compiler
