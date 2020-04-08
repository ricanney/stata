# create_generef
__Author:__ Richard Anney
__Edited:__ 31 Mar 2020
__Current Version:__ v1

__Overview__
This function creates the ``Homo_sapiens.GRCh37.`version'.gtf_gene.dta`` and ``Homo_sapiens.GRCh37.`version'.gtf_exon.dta`` files for use in functions ``graphgene``. The package runs on the linux and windows 10 (with ubuntu bash activates) operating systems. ``create_generef`` uses the ``wget`` and ``gunzip`` commands to download and unpack the raw data.

Optional: module 5 can be unlocked -``remove /* /*`` to generate ``*.hg-list`` plain text files are created for use in ``bedtools`` and ``plink`` - these are split by ``biotype``

__Dependencies__
The function has the following dependencies are ;
``checkfile``

__Input Files__
None

__Syntax__
```
create_generef
```
or
```
create_generef, version(87)
* this allows specific builds to be downloaded from ftp://ftp.ensembl.org/pub/grch37/release-`version'/gtf/homo_sapiens/Homo_sapiens.GRCh37.`version'.gtf.gz
```

__Output Files__
1. ``Homo_sapiens.GRCh37.`version'.gtf_gene.dta``
2. ``Homo_sapiens.GRCh37.`version'.gtf_exon.dta``

The file is deposited in the root directory in a sub-directory ``/create_generef/Homo_sapiens.GRCh37.`version'.gtf``. 

__Homo_sapiens.GRCh37.`version'.gtf_gene.dta file format__

| Variable Name | Storage Type | Variable Label | Notes |
|--------|--------|--------|--------|
|ensembl_geneID|str2|gee_id (ensembl): The stable identifier for the gene|.|
|chr|byte|chromosome|.|
|start|long|start position of the feature, with sequence numbering starting at 1 |.|
|end|long|end position of the feature, with sequence numbering starting at 1|.|
|strand|str1|defined as + (forward) or - (reverse)|.|
|symbol|str22|official symbol of this gene|.|
|biotype|str24|biotype of this gene|.|
|label|str22|symbol + strand label (+ = >; - = <)|.|

__Homo_sapiens.GRCh37.`version'.gtf_exon.dta file format__

| Variable Name | Storage Type | Variable Label | Notes |
|--------|--------|--------|--------|
|ensembl_geneID|str15|gene_id (ensembl): The stable identifier for the gene|.|
|chr|byte|chromosome|.|
|start|long|start position of the gene, with sequence numbering starting at 1|.|
|end|long|end position of the gene, with sequence numbering starting at 1|.|
|symbol|str22|official symbol of this gene|.|
|biotype|str24|biotype of this gene|.|
|ensembl_geneID2|str26|gene_id (ensembl): The stable identifier for the gene + version|.|
|ensembl_txID|str15|The stable identifier for this transcript|.|
|exonCount|str3|exon count|.|
|txstart|long|start position of the exon, with sequence numbering starting at 1|.|
|txend|long|end position of the exon, with sequence numbering starting at 1|.|
|label|str24|symbol + strand label (+ = >; - = <)|.|


