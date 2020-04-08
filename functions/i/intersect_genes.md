# intersect_genes
__Author:__ Richard Anney
__Edited:__ 2 Apr 2020
__Current Version:__ v1

__Overview__
This function was created to mimic the ``bedtools`` script ``-intersect``
Using either a single ``range`` or a list of ``range`` co-ordinates, this function identifies intersections of the co-ordinates with known genes. ``range`` is formatted as``chr<#>:<from#>..<to#>`` e.g. chr20:1000000..3000000

The gene co-ordinate file ``Homo_sapiens.GRCh37._version_.gtf_exon.dta`` is used as the ``generef`` - this file can be created using ``graphgene_create_ref``

__Dependencies__
The function has the following dependencies are ;
``checkfile`` ``ref_path_to_short``

__Input Files__
The gene co-ordinate file ``Homo_sapiens.GRCh37._version_.gtf_exon.dta`` - this file can be created using ``graphgene_create_ref``

__Syntax__
```
# single range 
intersect_genes, range(chr20:1000000..3000000) generef(${generef}) 

# multiple ranges
qui { // create range file intersect_genes_range.dta
 clear
 set obs 2
 gen intersect_genes_range = ""
 replace intersect_genes_range = "chr20:1000000..3000000" in 1
 replace intersect_genes_range = "chr20:2000000..4000000" in 2
 save ${init_root}/intersect_genes_range.dta, replace
 }
intersect_genes, range_list(${init_root}/intersect_genes_range.dta) generef(${generef}) 
```

__Output Files__
The files are deposited in the root directory in a sub-directory ``/intersect_genes/Homo_sapiens.GRCh37.87.gtf_exon``. 

1. For the single range script the file is named according to the co-ordinates
```
intersect_genes, range(chr20:1000000..3000000) generef(${generef})
```
The above script would produce the file ``Homo_sapiens.GRCh37.87.gtf_exon-with-chr20_1000000_3000000_intersect_genes.dta``
2. For the multiple ranges - the file is named according to the range-file naming convention
```
intersect_genes, range_list(${init_root}/intersect_genes_range.dta) generef(${generef}) 
```
The above script will produce the file ``Homo_sapiens.GRCh37.87.gtf_exon_intersect_genes_range_intersect_genes.dta``

__*intersect_genes.dta format__

| Variable Name | Storage Type | Variable Label | Notes |
|--------|--------|--------|--------|
|ensembl_geneID|str15|gene_id (ensembl): The stable identifier for the gene|.
|chr|byte|chromosome|.
|start|long|start position of the gene, with sequence numbering starting at 1|.
|end|long|end position of the gene, with sequence numbering starting at 1|.
|symbol|str13|official symbol of this gene|.
|biotype|str14|biotype of this gene|.
|ensembl_geneID2|str17|gene_id (ensembl): The stable identifier for the gene + version|.
|ensembl_txID|str15|The stable identifier for this transcript|.
|exonCount|str2|exon count|.
|txstart|long|start position of the exon, with sequence numbering starting at 1|.
|txend|long|end position of the exon, with sequence numbering starting at 1|.
|label|str15|symbol + strand label (+ = >; - = <)|.
|intersect|str22|overlapping range| formatted as ``range``


