# bim2rymk
__Author:__ Richard Anney
__Edited:__ 31 Mar 2020
__Current Version:__ v1

__Overview__
This function runs the ``bim2dta`` package and limits the SNPs to non-ambiguous set, specifically those with the genotypes ``R = G|A`` ``Y=C|T`` ``M=C|A`` ``K=G|T``.

__Dependencies__
The function has the following dependencies are ;
``bim_path_to_short`` ``bim2count`` ``bim2dta`` ``checkfile`` ``files2dta`` 

__Input Files__
For more details on input files see [input files](https://github.com/ricanney/stata/blob/master/input-files.md).
This function requires the ``*.bim`` file

__Syntax__

```
bim2rymk, bim(C:/data/plink_binary)
```
or
```
local bim C:/data/plink_binary
bim2rymk, bim(`bim') pc(5)
```

__Output Files__

The resulting script creates two files;
1. ``*_rymk.dta`` 

The file is deposited in the root directory in a sub-directory ``/bim2rymk/<filename>``. 
For ``plink_binary.bim`` the subdirectory will be ``<root>/bim2rymk/plink_binary``.

__*_rymk.dta file format__

| Variable Name | Storage Type | Variable Label | Notes |
|--------|--------|--------|--------|
|``snp``|str10|Marker Name| |
|``a1``|str1|Allele 1 (``ACGT(D/I)``)| convert long indel naming to ``I`` or ``D``
|``a2``|str1|Allele 2 (``ACGT(D/I)``)| convert long indel naming to ``I`` or ``D``
|``gt``|str2|Genotype (``RYMKSW(ID))``| from ``recodegenotypes``


__example log file__
```
#########################################################################
# > .................. bim2rymk
# > ........... author Richard Anney
# > ... report bugs to anneyr@cardiff.ac.uk
# >  code available at https://github.com/ricanney
# > .......... license CC BY-NC (CC Attribution-NonCommercial)
# > .................. all functions considered beta "caveat emptor"
#########################################################################
# Started:             31 Mar 2020 11:53:46 
# Username:            Richard Anney
# Operating System:    Windows
#########################################################################
# > ..... project name neonate-ahe-qc-v9
# > .... plink version H:/software/plink/bin/v1.9/plink_1.9b6.16-19feb2020.exe
# > . SNPs in original 679251
# >  SNPs in processed 677236
#########################################################################
# > overview of bim2rymk
#########################################################################
# > 1. use / create _bim.dta 
# > 2. remove on R Y M K genotypes
#########################################################################
```