# bim2frq
__Author:__ Richard Anney
__Edited:__ 30 Mar 2020
__Current Version:__ v3

__Overview__
This function applies the ``plink`` function ``--freq`` to calculate the allele frequency within the defined plink binaries. The resulting ``*.frq`` file is converted into  into a Stata format ``.dta`` file. Variables are named, formatted and additional naming and genotype calls are added.

__Dependencies__
The function has the following dependencies are ;
``bim_path_to_short`` ``checkfile`` ``files2dta`` ``recodegenotype``
``checktabbed`` - this function requires conversion of ``plink`` pretty format to useable tab-delimited.

__Input Files__
For more details on input files see [input files](https://github.com/ricanney/stata/blob/master/input-files.md).
This function requires the ``*.bed`` ``*.bim`` and ``*.fam`` file

__Syntax__
```
bim2frq, bim(C:/data/plink_binary)
```
or
```
local bim C:/data/plink_binary
bim2frq, bim(`bim')
```

You define the plink binaries (with common name) in the script using ``bim()`` - however, you ommit the ``*.bim`` from the file name. The above example is for the file path ``C:/data/plink_binary.bim``

__Output Files__
The resulting file is created with the ``_frq.dta`` suffix 
i.e. ``plink_binary.bim`` creates ``plink_binary_frq.dta``.

The file is deposited in the root directory in a sub-directory ``/bim2frq/<filename>``. 
For ``plink_binary.bim`` the subdirectory will be ``<root>/bim2frq/plink_binary``.

__frq.dta file format__

The ``_frq.dta`` file contains the following  variables.

| Variable Name | Storage Type | Variable Label | Notes |
|--------|--------|--------|--------|
|``snp``|str10|Marker Name| |
|``a1``|str1|Allele 1 (``ACGT(D/I) - minor``)| convert long indel naming to ``I`` or ``D``
|``a2``|str1|Allele 2 (``ACGT(D/I)``)| convert long indel naming to ``I`` or ``D``
|``maf``|double|Minor Allele Frequency | based on ``a1``
|``gt``|str2|Genotype (``RYMKSW(ID))``| from ``recodegenotypes``
|``var``|double|Variance (``((maf*(1-maf))/nchr)``)| derived variance for next-generation based on maf and number of chromosomes

__example log file__
```
#########################################################################
# > .................. bim2frq
# > ........... author Richard Anney
# > ... report bugs to anneyr@cardiff.ac.uk
# >  code available at https://github.com/ricanney
# > .......... license CC BY-NC (CC Attribution-NonCommercial)
# > .................. all functions considered beta "caveat emptor"
#########################################################################
# Started:             30 Mar 2020 11:07:29 
# Username:            Richard Anney
# Operating System:    Windows
#########################################################################
# > ..... project name hapmap3_ceu_tsi
# > .... plink version H:/software/plink/bin/v1.9/plink_1.9b6.16-19feb2020.exe
# > ...... SNPs in bim 1464434
# > ...... SNPs in frq 1464434
# > ....... saved data hapmap3_ceu_tsi_frq.dta
#########################################################################
# > overview of bim2frq
#########################################################################
# > 1. identify chromosome range from bim
# > 2. calculate frq by chromosome
# > 3. collate / append files
#########################################################################

```