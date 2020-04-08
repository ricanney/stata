# bim2dta
__Author:__ Richard Anney
__Edited:__ 27 Mar 2020
__Current Version:__ v3

__Overview__
This function converts the plink binary ``.bim`` into a Stata format ``.dta`` file. Variables are named, formatted and additional naming and genotype calls are added.

__Dependencies__
The function has the following dependencies are ;
``bim_path_to_short``  ``bim2count`` ``checkfile`` ``checkloc_name`` ``files2dta`` ``recodegenotype``

__Input Files__
For more details on input files see [input files](https://github.com/ricanney/stata/blob/master/input-files.md).
This function requires the ``*.bim`` file

__Syntax__
```
bim2dta, bim(C:/data/plink_binary)
```
or
```
local bim C:/data/plink_binary
bim2dta, bim(`bim')
```

You define the ``.bim`` file in the script using ``bim()`` - however, you ommit the ``*.bim`` from the file name. The above example is for the file path ``C:/data/plink_binary.bim``


__Output Files__
The resulting file is created with the ``_bim.dta`` suffix 
i.e. ``plink_binary.bim`` creates ``plink_binary_bim.dta``.

The file is deposited in the root directory in a sub-directory ``/bim2dta/<filename>``. 
For ``plink_binary.bim`` the subdirectory will be ``<root>/bim2dta/plink_binary``.

__bim.dta file format__

The ``_bim.dta`` file contains the following  variables.

| Variable Name | Storage Type | Variable Label | Notes |
|--------|--------|--------|--------|
|``chr``|str2|Chromosome|Chromosomes are numbered 1-25|
|``snp``|str10|Marker Name| |
|``bp``|str9|Physical Location| |
|``a1``|str1|Allele 1 (``ACGT(D/I)``)| convert long indel naming to ``I`` or ``D``
|``a2``|str1|Allele 2 (``ACGT(D/I)``)| convert long indel naming to ``I`` or ``D``
|``gt``|str2|Genotype (``RYMKSW(ID))``| from ``recodegenotypes``
|``loc_name``|str18|Marker loc_name (chr<#>:<pos>-<gt>)|from ``checkloc_name``

__example log file__
```
#########################################################################
# > .................. bim2dta
# > ........... author Richard Anney
# > ... report bugs to anneyr@cardiff.ac.uk
# >  code available at https://github.com/ricanney
# > .......... license CC BY-NC (CC Attribution-NonCommercial)
# > .................. all functions considered beta "caveat emptor"
#########################################################################
# Started:             27 Mar 2020 16:05:22 
# Username:            Richard Anney
# Operating System:    Windows
#########################################################################
# > ..... project name hapmap3_ceu_tsi
# > ..... SNPs in file 1464434
# > ...... saved graph hapmap3_ceu_tsi_bim.dta
#########################################################################
# > overview of bim2dta
#########################################################################
# > 1. import bim file 
# > 2. create gt variable
# > 3. convert indels to I/D
# > 4. convert chromosomes to numeric codes
# > 5. create loc_name variable
#########################################################################
```