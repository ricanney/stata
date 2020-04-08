# fam2dta
__Author:__ Richard Anney
__Edited:__ 30 Mar 2020
__Current Version:__ v2

__Overview__
This function converts the plink binary ``.fam`` into a Stata format ``.dta`` file. Variables are named and formatted.

__Dependencies__
The function has the following dependencies are ;
``bim_path_to_short``  ``bim2count`` ``checkfile`` ``files2dta`` 

__Input Files__
For more details on input files see [input files](https://github.com/ricanney/stata/blob/master/input-files.md).
This function requires the ``*.fam`` file

__Syntax__
```
fam2dta, fam(C:/data/plink_binary)
```
or
```
local fam C:/data/plink_binary
fam2dta, fam(`fam')
```

You define the ``.fam`` file in the script using ``fam()`` - however, you ommit the ``*.fam`` from the file name. The above example is for the file path ``C:/data/plink_binary.fam``


__Output Files__
The resulting file is created with the ``_fam.dta`` suffix 
i.e. ``plink_binary.fam`` creates ``plink_binary_fam.dta``.

The file is deposited in the root directory in a sub-directory ``/fam2dta/<filename>``. 
For ``plink_binary.fam`` the subdirectory will be ``<root>/fam2dta/plink_binary``.

__fam.dta file format__

The ``_fam.dta`` file contains the following  variables.

| Variable Name | Storage Type | Variable Label | Notes |
|--------|--------|--------|--------|
|``fid``|str#|Family Identifier|.
|``iid``|str#|Individual Identifier|.|
|``fatid``|str#|Father's (within family) Identifier|.|
|``motid``|str#|Mother's (within family) Identifier|.|
|``sex``|byte|Sex|Value label added ``1=Male, 2=Female, -9=Missing``|
|``pheno``|str#|Phenotype Code|.|

__example log file__
```
#########################################################################
# > .................. fam2dta
# > ........... author Richard Anney
# > ... report bugs to anneyr@cardiff.ac.uk
# >  code available at https://github.com/ricanney
# > .......... license CC BY-NC (CC Attribution-NonCommercial)
# > .................. all functions considered beta "caveat emptor"
#########################################################################
# Started:             30 Mar 2020 16:53:30 
# Username:            Richard Anney
# Operating System:    Windows
#########################################################################
# > ..... project name controls-combined-wtccc1
# > ... indiv. in file 3004
# > ....... saved data controls-combined-wtccc1_fam.dta
#########################################################################
# > overview of fam2dta
#########################################################################
# > 1. import fam file
# > 2. rename variable
#########################################################################

```