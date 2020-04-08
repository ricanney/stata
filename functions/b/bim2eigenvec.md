# bim2eigenvec
__Author:__ Richard Anney
__Edited:__ 31 Mar 2020
__Current Version:__ v2

__Overview__
This function runs a set of modules to calculate the top 10 principal components (PC's) based on ancestry for a set of ``plink`` binaries. 

__Dependencies__
The function has the following dependencies are ;
``bim_path_to_short`` ``bim2count`` ``bim2ldexclude``  ``checkfile`` ``checktabbed``  ``files2dta`` 

__Non Standard Stata Packages__
``dropmiss`` - to install type;

```
net install dm89_1.pkg
```

__Input Files__
For more details on input files see [input files](https://github.com/ricanney/stata/blob/master/input-files.md).
This function requires the ``*.bim`` file

__Syntax__

```
bim2eigenvec, bim(C:/data/plink_binary)
```
or
```
local bim C:/data/plink_binary
bim2eigenvec, bim(`bim') pc(5)
```

__Output Files__

The resulting script creates two files;
1. ``*_eigenval.dta`` 
2. ``*_eigenvec.dta``

The file is deposited in the root directory in a sub-directory ``/bim2eigenvec/<filename>``. 
For ``plink_binary.bim`` the subdirectory will be ``<root>/bim2eigenvec/plink_binary``.

__*_eigenval.dta file format__

| Variable Name | Storage Type | Variable Label | Notes |
|--------|--------|--------|--------|
|``pc``|double|Principal Component| .|
|``eigenval``|double|Eigenvalue|.|



__*_eigenev.dta file format__

| Variable Name | Storage Type | Variable Label | Notes |
|--------|--------|--------|--------|
|``fid``|str12|Family Identifier| .|
|``iid``|str12|Individual Identifier|.|
|``pc1``|double|Ancestry Principal Component 1 |.
|...|double| |.
|``pc10``|double|Ancestry Principal Component 10 |.


__example log file__
```
#########################################################################
# > .................. bim2eigenvec
# > ........... author Richard Anney
# > ... report bugs to anneyr@cardiff.ac.uk
# >  code available at https://github.com/ricanney
# > .......... license CC BY-NC (CC Attribution-NonCommercial)
# > .................. all functions considered beta "caveat emptor"
#########################################################################
# Started:             31 Mar 2020 11:11:40 
# Username:            Richard Anney
# Operating System:    Windows
#########################################################################
# > ..... project name controls-combined-wtccc1
# > .... plink version H:/software/plink/bin/v1.9/plink_1.9b6.16-19feb2020.exe
# > .... plink version H:/software/plink/bin/v2.0/plink_2_dev2-20mar2017.exe
# > ............ N PCs 10
# > .... indiv. in bim 3004
# > indiv. in eigenvec 3004
# > ....... saved data controls-combined-wtccc1_eigenvec.dta
# > ....... saved data controls-combined-wtccc1_eigenval.dta
#########################################################################
# > overview of bim2eigenvec
#########################################################################
# > 1. define ld pruned subset 
# > 2. calculate pc
# > 3. process files to eigenvec and eigenval dta files
#########################################################################

```