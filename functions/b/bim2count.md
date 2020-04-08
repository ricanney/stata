# bim2count
__Author:__ Richard Anney
__Edited:__ 27 Mar 2020
__Current Version:__ v3

__Overview__
This function counts the number of SNPs and individuals present in a plink binary ``.bim`` and ``.fam`` file.

__Dependencies__
The function has the following dependencies are ;
``bim_path_to_short`` ``checkfile`` 

__Input Files__
For more details on input files see [input files](https://github.com/ricanney/stata/blob/master/input-files.md).
This function requires the ``*.bim`` file

__Syntax__
```
bim2count, bim(C:/data/plink_binary)
```
or
```
local bim C:/data/plink_binary
bim2count, bim(`bim')
```

You define the ``.bim`` file in the script using ``bim()`` - however, you ommit the ``*.bim`` from the file name. The above example is for the file path ``C:/data/plink_binary.bim``

__Output Files__
This function does not create an output file. The funcion creates two ``global`` links;
``{$bim2count_snp}`` = number of SNPs in the ``*.bim`` file
``{$bim2count_ind}`` = number of individuals in the ``*.fam`` file

__example log file__
The log file is deposited in the root directory in a sub-directory ``/bim2count/<filename>``. 
For ``plink_binary.bim`` the subdirectory will be ``<root>/bim2count/plink_binary``.

```
#########################################################################
# > .................. bim2count
# > ........... author Richard Anney
# > ... report bugs to anneyr@cardiff.ac.uk
# >  code available at https://github.com/ricanney
# > .......... license CC BY-NC (CC Attribution-NonCommercial)
# > .................. all functions considered beta "caveat emptor"
#########################################################################
# Started:             27 Mar 2020 16:05:06 
# Username:            Richard Anney
# Operating System:    Windows
#########################################################################
# > ..... project name hapmap3_ceu_tsi
# > ..... SNPs in file 1464434
# > .... indiv in file 267
#########################################################################
```
