# bim2ldexclude
__Author:__ Richard Anney
__Edited:__ 31 Mar 2020
__Current Version:__ v2

__Overview__
This function creates a ``plink`` formatted  ``*.exclude`` file of SNPs that fall within known regions of high linkage disequilibrium.
The function uses hard coded chromosome positions based on ``hg19|GRCH Build 37`` co-ordinates. The co-ordinates were extracted from 

1. Price et al. (2008) Long-Range LD Can Confound Genome Scans in Admixed Populations. Am. J. Hum. Genet. 86, 127-147
2. Weale M. (2010) Quality Control for Genome-Wide Association Studies from Michael R. Barnes and Gerome Breen (eds.), Genetic Variation: Methods and Protocols, Methods in Molecular Biology, vol. 628, DOI 10.1007/978-1-60327-367-1_19, © Springer Science+Business Media, LLC 2010

via https://genome.sph.umich.edu/wiki/Regions_of_high_linkage_disequilibrium_(LD):

Note that this can also be achieved using ``plink`` directly;

```
plink --file mydata --make-set high-ld.txt --write-set --out hild
plink --file mydata --exclude hild.set --recode --out mydatatrimmed
```

__Dependencies__
The function has the following dependencies are ;
``bim_path_to_short`` ``checkfile`` ``files2dta`` 

__Input Files__
For more details on input files see [input files](https://github.com/ricanney/stata/blob/master/input-files.md).
This function requires the ``*.bim`` file

__Syntax__

```
bim2ldexclude, bim(C:/data/plink_binary)
```
or
```
local bim C:/data/plink_binary
bim2ldexclude, bim(`bim')
```

__Output Files__

The resulting script creates two files;
1. ``*.exclude`` - list of excludable SNPs to use in a ``plink`` ``--exclude`` routine

__example log file__
```
#########################################################################
# > .................. bim2ldexclude
# > ........... author Richard Anney
# > ... report bugs to anneyr@cardiff.ac.uk
# >  code available at https://github.com/ricanney
# > .......... license CC BY-NC (CC Attribution-NonCommercial)
# > .................. all functions considered beta "caveat emptor"
#########################################################################
# Started:             31 Mar 2020 14:46:37 
# Username:            Richard Anney
# Operating System:    Windows
#########################################################################
# > ..... project name neonate-ahe-qc-v9
# > .... plink version H:/software/plink/bin/v1.9/plink_1.9b6.16-19feb2020.exe
# > ....... saved data neonate-ahe-qc-v9-bim2ldexclude.exclude
#########################################################################
# > overview of bim2ldexclude
#########################################################################
# > 1. import markers
# > 2. assume hg19 - drop regions of known high LD
#########################################################################
# References:
# > 1. Price et al. (2008) Long-Range LD Can Confound Genome Scans in Admixed
# >    Populations. Am. J. Hum. Genet. 86, 127-147
# > 2. Weale M. (2010) Quality Control for Genome-Wide Association Studies 
# >    from Michael R. Barnes and Gerome Breen (eds.), Genetic Variation: Methods
# >    and Protocols, Methods in Molecular Biology, vol. 628, 
# >    DOI 10.1007/978-1-60327-367-1_19, © Springer Science+Business Media, LLC 
# >    2010
# > from "https://genome.sph.umich.edu/wiki/Regions_of_high_linkage_disequilibrium_(LD)
#########################################################################

```