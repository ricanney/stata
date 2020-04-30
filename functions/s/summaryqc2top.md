# summaryqc2top
__Author:__ Richard Anney
__Edited:__ 30 Apr 2020
__Current Version:__ v2

__Overview__
This function processes the GWAS summary statistics in  ``summaryqc`` and limits to the top LD independent markers (P < 1e-5; adjustable). LD independent markers are identified through a ``PLINK`` ``--clump`` routine. Summary data is reported as ``summaryqc`` with additional variables ``ld_range`` indicating the boundaries of the ld-pruning.

__Dependencies__
The function has the following dependencies are ;
``bim2dta`` ``checkfile`` ``checktabbed`` ``files2dta``  ``ref_path_to_short`` 

__Input Files__
1. ``summaryqc`` formated input file
2. ``plink`` binaries from an ancestrally related reference set to determine the linkage disequilibrium structure in the region - 1000 genomes reference can be created via ``create_1000genomes``.

__Syntax__
```
global ldref      ${init_root}/create_1000genomes/ftp.1000genomes.ebi.ac.uk/gbr-1000g-phase3-chrall-mac5-v2
global summaryqc  ${init_root}/summaryqc/cardiff-alspac-genetics-2020-p25DAWBA/cardiff-alspac-genetics-2020-p25DAWBA
summaryqc2top, summaryqc(${summaryqc}) ldref(${ldref})
```

__Output Files__
The resulting file is created with the ``<summaryqc>-summaryqc2top.dta`` suffix 
The file is deposited in the root directory in a sub-directory ``/summaryqc2top/<summaryqc>``. 

| Variable Name | Storage Type | Variable Label | Notes |
|--------|--------|--------|--------|
|``snp``  |str |Marker Name| |
|``position``  |str |Chromosome: Physical Location
|``a1``   |str |Allele 1 (``ACGT(D/I)``)| convert long indel naming to ``I`` or ``D``
|``a2``   |str |Allele 2 (``ACGT(D/I)``)| convert long indel naming to ``I`` or ``D``
|``beta`` |double | Beta Coefficient| Beta coefficient from the GWAS 
|``se``   |double | Standard Error | Standard Error from the GWAS 
|``z``    |double | Z-statistic | Standardised Z-score  ``beta/se``
|``or``   |double | Odds Ratio |.
|``p``    |double | P-value|.
|``ld_range``|str|Clumped Region|.


__example log file__
```
#########################################################################
# > .................. summaryqc2top
# > ........... author Richard Anney
# > ... report bugs to anneyr@cardiff.ac.uk
# >  code available at https://github.com/ricanney
# > .......... license CC BY-NC (CC Attribution-NonCommercial)
# > .................. all functions considered beta "caveat emptor"
#########################################################################
# Started:             30 Apr 2020 14:33:42 
# Username:            Richard Anney
# Operating System:    Windows
#########################################################################
# > ....... summaryqc cardiff-alspac-genetics-2020-p25DAWBA
# > ....... reference gbr-1000g-phase3-chrall-mac5-v2
# > .. min-P reported .00001
# > ...... saved data cardiff-alspac-genetics-2020-p25DAWBA-summaryqc2top.dta
#########################################################################
# > overview of summaryqc2top
#########################################################################
# > 1. import top N SNPs
# > 2. clump versus reference
# > 3. prepare data for publication
#########################################################################
```