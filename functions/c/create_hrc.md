# create_hrc
__Author:__ Richard Anney
__Edited:__ 6 Apr 2020
__Current Version:__ v1

__Overview__
This function creates a ``*.bim`` file for the ``HRC.r1-1.GRCh37.wgs.mac5`` data. The resulting dataset is renamed ``all-hrc-1.1-chrall-mac5``

The original ``*.vcf`` files do not contain individual genotype information, but a ``*.bim`` file can be created and used other applications such as ``genotypeqc``. The main use of this file is for data harmonisation.

The package runs on the linux and windows 10 (with ubuntu bash activates) operating systems. ``create_hrc`` uses the ``wget`` ``zcat`` ``cut``  commands to download and process raw data. The raw data was downloaded from ``ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.vcf.gz``

__Dependencies__
The function has the following dependencies are ;
``bim2dta`` ``checkfile`` ``checkloc_name`` ``files2dta`` ``recodegenotype`` 

__Input Files__
None

__Syntax__
```
create_hrc
```

__Output Files__
1. ``HRC.r1-1.GRCh37.wgs.mac5.sites.core.chrall.dta``
1. ``all-hrc-1.1-chrall-mac5.bim``
1. ``\bim2dta\all-hrc-1.1-chrall-mac5_bim.dta``
1. ``\bim2dta\all-hrc-1.1-chrall-mac5_bim_noALLELE.dta``

The file is deposited in the root directory in a sub-directory ``/create_hrc/HRC.r1-1.GRCh37.wgs.mac5/``. 

__HRC.r1-1.GRCh37.wgs.mac5.sites.core.chrall.dta file format__

| Variable Name | Storage Type | Variable Label | Notes |
|--------|--------|--------|--------|
|``snp``|str10|Marker Name| |
|``chr``|str2|Chromosome|Chromosomes are numbered 1-25|
|``bp``|str9|Physical Location| |
|``a2``|str1|Allele 2 (``ACGT(D/I)``)| convert long indel naming to ``I`` or ``D``
|``a1``|str1|Allele 1 (``ACGT(D/I)``)| convert long indel naming to ``I`` or ``D``
|``af``|double|Allele Frequency (``a1``)|.
|``gt``|str2|Genotype (``RYMKSW(ID))``| from ``recodegenotypes``
|``loc_name``|str18|Marker loc_name (``chr<#>:<pos>-<gt>``)|from ``checkloc_name``






