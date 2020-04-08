# create_1000genomes
__Author:__ Richard Anney
__Edited:__ 6 Apr 2020
__Current Version:__ v1

__Overview__
This function creates the ``all-1000g-phase3-chrall-mac5-v2`` reference genotype files from the 1000 genome project phase 3 ``*.vcf`` release.
The original ``*.vcf`` files are processed to a usable set of ``plink`` binaries. These ``plink`` binaries are essential in other applications such as ``graphlocuszoom``.

The function also creates;
1. ``all-1000g-phase3-chrall-mac5-v2.aims`` - list of the ancestry informative markers
2. ``all-1000g-phase3-chrall-mac5-v2.population`` - linker of identifiers to 1000 genomes population groups
3. ``all-1000g-phase3-chrall-mac5-v2.super-population`` - linker of identifiers to 1000 genomes super-population groups
4. ``gbr-1000g-phase3-chrall-mac5-v2.*`` - ``plink`` binaries for the ``GBR`` ancestry individuals only

The package runs on the linux and windows 10 (with ubuntu bash activates) operating systems. ``create_1000genomes`` uses the ``wget`` commands to download raw data. The raw data was originally released as part of the 1000 genomes (phase3) release at
http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/

__Dependencies__
The function has the following dependencies are ;
``bim2dta`` ``checkfile`` ``checkloc_name`` ``files2dta`` ``recodegenotype`` 

__Input Files__
None

__Syntax__
```
create_1000genomes
```

__Output Files__
1. ``all-1000g-phase3-chrall-mac5-v2.bim``
1. ``all-1000g-phase3-chrall-mac5-v2.bed``
1. ``all-1000g-phase3-chrall-mac5-v2.fam``
1. ``all-1000g-phase3-chrall-mac5-v2.log``
1. ``all-1000g-phase3-chrall-mac5-v2.aims``
1. ``all-1000g-phase3-chrall-mac5-v2.population``
1. ``all-1000g-phase3-chrall-mac5-v2.super-population``
1. ``all-1000g-phase3-chrall-mac5-v2-panel.dta``
1. ``gbr-1000g-phase3-chrall-mac5-v2.bim``
1. ``gbr-1000g-phase3-chrall-mac5-v2.bed``
1. ``gbr-1000g-phase3-chrall-mac5-v2.fam``
1. ``gbr-1000g-phase3-chrall-mac5-v2.log``

The file is deposited in the root directory in a sub-directory ``/create_1000genomes/ftp.1000genomes.ebi.ac.uk/``. 

__all-1000g-phase3-chrall-mac5-v2-panel.dta file format__

| Variable Name | Storage Type | Variable Label | Notes |
|--------|--------|--------|--------|
|fid|str7|Family Identifier|.|
|iid|str7|Individual Identifier|.|
|pop|str3|Population Code|.|
|super_pop|str3|Super Population Code|.|
|sex|byte|Family Identifier|labelled as; 1"Male" 2"Female" -9"Unknown"|






