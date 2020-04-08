# create_recombref
__Author:__ Richard Anney
__Edited:__ 2 Apr 2020
__Current Version:__ v1

__Overview__
This function creates the ``genetic_map_HapMapII_GRCh37-recombination-rate.dta``  files for use in functions ``graphlocuszoom``. The package runs on the linux and windows 10 (with ubuntu bash activates) operating systems. ``create_recombref`` uses the ``wget`` and ``tar`` commands to download and unpack the raw data. The raw data was originally released as part of the HapMap (phase2) release at ftp://ftp.ncbi.nlm.nih.gov/hapmap/recombination/2011-01_phaseII_B37/ 

__Dependencies__
The function has the following dependencies are ;
``checkfile``

__Input Files__
None

__Syntax__
```
create_recombref
```

__Output Files__
``genetic_map_HapMapII_GRCh37-recombination-rate.dta``

The file is deposited in the root directory in a sub-directory ``/create_recombref/genetic_map_HapMapII_GRCh37``. 

__Homo_sapiens.GRCh37.`version'.gtf_gene.dta file format__

| Variable Name | Storage Type | Variable Label | Notes |
|--------|--------|--------|--------|
|chr|byte|chromosome|.|
|bp|long|Position(bp)|.
|recomb_rate|double|Rate(cM/Mb)|.





