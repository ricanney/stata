# snp2refid
__Author:__ Richard Anney
__Edited:__ 29 Apr 2020
__Current Version:__ v2

__Overview__
This function updates the naming of the ``snp`` variable to that of a reference file using the ``loc_name`` variable. The function requires an open file with ``snp`` ``chr`` and ``bp`` ``a1`` ``a2``. 

__Dependencies__
The function has the following dependencies are ;
``bim2dta`` ``checkfile`` ``checkloc_name`` ``ref_path_to_short`` 

__Input Files__
NA

_Reference data files_
The function utilises a pre-prepared file that are created via ``create_hrc``

__Syntax__

```
global snp2refid  ${init_root}/create_hrc/HRC.r1-1.GRCh37.wgs.mac5/all-hrc-1.1-chrall-mac5
use H:\software\stata\data\summaryqc\cardiff-alspac-genetics-2020-s25BAARS\cardiff-alspac-genetics-2020-s25BAARS-summaryqc.dta, clear
snp2refid, ref(${snp2refid})
```

