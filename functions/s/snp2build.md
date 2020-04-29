# snp2build
__Author:__ Richard Anney
__Edited:__ 29 Apr 2020
__Current Version:__ v2

__Overview__
This function assesses the build of an open file with a list of ``snp`` ``chr`` and ``bp``. The output is a global term ${snp2build_build}

__Dependencies__
The function has the following dependencies are ;
``checkfile`` ``checktabbed`` ``files2dta``  ``ref_path_to_short`` 

__Input Files__
NA

_Reference data files_
The function utilises a pre-prepared file ``bim2build.dta`` - available on request. 

_to do:_ create ``create_bim2build`` function to create file on the fly

__Syntax__

```
global bim2build  H:/databank/3-packages/bim2build/bim2build.dta
use H:\software\stata\data\summaryqc\test\test-summaryqc.dta, clear
snp2build, ref(${bim2build})
```

