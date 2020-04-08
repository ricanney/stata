# bim_path_to_short
__Author:__ Richard Anney
__Edited:__ 30 Mar 2020
__Current Version:__ v1

__Overview__
This function creates a global ``${bim_short}``, containing the filename,from the full path of the ``*.bim`` file.
For example; ``C:/data/plink_binary`` would be converted to ``plink_binary``

__Dependencies__
None

__Input Files__
None

__Syntax__
```
bim_path_to_short, bim(C:/data/plink_binary)
```
or
```
local bim C:/data/plink_binary
bim_path_to_short, bim(`bim')
```

You define the path of the plink binaries (with common name) in the script using ``bim()`` - however, you ommit the ``*.bim`` from the file name. The above example is for the file path ``C:/data/plink_binary.bim``

__Output Files__
None
This function creates the global link ``${bim_short}``

