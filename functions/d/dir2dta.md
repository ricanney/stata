# files2dta
__Author:__ Richard Anney
__Edited:__ 31 Mar 2020
__Current Version:__ v1

__Overview__
This function creates a ``*.dta`` listing the files in a directory

__Dependencies__
None

__Input Files__
None

__Syntax__

```
files2dta, dir(H:\software\stata\functions\f)
```

__Output Files__

The function creates a single temporary file called ``_files2dta.dta`` containing the filenames listed in the variable ``files``



___files2dta.dta file format__

| Variable Name | Storage Type | Variable Label | Notes |
|--------|--------|--------|--------|
|``files``|str#|Files Listed in Directorye|.|



__example log file__
This function does not create a log file
