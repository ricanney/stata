# Stata Cheat Sheet
__Author:__ Richard Anney
__Edited:__ 30 Mar 2020

[TOC]

## defining macros from a cell
To define a local macro ``new_macro`` from observation ``1`` for the variable ``var1`` use the following code;

```
local new_macro = var1[1]
```

To define a global macro ``new_macro`` from observation ``1`` for the variable ``var1`` use the following code;

```
global new_macro = var1[1]
```


