# recodegenotypes
__Author:__ Richard Anney
__Edited:__ 30 Mar 2020
__Current Version:__ v1

__Overview__
This function works directly within the code and converts variables ``a1`` and ``a2`` to IUPAC genotype codes ``M`` ``R`` ``W`` ``K`` ``Y`` ``S`` and non-standard ``ID`` for indels. Indels are crudely converted to ``I``and ``D`` coded alleles based on allele length using the code;
```
	noi di as text"# > .......... "as result"convert alleles to indel (I/D)"
	gen recodegenotype_counta1 = length(a1)
	gen recodegenotype_counta2 = length(a2)
	replace a1 = "I" if recodegenotype_counta1 > recodegenotype_counta2
	replace a2 = "I" if recodegenotype_counta2 > recodegenotype_counta1
	replace a1 = "D" if a2 == "I"
	replace a2 = "D" if a1 == "I"
	compress
```
Where non-standard (not defined as above) SNPs are observed an error code is created and two files are exported;
1. ``recodegenotype_error.exclude`` - contains list of ``snp`` containing errors
2. ``recodegenotype_error.update-alleles`` - contains a tab-delimited / no-header list of ``snp`` ``a1`` ``a2`` ``b1`` ``b2`` to update the errors. This error fix assumes that the allele codes differ at position 1. __this is a workaround for UKBB and is not universal - please check file before fixing using this file__ 

```
#example of recodegenotype_error.update-alleles
==========
rs2502408	CG	TG	C	T
rs7594227	GC	AC	G	A
rs9842095	AG	GG	A	G
rs7618587	TG	CG	T	C
rs1512006	AT	CT	A	C
rs4687002	TA	CA	T	C
rs1251144	ATT	CTT	A	C
rs2445871	GTATCT	ATATCT	G	A
==========
```

__Dependencies__
None

__Input Files__
None

__Syntax__
```
recodegenotypes
```

__Output Files__
None
This function creates the variable ``recodegenotype_gt``

| Variable Name | Storage Type | Variable Label | Notes |
|--------|--------|--------|--------|
|``recodegenotype_gt``|str2|Genotype (``RYMKSW(ID))``| NA


