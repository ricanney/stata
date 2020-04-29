# summaryqc
__Author:__ Richard Anney
__Edited:__ 29 Apr 2020
__Current Version:__ v3

__Overview__
This function creates processes GWAS summary statistics. The function;
1. removes SNPs with out of bounds P-values
2. removes SNP duplicates (by identifier)
3. removes SNPs with ``info`` > 2 (out of bounds) or ``info`` < .8
4. removes SNPs from meta-analysis where they were absent from >2 studies - via ``direction`` variable
5. check build is ``hg19+1`` using ``snp2build``
6. update name to ``hrc`` using ``snp2refid``
7. plot manhattan via ``graphmanhattan``

There are some pre-processing steps to prepare the file for ``summaryqc``. Specifically, the input file need to be imported into ``Stata`` and the following variables are needed to be generated;

_core_  ``chr`` ``bp`` ``snp`` ``a1`` ``a2`` ``beta`` ``se`` ``z`` ``or`` ``l95 `` ``u95`` ``p`` ``n`` 
_optional_ ``info`` INFO score / Imputation metric
_optional_ ``direction`` METAL direction / Measure of imputation direction as (``+-?``)

__Deriving missing variables__
Where the variable is missing, we can crudely derive the missing variable. _This is not ideal_, especially where ``se`` is missing. These variables are often test specific and therefore are a crude proxy.

```
* missing se
gen sd = sqrt(n) * p * (1-p)
gen mean = n * p
gen se = sd / mean

* missing beta
gen beta = ln(or)

* missing or
gen or = exp(beta)
gen l95 = or - (1.95 * se)
gen u95 = or + (1.95 * se)

* missing z
gen z = beta / se
```

__Dependencies__
The function has the following dependencies are ;
``checkfile`` ``checkloc_name`` ``files2dta`` ``graphmanhattan`` ``ref_path_to_short`` ``snp2build`` ``snp2refid``

__Input Files__
None

__Syntax__
```
global bim2build  H:/databank/3-packages/bim2build/bim2build.dta
global hrcref     ${init_root}/create_hrc/HRC.r1-1.GRCh37.wgs.mac5/all-hrc-1.1-chrall-mac5

qui { // bespoke pre-process dataset
	use H:\software\stata\data\bim2dta\alspac-ext-b2661-2-1000g-qc-v9-25yr-subset\alspac-ext-b2661-2-1000g-qc-v9-25yr-subset_bim.dta, clear
	keep snp a1 a2
	rename (a1 a2) (b1 b2)
	merge 1:1 snp using alspac-ext-b2661-2-1000g-qc-v9-25yr-subset-s25BAARS.dta
	gen a2 = b2
	replace a2 = b1 if a1 == a2
	drop b1 b2 _m
	gen n = 3122
	tostring bp, replace
	tostring chr, replace
	gen or = exp(beta)
	gen l95  = or - (1.96 * se)
	gen u95  = or + (1.96 * se)
	drop if chr == "23"
	drop if p == .
	}
noi summaryqc, input(alspac-ext-b2661-2-1000g-qc-v9-25yr-subset-s25BAARS) ///
               output(cardiff-alspac-genetics-2020-s25BAARS) ///
               id_ref(${hrcref}) ///
               build_ref(${bim2build})

```

__Output Files__
The resulting file is created with the ```output'-summaryqc.dta`` suffix 
The file is deposited in the root directory in a sub-directory ``/summaryqc/`output'``. 

| Variable Name | Storage Type | Variable Label | Notes |
|--------|--------|--------|--------|
|``chr``  |str |Chromosome|Chromosomes are numbered 1-25|
|``snp``  |str |Marker Name| |
|``bp``   |str |Physical Location| |
|``a1``   |str |Allele 1 (``ACGT(D/I)``)| convert long indel naming to ``I`` or ``D``
|``a2``   |str |Allele 2 (``ACGT(D/I)``)| convert long indel naming to ``I`` or ``D``
|``beta`` |double | Beta Coefficient| Beta coefficient from the GWAS 
|``se``   |double | Standard Error | Standard Error from the GWAS 
|``z``    |double | Z-statistic | Standardised Z-score  ``beta/se``
|``or``   |double | Odds Ratio |.
|``l95``  |double | 95% Confidence Interval (lower)|.
|``u95``  |double | 95% Confidence Interval (upper)|.
|``p``    |double | P-value|.
|``n``    |double | Sample Size|.
