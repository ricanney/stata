*! 1.0.1 Richard Anney 22jan2020

program  summaryqc2eqtl
local    function summaryqc2eqtl
local    `function'_version v1

syntax , summaryqc(string asis) ldref(string asis) glist_hg19(string asis) eqtl_list(string asis) [minp(real 1e-5) mineqtl(real 5e-8)]

noi di as text""
noi di as text"#########################################################################"
noi di as text"# > .................. `function'"
noi di as text"# > ........... author "as result"Richard Anney"
noi di as text"# > ....report bugs to "as result"anneyr@cardiff.ac.uk"
noi di as text"# >  code available at "as result"https://github.com/ricanney"
noi di as text"# > .......... license "as result"CC BY-NC (CC Attribution-NonCommercial)"
noi di as text"#########################################################################"
noi di as text"# Started: $S_DATE $S_TIME"
noi di as text"# Username:            `c(username)'"
noi di as text"# Version:             "as result"``function'_version'"
noi di as text"# Operating System:    `c(os)'"
noi di as text"#########################################################################"
noi di as text""	

preserve
local return `c(pwd)'
qui { // module 1 - prepare folders 
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 1: prepare folders"
	noi di as text"#########################################################################"
	ref_path_to_short, ref(`ldref')
	local `function'_ldref_short ${ref_short}
	ref_path_to_short, ref(`eqtl_list')
	local `function'_elist_short ${ref_short}
	ref_path_to_short, ref(`glist_hg19')
	local `function'_glist_short ${ref_short}
	ref_path_to_short, ref(`summaryqc')
	local `function'_summaryqc_short ${ref_short}
	cd ${init_root}	
	!mkdir `function'
	cd `function'
	noi di as text"# > ..... project name "as result"``function'_summaryqc_short'"
	!mkdir ``function'_summaryqc_short'
	cd ``function'_summaryqc_short'
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 2 - define files
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 2: define files"
	noi di as text"#########################################################################"
	noi checkfile, file(${plink})
	noi checkfile, file(${plink2})
	checktabbed
	noi di as text"# > ......... locating files to process"
	noi checkfile, file(`summaryqc'-summaryqc.dta)
	foreach file in bim bed fam {
		noi checkfile, file(`ldref'.`file')
		}
	noi checkfile, file(`eqtl_list')
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 3 - perform data clumping
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3: load top hits"
	noi di as text"#########################################################################"
	noi di as text"# > ..... loading from "as result"``function'_summaryqc_short'-summaryqc2top.dta"
	noi di as text"# > ..... limit to P < "as result"`minp'"
	local already_created ${init_root}/summaryqc2top/``function'_summaryqc_short'/``function'_summaryqc_short'-summaryqc2top.dta
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		}
	else {	
		noi summaryqc2top, summaryqc(`summaryqc') ldref(`ldref') glist_hg19(`glist_hg19')
		noi checkfile, file(`already_created')
		}
	noi di as text"# > define LD pairs of top hits r^2"
	noi di as text"# > ..... ld reference "as result"``function'_ldref_short'"
	local already_created ${init_root}/`function'/``function'_summaryqc_short'/``function'_summaryqc_short'-ld-pairs.dta
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		}
	else {
		use 	${init_root}/summaryqc2top/``function'_summaryqc_short'/``function'_summaryqc_short'-summaryqc2top.dta, clear
		keep SNP
		outsheet using ``function'_summaryqc_short'.snp-in, non noq replace
		noi di as text"# > mapping LD-pairs using plink"
		!plink2 --bfile `ldref' --r2 --ld-window 100000 --ld-window-kb 100000 --ld-snp-list ``function'_summaryqc_short'.snp-in --ld-window-r2 0.8 --out ``function'_summaryqc_short'.snp-out
		!$tabbed ``function'_summaryqc_short'.snp-out.ld
		import delim using ``function'_summaryqc_short'.snp-out.ld.tabbed, clear varnames(1)
		rename (snp_a snp_b) (SNP PROXY_SNP)
		for var chr_b bp_b: tostring X, replace
		gen PROXY_POSITION = "chr" + chr_b + ":" + bp_b
		rename r2 R2
		keep SNP PROXY_SNP PROXY_POSITION R2
		save ${init_root}/`function'/``function'_summaryqc_short'/``function'_summaryqc_short'-ld-pairs.dta, replace
		noi checkfile, file(`already_created')
		}
	noi di as text"# >  map cis-eQTL from "as result"``function'_elist_short'"
	local already_created ${init_root}/`function'/``function'_summaryqc_short'/``function'_summaryqc_short'-`function'.dta
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		}
	else {
		use if pval < 5e-8 using `eqtl_list', clear
		renvars, upper
		rename SNP PROXY_SNP
		keep PROXY_SNP SYMBOL TISSUE SUPER_TISSUE
		merge m:m PROXY_SNP using ``function'_summaryqc_short'-ld-pairs.dta
		keep if _m == 3
		sort SNP SYMBOL SUPER TISSUE
		keep SNP SYMBOL SUPER TISSUE
		duplicates drop
		egen tissue_count = seq(),by(SNP SYMBOL)
		gsort -tissue_count
		egen tissue_obs = seq(),by(SNP SYMBOL)
		keep if tissue_obs == 1
		gen count = tissue_count
		tostring count, replace
		gen EQTL = ""
		replace EQTL = "(multi)" if tissue_count > 1
		replace EQTL = "(fibroblast)" if TISSUE == "cells_transformed_fibroblasts" & EQTL == ""
		replace EQTL = "(cerebellum)" if TISSUE == "brain_cerebellum" & EQTL == ""
		replace EQTL = "(" + SUPER + ")" if EQTL == "" 
		keep SNP SYMBOL EQTL
		replace EQTL = SYMBOL + EQTL
		keep SNP EQTL
		sort SNP EQTL
		order SNP 
		egen x = seq(),by(SNP)
		reshape wide EQTL, i(SNP) j(x)
		gen EQTL99 = ""
		gen EQTL = ""
		for var EQTL1-EQTL99: replace EQTL = EQTL + " " +  X
		foreach i of num 1/99 {
			replace EQTL = subinstr(EQTL,"  "," ",.)
			}
		keep SNP EQTL
		merge 1:1 SNP using ${init_root}/summaryqc2top/``function'_summaryqc_short'/``function'_summaryqc_short'-summaryqc2top.dta
		drop _m
		order SNP POS A1 A2 OR P LD_RA EQTL
		replace EQTL = "NONE" if EQTL == ""
		sort P
		save ``function'_summaryqc_short'-`function'.dta,replace
		capture confirm file `already_created'
		}
	local already_created ${init_root}/`function'/``function'_summaryqc_short'/``function'_summaryqc_short'-`function'-detailed.dta
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		}
	else {
		use if pval < 5e-8 using `eqtl_list', clear
		renvars, upper
		rename SNP PROXY_SNP
		keep PROXY_SNP SYMBOL TISSUE SUPER_TISSUE
		merge m:m PROXY_SNP using ``function'_summaryqc_short'-ld-pairs.dta
		keep if _m == 3
		qui { // define proxy SNPs
			keep SNP SYMBOL SUPER TISSUE PROXY_SNP R2
			gen PROXY = "*" if PROXY_SNP != SNP
			sort SNP SYMBOL SUPER TISSUE 
			keep SNP SYMBOL SUPER TISSUE PROXY PROXY_SNP R2
			gen TOP = ""
			replace TOP = "Y" if SNP == PROXY_SNP
			sort TOP R2
			egen x = seq(),by(TISSUE SYMBOL SNP)
			keep if x == 1
			keep SNP SYMBOL TISSUE SUPER PROXY_SNP R2
			order SNP SYMBOL TISSUE SUPER PROXY_SNP R2
			sort SNP SYMBOL TISSUE SUPER PROXY_SNP R2
			}
		merge m:1 SNP using ${init_root}/summaryqc2top/``function'_summaryqc_short'/``function'_summaryqc_short'-summaryqc2top.dta
		drop _m
		for var SYMBOL PROXY : replace X = "NONE" if X == ""
		for var TISSUE SUPER : replace X = "NA" if X == ""
		sort P SNP SUPER TISSUE SYMBOL
		save ``function'_summaryqc_short'-`function'-detailed.dta,replace
		capture confirm file `already_created'
		}
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 4 - write log file
	file open myfile using "``function'_summaryqc_short'-`function'.log", write replace
	file write myfile`"#########################################################################"' _n
	file write myfile`"# > .................. `function'"' _n
	file write myfile`"# > ........... author Richard Anney"' _n
	file write myfile`"# > ... report bugs to anneyr@cardiff.ac.uk"' _n
	file write myfile`"# >  code available at https://github.com/ricanney"' _n
	file write myfile`"# > .......... license CC BY-NC (CC Attribution-NonCommercial)"' _n
	file write myfile`"# > .................. all functions considered beta "caveat emptor""' _n
	file write myfile`"#########################################################################"' _n
	file write myfile`"# Started:             $S_DATE $S_TIME "' _n
	file write myfile`"# Username:            `c(username)'"' _n
	file write myfile`"# Operating System:    `c(os)'"' _n
	file write myfile`"#########################################################################"' _n
	file write myfile `"# > ....... summaryqc ``function'_summaryqc_short'"' _n
	file write myfile `"# > ....... reference ``function'_ldref_short'"' _n
	file write myfile `"# > ....... gene-list ``function'_glist_hg19_short'"' _n
	file write myfile `"# > .. eQTL reference ``function'_elist_short'"' _n
	file write myfile `"# > .. min-P reported `minp'"' _n
	file write myfile `"# > ...... saved data ``function'_summaryqc_short'-`function'.dta"' _n
	file write myfile `"# > ...... saved data ``function'_summaryqc_short'-`function'-all.dta"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > overview of `function'"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > 1. identify top N SNPs"' _n
	file write myfile `"# > 2. extend LD to capture proxy_SNPs"' _n
	file write myfile `"# > 3. map eQTL pairs (abridged to *.dta)"' _n
	file write myfile `"# > 4. map eQTL pairs (detailed to detailed.dta"' _n
	file write myfile `"#########################################################################"' _n
	file close myfile	
	}
qui { // module 5 - clean folder
	files2dta, dir(`c(pwd)')
	gen keep = .
	replace keep = 1 if file == "``function'_summaryqc_short'-`function'.log"
	replace keep = 1 if file == "``function'_summaryqc_short'-`function'.dta"
	replace keep = 1 if file == "``function'_summaryqc_short'-`function'-detailed.dta"
	drop if keep == 1
	gen script = "erase " + file
	outsheet script using temp.do, non noq replace 
	do temp.do
	erase temp.do
	clear
	}
qui cd `return'
restore
noi di as text"#########################################################################"
noi di as text"# Completed: $S_DATE $S_TIME"
noi di as text"#########################################################################"
end;

   
