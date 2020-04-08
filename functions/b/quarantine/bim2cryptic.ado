*! 1.0.1 Richard Anney 21dec2018

program  bim2cryptic
local    function bim2cryptic
local    `function'_version v2
syntax , bim(string asis)  

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
	bim_path_to_short, bim(`bim')
	local `function'_bim_short ${bim_short}
	cd ${init_root}	
	!mkdir `function'
	cd `function'
	noi di as text"# > ..... project name "as result"``function'_bim_short'"
	!mkdir ``function'_bim_short'
	cd ``function'_bim_short'
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 2 - define files
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 2: define files"
	noi di as text"#########################################################################"
	noi checkfile, file(${plink2})
	foreach file in bim bed fam {
		noi checkfile, file(`bim'.`file')
		}
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 3 - define kinship matrix
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3: define kinship matrix"
	noi di as text"#########################################################################"		
	bim2count, bim(`bim')
	noi di as text"# > ..... SNPs in file "as result "${bim2count_snp}"
	noi di as text"# > .... indiv in file "as result "${bim2count_ind}"
	noi di as text"# > .................. create kinship matrix "as result "${bim2count_ind}"as text" x "as result "${bim2count_ind}"
	noi di as text"# > .................. define ld pruned subset"
	local already_created ${init_root}/bim2ld_subset/``function'_bim_short'/``function'_bim_short'-bim2ld_subset.extract
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		}
	else {
		noi di as text"# > .................. create "as result "``function'_bim_short'-bim2ld_subset.extract"
		noi bim2ld_subset, bim(`bim')
		noi checkfile, file(`already_created')
		}
	noi di as text"# > .................. create kinship square"
	!$plink2 --bfile `bim' --extract `already_created' --make-king square --out ``function'_bim_short'
	noi di as text"# > .................. importing kinship matrix"
	import delim using ``function'_bim_short'.king, clear case(lower)
	count
	forvalues i=1/ `r(N)' {
		replace v`i' = . in `i'
		replace v`i' = 0 if v`i' < 0
		}
	gen obs = _n
	aorder
	dropmiss, force
	save ``function'_bim_short'-king.dta, replace
	import delim using ``function'_bim_short'.king.id, clear case(lower)
	rename (v1 v2) (fid iid)
	for var fid iid: tostring X, replace
	gen obs = _n
	aorder
	merge 1:1 obs using ``function'_bim_short'-king.dta, update
	drop _m
	qui { // blank out known family identifiers
		noi di as text"# > .................. blank out know related individuals (via. FID)"
		gen a = ""
		tostring obs, replace
		replace a = "replace v" + obs + `" = . if fid ==""' + fid + `"""'
		outsheet a using tmp.do, non noq replace
		do tmp.do
		erase tmp.do
		drop a obs
		}
	qui { // calculate by-individual metrics
		noi di as text"# > .................. calculate by-individual metrics"
		count
		egen rm = rowmean(v1-v`r(N)')
		count
		egen rx = rowmax(v1-v`r(N)')	
		keep fid iid rm
		gen xs_relate = .
		sum rm
		foreach i of num 1/5 {
			sum rm
			replace xs_relate = `i' if rm > `r(mean)' + (`i' * `r(sd)')
			}
		}
	qui { // define excessive relatedness
		noi di as text"# > .................. excessive relatedness defined as "as result"3 "as text"x population standard deviation"
		count if xs > 3 & xs != .
		local `function'_N `r(N)'

		noi di as text"# > .................. "as result `r(N)' as text" individuals defined as having excessive relatedness"
		if `r(N)' != 0 {
			outsheet fid iid if xs > 3 & xs != . using ``function'_bim_short'-`function'.remove, replace non noq
			}
		else {
			clear
			set obs 1
			outsheet using ``function'_bim_short'-`function'.remove, replace non noq
			}
		noi di as text"# > ......... saved to "as result"``function'_bim_short'-function.remove"
		}
	}
qui { // module 4 - write log file
	file open myfile using "``function'_bim_short'-`function'.log", write replace
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
	file write myfile `"# > ..... project name ``function'_bim_short'"' _n
	file write myfile `"# > .... plink version ${plink}"' _n
	file write myfile `"# > .... plink version ${plink2}"' _n
	file write myfile `"# > .... indiv. in bim ${bim2count_ind}"' _n
	file write myfile `"# > .................. ``function'_N' defined as having excessive relatedness"' _n
	file write myfile `"# > ....... saved data ``function'_bim_short'-`function'.remove"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > overview of `function'"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > 1. define ld pruned subset "' _n
	file write myfile `"# > 2. create square relatedness matrix "' _n
	file write myfile `"# > 3. exclude calculations with know relatives (i.e. common FID)"' _n
	file write myfile `"# > 3. calculate mean relatedness"' _n
	file write myfile `"# > 4. identify individuals with > 3x sd of the mean"' _n
	file write myfile `"#########################################################################"' _n
	file close myfile	
	}
qui { // module 5 - clean folder
	files2dta, dir(`c(pwd)')
	gen keep = .
	replace keep = 1 if file == "``function'_bim_short'-`function'.remove"
	replace keep = 1 if file == "``function'_bim_short'-`function'.log"
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
