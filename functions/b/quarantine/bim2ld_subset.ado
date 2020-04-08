*! 1.0.2 Richard Anney 24jan2019

program  bim2ld_subset
local    function bim2ld_subset
local    `function'_version v2
syntax , bim(string asis) [n(real 50000)]

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
	foreach file in bim bed fam {
		noi checkfile, file(`bim'.`file')
		}
	noi checkfile, file(${plink})
	noi checkfile, file(${plink2})
	checktabbed
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 3 - extract subset of ld independent SNPs
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3: extract subset of ld independent SNPs"
	noi di as text"#########################################################################"
	noi di as text"# > .................. apply bim2ldexclude"
	noi di as text"# > .................. exclude high ld regions"
	local already_created ${init_root}/bim2ldexclude/``function'_bim_short'/``function'_bim_short'-bim2ldexclude.exclude
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		}
	else {
		noi di as text"# > .................. create "as result "``function'_bim_short'-bim2ldexclude.exclude"
		noi bim2ldexclude, bim(`bim')
		noi checkfile, file(`already_created')
		}
	!$plink  --bfile `bim' --maf 0.05 --exclude `already_created' --make-bed --out ``function'_bim_short'-`function'
	noi di as text"# > .................. ld prune data (window=1000, slide=5, r2=.2)"
	!$plink  --bfile ``function'_bim_short'-`function' --indep-pairwise 1000 5 0.2  --out ``function'_bim_short'-`function'
	import delim using ``function'_bim_short'-`function'.prune.in, clear 
	gen x = uniform()
	sort x
	drop if _n > `n'
	noi di as text"# > .................. limit to maximum of `n' SNPS"
	outsheet v1 using ``function'_bim_short'-`function'.extract, non noq replace
	count
	local `function'_N `r(N)'
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
	file write myfile `"# > ... SNPs in subset ``function'_N'"' _n
	file write myfile `"# > ....... saved data ``function'_bim_short'-`function'.extract"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > overview of `function'"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > 1. limit bim to maf .05"' _n
	file write myfile `"# > 2. exclude high ld markers"' _n
	file write myfile `"# > 3. perform ld pruning (window=1000, slide=5, r2=.2)"' _n
	file write myfile `"# > 4. limit bim to `n' SNPs"' _n
	file write myfile `"#########################################################################"' _n
	file close myfile	
	}
qui { // module 8 - clean folder
	files2dta, dir(`c(pwd)')
	gen keep = .
	replace keep = 1 if file == "``function'_bim_short'-`function'.extract"
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
