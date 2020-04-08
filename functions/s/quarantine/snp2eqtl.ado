*! 1.0.1 Richard Anney 26jun2019

program  snp2eqtl
local    function snp2eqtl
local    `function'_version v1
syntax, index(string asis) gtex(string asis) ldref(string asis)

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
*
qui { // module 1 - prepare folders 
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 1: prepare folders"
	noi di as text"#########################################################################"
	cd ${init_root}	
	!mkdir `function'
	cd `function'
	noi di as text"# > ........ index snp "as result"`index'"
	!mkdir `index'
	cd `index'
	ref_path_to_short, ref(`ldref')
	local `function'_ldref_short ${ref_short}
	ref_path_to_short, ref(`gtex')
	local `function'_gtex_short ${ref_short}
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 2 - define files
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 2: define files"
	noi di as text"#########################################################################"
	foreach file in bim bed fam {
		noi checkfile, file(`ldref'.`file')
		}
	noi checkfile, file(`gtex'.dta)
	noi di as text"#########################################################################"
	noi di as text""
	}
qui { // module 3 - define LD partners in reference
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3: define LD partners in reference"
	noi di as text"#########################################################################"
	noi di as text"# > .................. identifying ld-partners"
	noi di as text"# > ... reference data "as result"``function'_ldref_short'"
	local already_created ${init_root}/`function'/`index'/`function'_`index'-ld.dta
	capture confirm file `already_created'
	if !_rc {
		use `function'_`index'-ld.dta, clear
		checkfile, file(`already_created')
		}
	else {
		!$plink --bfile `ldref' --ld-snp `index' --r2 --out `function'_`index'
		!$tabbed `function'_`index'.ld
		import delim using `function'_`index'.ld.tabbed, clear
		keep snp_b r2
		rename snp_b snp
		tostring r2, replace
		duplicates drop
		replace r2 = "index" if snp == "`index'"
		save `function'_`index'-ld.dta, replace
		checkfile, file(`already_created')
		}
	noi di as text"#########################################################################"
	noi di as text""
	}
qui { // module 4 - merge with gtex
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 4: merge with gtex"
	noi di as text"#########################################################################"
	noi di as text"# > .................. merge against gtex"
	noi di as text"# > ........ gtex name "as result"``function'_gtex_short'"
	local already_created ${init_root}/`function'/`index'/`function'_`index'-by-``function'_gtex_short'.dta
	capture confirm file `already_created'
	if !_rc {
		checkfile, file(`already_created')
		}
	else {
		merge 1:m snp using $snp2eqtl_gtex
		keep if _m == 3
		drop _m
		sort symbol tissue pval
		compress
		save `function'_`index'-by-``function'_gtex_short'.dta, replace
		checkfile, file(`already_created')
		}
	noi di as text"#########################################################################"
	noi di as text""
	}
qui { // module 5 - write log file
	file open myfile using "`function'_`index'-by-``function'_gtex_short'.log", write replace
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
	file write myfile `"# > ........ index SNP `index'"' _n
	file write myfile `"# > ........ reference ``function'_ldref_short'"' _n
	file write myfile `"# > ............. gtex ``function'_gtex_short'"' _n
	file write myfile `"# > ....... saved data `function'_`index'-by-``function'_gtex_short'"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > overview of `function'"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > 1. define LD partners with index SNP"' _n
	file write myfile `"# > 2. merge snp-list against gtex data "' _n
	file write myfile `"# > 3. report in table"' _n
	file write myfile `"#########################################################################"' _n
	file close myfile	
	}
qui { // module 6 - clean folder
	files2dta, dir(`c(pwd)')
	gen keep = .
	replace keep = 1 if file == "`function'_`index'-by-``function'_gtex_short'.dta"
	replace keep = 1 if file == "`function'_`index'-by-``function'_gtex_short'.log"
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

	
