program  snp2refid
local    function snp2refid
local    `function'_version v2
syntax , ref(string asis) project(string asis)

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

local return `c(pwd)'
qui { // module 1 - prepare folders 
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 1: prepare folders"
	noi di as text"#########################################################################"
	ref_path_to_short, ref(`ref')
	local `function'_ref_short ${ref_short}
	cd ${init_root}	
	!mkdir `function'
	cd `function'
	noi di as text"# > ..... project name "as result"`project'-``function'_ref_short'"
	!mkdir `project'-``function'_ref_short'
	cd `project'-``function'_ref_short'
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 2 - create _bim.dta
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 2: create _bim.dta"
	noi di as text"#########################################################################"
	local already_created ${init_root}/bim2dta/``function'_ref_short'/``function'_ref_short'_bim.dta
	capture confirm file `already_created'
	if !_rc {
				noi checkfile, file(`already_created')
				}
	else {	
		bim2dta, bim(`ref')
		noi checkfile, file(`already_created')
		}
	noi di as text"#########################################################################"
	noi di as text""
	}
qui { // module 3 - merge data
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3: merge files"
	noi di as text"#########################################################################"
	count
	local `function'_Nin `r(N)' 
	noi di as text"# > ..... SNPs in file "as result "``function'_Nin'"
	noi di as text"# > ... updating using "as result"``function'_ref_short'"
	compress
	for var chr bp: tostring X, replace force
	replace chr = subinstr(chr,"chr","",.)
	replace chr = strupper(chr)	
	replace chr = "23" if chr == "X"
	replace chr = "23" if chr == "X_nonPAR"
	replace chr = "24" if chr == "Y"
	replace chr = "25" if chr == "XY"
	replace chr = "26" if chr == "MT"
	noi di as text"# > .................. generating loc_name"
	checkloc_name
	gen snp999 = ""
	aorder
	drop snp-snp999 chr bp
	noi di as text"# > .................. merge with reference"
	merge m:1 loc_name using `already_created'
	keep if _m == 3
	compress
	qui { // remove duplicates
		noi di as text"# > .................. remove duplicates"
		duplicates tag snp, gen(tag)
		egen keep = seq(),by(snp tag)
		keep if keep == 1
		drop _m gt loc_name  tag keep 
		}
	order chr bp snp
	for var chr bp: destring X, replace
	sort chr bp
	count
	local `function'_Nout `r(N)' 
	noi di as text"# > ..... SNPs in file "as result "``function'_Nout'"
	save `project'-``function'_ref_short'-`function'.dta, replace
	noi di as text"#########################################################################"
	noi di as text""
	}
qui { // module 4 - write log file
	file open myfile using "`project'-``function'_ref_short'-`function'.log", write replace
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
	file write myfile `"# > ..... project name `project'"' _n
	file write myfile `"# > .... SNPs in input ``function'_Nin'"' _n
	file write myfile `"# > ........ reference ``function'_ref_short'"' _n
	file write myfile `"# > .... SNPs in input ``function'_Nout'"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > overview of `function'"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > 1. rename chromosomes"' _n
	file write myfile `"# > 2. merge with ref"' _n
	file write myfile `"# > 3. remove duplicates"' _n
	file write myfile `"#########################################################################"' _n
	file close myfile	
	}	
qui { // module 5 - clean folder
	files2dta, dir(`c(pwd)')
	gen keep = .
	replace keep = 1 if file == "`project'-``function'_ref_short'-`function'.log"
	replace keep = 1 if file == "`project'-``function'_ref_short'-`function'.dta"
	drop if keep == 1
	gen script = "erase " + file
	outsheet script using temp.do, non noq replace 
	do temp.do
	erase temp.do
	clear
	use 		`project'-``function'_ref_short'-`function'.dta, clear
	erase 	`project'-``function'_ref_short'-`function'.dta
	}
qui cd `return'
noi di as text"#########################################################################"
noi di as text"# Completed: $S_DATE $S_TIME"
noi di as text"#########################################################################"
end; 
