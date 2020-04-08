*! 1.0.2 Richard Anney 23apr2019
* clean up-front functions
*! 1.0.1 Richard Anney 21dec2018

program  fam2dta
local    function fam2dta
local    `function'_version v2
syntax , fam(string asis)  

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
	bim_path_to_short, bim(`fam')
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
	noi checkfile, file(`fam'.fam)
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 3 - importing bim file (checking for delimit format)
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3: update variable names"
	noi di as text"#########################################################################"
	bim2count, bim(`fam')
	noi di as text"# > ... indiv. in file "as result "${bim2count_ind}"
	noi di as text"# > ........... import "as result "``function'_bim_short'"
	import delim  using `fam'.fam, clear 
	capture confirm variable v2
		if !_rc { 
			}
		else {
			rename v1 v
			split v,p(" ")
			drop v
			}
	rename (v1-v6) (fid iid fatid motid sex pheno)
	for var fid iid fatid motid: tostring X, replace
	qui { // standardise sex coding to 1 2 -9 
		destring sex, replace
		gen x = .
		replace x = 1 if sex == 1
		replace x = 1 if sex == 2
		replace sex = -9 if x == .
		drop x
		}
	noi di as text"#########################################################################"
	}
qui { // module 4 - saving file
	noi di as text""	
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 4: saving file"
	noi di as text"#########################################################################"	
	order fid iid fatid motid sex pheno
	keep  fid iid fatid motid sex pheno
	compress
	noi di as text"# > ... saving file as "as result"``function'_bim_short'_fam.dta"
	lab var fid 			"Family Identifier"
	lab var iid				"Individual Identifier"
	lab var fatid			"Father's (within family) Identifier"
	lab var motid			"Mother's (within family) Identifier"
	lab var sex				"Sex"
	lab var pheno			"Phenotype"
	lab define sex   1"Male" 2"Female" -9"Unknown"
	lab val sex sex
	save ``function'_bim_short'_fam.dta, replace
	}
qui { // module 5 - write log file
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
	file write myfile `"# > ... indiv. in file ${bim2count_ind}"' _n
	file write myfile `"# > ....... saved data ``function'_bim_short'_fam.dta"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > overview of `function'"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > 1. import fam file"' _n
	file write myfile `"# > 2. rename variable"' _n
	file write myfile `"#########################################################################"' _n
	file close myfile	
	}
qui { // module 6 - clean folder
	files2dta, dir(`c(pwd)')
	gen keep = .
	replace keep = 1 if file == "``function'_bim_short'_fam.dta"
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

	
