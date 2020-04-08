*! 1.0.3 Richard Anney 16apr2019
* clean up-front functions
*! 1.0.2 Richard Anney 21jan2019

program  bim2build
local    function bim2build
local    `function'_version v3
syntax , bim(string asis) ref(string asis) 

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
	noi di as text"# SECTION - 2: check files"
	noi di as text"#########################################################################"
	noi checkfile, file(`bim'.bim)
	noi checkfile, file(`ref')
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 3 - measure overlap
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3: measure overlap"
	noi di as text"#########################################################################"
	noi di as text"# > ........... import "as result"``function'_bim_short'.bim"
	local already_created ${init_root}/bim2dta/``function'_bim_short'/``function'_bim_short'_bim.dta
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		}
	else {
		bim2dta, bim(`bim')
		noi checkfile, file(`already_created')
		}
	use `already_created', clear
	keep snp chr bp
	duplicates drop snp, force
	noi di as text"# > ....... merge with "as result"`ref'"
	merge 1:1 snp using `ref'
	keep if _m == 3
	tostring bp, replace
	compress
	qui { // measure overlap with reference
		foreach i in 17 18 19 { 
			gen hg`i'_0 = .
			gen hg`i'_1 = .
			replace hg`i'_0 = 1 if bp == hg`i'_chromStart 
			replace hg`i'_1 = 1 if bp == hg`i'_chromEnd 
			}		
		count
		gen all = `r(N)'
		foreach i in 17 18 19 { 
			foreach j in 0 1 { 
				sum hg`i'_`j'
				gen phg`i'_`j' = r(N) / all
				}
			}
		keep in 1
		keep phg17_0 - phg19_1
		xpose, clear v
		rename v1 percentMatched
		rename _v build
		replace build = "hg17 +0" if build == "phg17_0"
		replace build = "hg17 +1" if build == "phg17_1"
		replace build = "hg18 +0" if build == "phg18_0"
		replace build = "hg18 +1" if build == "phg18_1"
		replace build = "hg19 +0" if build == "phg19_0"
		replace build = "hg19 +1" if build == "phg19_1"
		gsort -p
		gen MostLikely = "+++" in 1
		replace MostLikely = "++" if p > 0.9 & MostLikely == ""
		replace MostLikely = "+" if p > 0.8 & MostLikely == ""
		keep in 1 
		tostring per, replace force
		gen a = ""
		replace a = "global bim2build " + build
		outsheet a using bim2build_tmp.do, non noq replace
		do bim2build_tmp.do
		noi di as text"# > ......... build is " as result"${bim2build}"as text" for "as result"``function'_bim_short'.bim"
		erase bim2build_tmp.do
		}
	}
qui { // module 4 - write log file
	file open myfile using "``function'_bim_short'.build", write replace
	file write myfile `"build"' _n
	file write myfile `"${bim2build}"' _n
	file close myfile		

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
	file write myfile `"# >  most likely build ${bim2build}"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > overview of `function'"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > 1. import markers"' _n
	file write myfile `"# > 2. merge snp against known reference"' _n
	file write myfile `"# > 3. compare chr/bp versus reference"' _n
	file write myfile `"# > 4. define build as that with highest overlap"' _n
	file write myfile `"#########################################################################"' _n
	file close myfile	
	}
qui { // module 5 - clean folder
	files2dta, dir(`c(pwd)')
	gen keep = .
	replace keep = 1 if file == "``function'_bim_short'.build"
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
