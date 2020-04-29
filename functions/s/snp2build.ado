*! 1.0.3 Richard Anney 16apr2019
* clean up-front functions
*! 1.0.2 Richard Anney 29apr2020

program  snp2build
local    function snp2build
local    `function'_version v2
syntax , ref(string asis) 

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
qui { // module 2 - define files
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 1: check files"
	noi di as text"#########################################################################"
	noi checkfile, file(`ref')
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 3 - measure overlap
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 2: measure overlap"
	noi di as text"#########################################################################"
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
		gsort -percentMatched
		global `function'_build 	= build[1]
		noi di as text"# > ......... build is " as result"${`function'_build}"
	noi di as text"#########################################################################"
	noi di as text""	
	}
	}
restore
noi di as text"#########################################################################"
noi di as text"# Completed: $S_DATE $S_TIME"
noi di as text"#########################################################################"
end;	
