*! 1.0.1 Richard Anney 17apr2019

program  bim2hwe
local    function bim2hwe
local    `function'_version v1
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
	foreach file in bim bed fam {
		noi checkfile, file(`bim'.`file')
		}
	noi checkfile, file(${plink})
	checktabbed
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 3 - calculate hwe
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3: calculate marker hardy-weinberg statistics"
	noi di as text"#########################################################################"
	noi di as text"# > .................. calculating hwe from "as result"``function'_bim_short'"
	noi di as text"# > .................. split by chromosome "
	noi di as text"# > .................. define chromosomes in model "
	import delim `bim'.bim, colr(1:1) clear
	count
	local `function'_bimN `r(N)'	
	duplicates drop 
	rename v1 chr
	save ``function'_bim_short'-chr.dta, replace
	noi di as text"# > .................. measure hwe by chromosome"
	use ``function'_bim_short'-chr.dta, clear
	expand 12
	tostring chr, replace
	egen obs = seq(),by(chr)
	sort chr obs
	gen a = ""
	replace a = "!${plink} --bfile `bim' --chr " + chr + " --hardy --out ``function'_bim_short'-chr" + chr if obs == 1
	replace a = "!${tabbed}         ``function'_bim_short'-chr" + chr + ".hwe"                             if obs == 2
	replace a = "import delim using ``function'_bim_short'-chr" + chr + ".hwe.tabbed, clear case(lower)"   if obs == 3
	replace a = "for var p : destring X, replace force"                                                    if obs == 4
	replace a = `"replace test =  "ALL" if test == "ALL(NP)""'                                             if obs == 5
	replace a = `"keep if test == "ALL""'                                                                  if obs == 6
	replace a = `"for var p : lab var X "HWE (p)""'                                                        if obs == 7
	replace a = `"dropmiss, force"'                                                                        if obs == 8
	replace a = "save ``function'_bim_short'-chr" + chr + ".hwe.dta, replace"                              if obs == 9
	outsheet a using temp.do, non noq replace
	do temp.do
	erase temp.do
	noi di as text"# > .................. join files"
	use ``function'_bim_short'-chr.dta, clear	
	tostring chr, replace
	sort chr
	gen a = ""	
	replace a = "use ``function'_bim_short'-chr" + chr +".hwe.dta, clear" if _n == 1
	replace a = "append using ``function'_bim_short'-chr" + chr +".hwe.dta" if _n != 1
	outsheet a using temp.do, non noq replace
	do temp.do
	erase temp.do
	noi di as text"# > .................. saving file as "as result"``function'_bim_short'_hwe.dta"
	count
	local `function'_hweN `r(N)'
	save ``function'_bim_short'_hwe.dta, replace
	noi di as text"#########################################################################"
	noi di as text""	
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
	file write myfile `"# > ...... SNPs in bim ``function'_bimN'"' _n
	file write myfile `"# > ...... SNPs in hwe ``function'_hweN'"' _n
	file write myfile `"# > ....... saved data ``function'_bim_short'_hwe.dta"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > overview of `function'"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > 1. calculate hwe in SNPs"' _n
	file write myfile `"#########################################################################"' _n
	file close myfile	
	}
qui { // module 5 - clean folder
	files2dta, dir(`c(pwd)')
	gen keep = .
	replace keep = 1 if file == "``function'_bim_short'_hwe.dta"
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
