*! 1.0.1 Richard Anney 17apr2019

program  bim2miss
local    function bim2miss
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
qui { // module 3 - calculate missingness	
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3: calculate missingness"
	noi di as text"#########################################################################"
	noi di as text"# > .................. calculating missingness from "as result"``function'_bim_short'"
	!$plink  --bfile `bim' --missing          --out ``function'_bim_short'
	qui { // process missingness by individual
		!$tabbed ``function'_bim_short'.imiss
		import delim using ``function'_bim_short'.imiss.tabbed, clear case(lower)
		for var fid iid: tostring X, replace force
		for var f_miss : destring X, replace force
		for var f_miss : lab var X "Frequency of missing genotypes per individual"
		dropmiss, force
		save ``function'_bim_short'_imiss.dta, replace
		noi di as text"# > ......... saved to "as result"``function'_bim_short'_imiss.dta"
		count
		local `function'_Nind `r(N)'
		noi di as text"# > .... indiv. in fam "as result ``function'_Nind'
		}
	qui { // process missingness by locus
		!$tabbed ``function'_bim_short'.lmiss
		import delim using ``function'_bim_short'.lmiss.tabbed, clear case(lower)
		for var f_miss : destring X, replace force
		for var f_miss : lab var X "Frequency of missing genotypes per SNP"
		dropmiss, force
		save ``function'_bim_short'_lmiss.dta, replace
		noi di as text"# > ......... saved to "as result"``function'_bim_short'_lmiss.dta"
		count
		local `function'_Nsnp `r(N)'
		noi di as text"# > ...... SNPs in bim "as result ``function'_Nsnp'
		replace f_miss = 0.1 if f_miss >0.1 & f_miss !=.	
		}	
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
	file write myfile `"# > ...... SNPs in bim ``function'_Nsnp'"' _n
	file write myfile `"# > .... indiv. in fam ``function'_Nind'"' _n
	file write myfile `"# > ....... saved data ``function'_bim_short'_imiss.dta"' _n
	file write myfile `"# > ....... saved data ``function'_bim_short'_lmiss.dta"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > overview of `function'"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > 1. calculate missingness using PLINK"' _n
	file write myfile `"#########################################################################"' _n
	file close myfile	
	}
qui { // module 5 - clean folder
	files2dta, dir(`c(pwd)')
	gen keep = .
	replace keep = 1 if file == "``function'_bim_short'_imiss.dta"
	replace keep = 1 if file == "``function'_bim_short'_lmiss.dta"
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
