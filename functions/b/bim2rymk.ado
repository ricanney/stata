*! 1.0.2 Richard Anney 15apr2019

program  bim2rymk
local    function bim2rymk
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
	noi checkfile, file(`bim'.bim)
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 3 - use _bim.dta
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3: use _bim.dta"
	noi di as text"#########################################################################"	
	bim2count, bim(`bim')
	local `function'_Nin ${bim2count_snp}
	noi di as text"# > ..... SNPs in file "as result "${bim2count_snp}"
	noi di as text"# > ........... import "as result "``function'_bim_short'"
	local already_created ${init_root}/bim2dta/``function'_bim_short'/``function'_bim_short'_bim.dta
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		}
	else {
		noi di as text"# > .................. create "as result "``function'_bim_short'_bim.dta"
		noi bim2dta, bim(`bim')
		noi checkfile, file(`already_created')
		}
	noi di as text"#########################################################################"
	}	
qui { // module 4 - keeping RYMK
	noi di as text""	
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3: keeping RYMK"
	noi di as text"#########################################################################"
	use ${init_root}/bim2dta/``function'_bim_short'/``function'_bim_short'_bim.dta, clear
	keep snp a1 a2 gt
	gen keep = .
	foreach gt in R Y M K {
		replace keep = 1 if gt ==  "`gt'"
		}
	keep if keep == 1
	drop keep
	foreach allele in a1 a2 {
		gen drop = 1
		foreach a in A C G T {
			replace drop = . if `allele' ==  "`a'"
			}
		drop if drop == 1
		drop drop
		}
	order snp a1 a2 gt
	keep  snp a1 a2 gt
	compress
	count
	local `function'_Nout `r(N)'
	noi di as text"# > ..... SNPs in file "as result "``function'_Nout'"
	noi di as text"# > .......... save to "as result"``function'_bim_short'_rymk.dta"
	save ``function'_bim_short'_rymk.dta,replace
	noi di as text"#########################################################################"
	noi di as text""
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
	file write myfile `"# > .... plink version ${plink}"' _n
	file write myfile `"# > . SNPs in original ``function'_Nin'"' _n
	file write myfile `"# >  SNPs in processed ``function'_Nout'"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > overview of `function'"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > 1. use / create _bim.dta "' _n
	file write myfile `"# > 2. remove on R Y M K genotypes"' _n
	file write myfile `"#########################################################################"' _n
	file close myfile	
	}
qui { // module 6 - clean folder
	files2dta, dir(`c(pwd)')
	gen keep = .
	replace keep = 1 if file == "``function'_bim_short'_rymk.dta"
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
