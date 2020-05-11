program  summaryqc2sumstats
local    function summaryqc2sumstats
local    `function'_version v2
syntax ,  summaryqc(string asis) w_hm3(string asis) munge_sumstats(string asis) 

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
	bim_path_to_short, bim(`summaryqc')
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
	qui { // check files
		noi checkfile, file(`summaryqc'-summaryqc.dta)
		noi checkfile, file(`w_hm3'w_hm3.dta)
		noi checkfile, file(`w_hm3'w_hm3.snplist)
		noi checkfile, file(`munge_sumstats')
		noi di as text"#########################################################################"
		noi di as text""
		}
	}
qui { // module 3 - convert summaryqc to sumstats
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3: apply `munge_sumstats'"
	noi di as text"#########################################################################"
	use `summaryqc'-summaryqc.dta, clear
	merge 1:1 snp using `w_hm3'w_hm3.dta
	keep if _m == 3
	sort snp
	keep snp a1 a2 z p n
	renvars, upper
	missings dropobs, force
	for var  Z P N : tostring X, replace
	outsheet SNP A1 A2 Z P N     using ``function'_bim_short'.tmp, noq replace
	!python "`munge_sumstats'" --sumstats ``function'_bim_short'.tmp --out ``function'_bim_short'-`function'-hw3 --merge-alleles `w_hm3'w_hm3.snplist
	noi checkfile, file(``function'_bim_short'-`function'-hw3.sumstats)
	noi checkfile, file(``function'_bim_short'-`function'-hw3.sumstats.gz)
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
	file write myfile `"# > ....... saved data ``function'_bim_short'-`function'-hw3.sumstats.gz"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > overview of `function'"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > 1. wrapper for `munge_sumstats'"' _n
	file write myfile `"#########################################################################"' _n
	file close myfile	
	}
qui { // module 5 - clean folder
	files2dta, dir(`c(pwd)')
	gen keep = .
	replace keep = 1 if file == "``function'_bim_short'-`function'-hw3.log"
	replace keep = 1 if file == "``function'_bim_short'-`function'-hw3.sumstats.gz"
	replace keep = 1 if file == "``function'_bim_short'-`function'-hw3.sumstats"
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

	
