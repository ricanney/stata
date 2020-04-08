program  summaryqc2rg
local    function summaryqc2rg
local    `function'_version v1
syntax ,  summaryqc(string asis) w_hm3(string asis) ldsc(string asis) [munge_sumstats(string asis)]

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
	clear
	set obs 1
	gen summaryqc = "`summaryqc'"
	split summaryqc,p(" "",")
	drop summaryqc
	sxpose, clear
	gen obs = _n
	sum obs
	local `function'_summaryqc_N `r(max)'
	foreach i of num 1/ ``function'_summaryqc_N' {
		noi bim_path_to_short, bim(`=_var1[`i']')
		local `function'_summaryqc`i'_short ${bim_short}
		}
	cd ${init_root}	
	!mkdir `function'
	cd `function'
	noi di as text"# > ..... project name "as result"``function'_summaryqc1_short'"
	!mkdir ``function'_summaryqc1_short'
	cd ``function'_summaryqc1_short'
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 2 - define files
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 2: define files"
	noi di as text"#########################################################################"
	qui { // check files
		clear
		set obs 1
		gen summaryqc = "`summaryqc'"
		split summaryqc,p(" "",")
		drop summaryqc
		sxpose, clear
		foreach i of num 1/ ``function'_summaryqc_N' {
			noi checkfile, file(`=_var1[`i']'-summaryqc.dta)
			}
		noi checkfile, file(`w_hm3'w_hm3.dta)
		noi checkfile, file(`w_hm3'w_hm3.snplist)
		noi checkfile, file(`ldsc')
		noi di as text"#########################################################################"
		noi di as text""
		}
	}
qui { // module 3 - convert summaryqc to sumstats
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3: apply munge_sumstats"
	noi di as text"#########################################################################"
	foreach i of num 1/ ``function'_summaryqc_N' {
		local already_created ${init_root}/summaryqc2sumstats/``function'_summaryqc`i'_short'/``function'_summaryqc`i'_short'-summaryqc2sumstats-hw3.sumstats.gz
		capture confirm file `already_created'
		if !_rc {
			noi checkfile, file(`already_created')
			}
		else {
			clear
			set obs 1
			gen summaryqc = "`summaryqc'"
			split summaryqc,p(" "",")
			drop summaryqc
			sxpose, clear
			summaryqc2sumstats, summaryqc(`=_var1[`i']') w_hm3(`w_hm3') munge_sumstats(`munge_sumstats')
			noi checkfile, file(`already_created')
			}
		}
	noi di as text"#########################################################################"
	noi di as text""
	}
qui { // module 4 - define rg
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 4: apply `ldsc'"
	noi di as text"#########################################################################"
	noi di as text"# > ..... project name "as result"``function'_summaryqc1_short'"
	noi di as text"#########################################################################"
	noi di as text"# > ...... estimate rg"
	clear
	set obs ``function'_summaryqc_N'
	gen summaryqc = ""
	foreach i of num 1/ ``function'_summaryqc_N' {
			replace summaryqc = "``function'_summaryqc`i'_short'" in  `i'
			}
	gen sumstat = "${init_root}/summaryqc2sumstats/" + summaryqc + "/" + summaryqc + "-summaryqc2sumstats-hw3.sumstats.gz,"
	replace sumstat = "${init_root}/summaryqc2sumstats/" + summaryqc + "/" + summaryqc + "-summaryqc2sumstats-hw3.sumstats.gz" in ``function'_summaryqc_N'
	sxpose, clear
	drop in 1
	gen _var999 = ""
	gen _var = ""
	for var _var1-_var999: replace _var = _var + X
	!python `ldsc' --rg `=_var[1]' --ref-ld `w_hm3'w_hm3/w_hm3 --w-ld `w_hm3'w_hm3/w_hm3 --out ``function'_summaryqc1_short'-`function'-hw3-rg
	import delim using  ``function'_summaryqc1_short'-`function'-hw3-rg.log, varnames(nonames) clear		
	split v1,p("/home/mdnra/software/stata/data/summaryqc2sumstats/")
	keep if v13 !=""
	keep v1
	split v1,p("/")
	rename (v18 v116)(p1 p2)
	keep p1 p2 v1
	gen v2 = v1
	foreach i of num 1/10 {
			replace v2 = subinstr(v2,"  "," ",.)
			}
	split v2,p(" ")
	rename (v23-v212) (rg se z p h2_obs h2_obs_se h2_int h2_int_se gcov_int gcov_int_se)
	keep p1 p2 rg se z p h2_obs h2_obs_se h2_int h2_int_se gcov_int gcov_int_se
	set obs ``function'_summaryqc_N'
	noi di as text"# > ...... index study "as result`"``function'_summaryqc1_short'"'
	noi di as text"#########################################################################"
	foreach i of num 2 / ``function'_summaryqc_N' {
			noi di as text"# > . comparison study "as result`"`=p2[`i']'"'
			noi di as text"# > ............... rg "as result`"`=rg[`i']' (`=se[`i']')"'
			noi di as text"# > ................ z "as result`"`=z[`i']'"'
			noi di as text"# > ................ p "as result`"`=p[`i']'"'
			noi di as text"# > ... gcov intercept "as result`"`=gcov_int[`i']' (`=gcov_int_se[`i']')"'
			noi di as text"#########################################################################"
			}
	save ``function'_summaryqc1_short'-`function'.dta, replace
	}
qui { // module 4 - write log file
	file open myfile using "``function'_summaryqc1_short'-`function'.log", write replace
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
	file write myfile `"# > ...... index study ``function'_summaryqc1_short'"' _n
	foreach i of num 2 / ``function'_summaryqc_N' {
		file write myfile `"# > . comparison study `=p2[`i']'"' _n
		file write myfile `"# > ............... rg `=rg[`i']' (`=se[`i']')"' _n
		file write myfile `"# > ................ z `=z[`i']'"' _n
		file write myfile `"# > ................ p `=p[`i']'"' _n
		file write myfile `"# > ... gcov intercept `=gcov_int[`i']' (`=gcov_int_se[`i']')"' _n
		file write myfile`"#########################################################################"' _n
		}
	file write myfile `"# > ........ saved log ``function'_summaryqc1_short'-`function'-hw3-rg.log"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > overview of `function'"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > 1. wrapper for `ldsc'"' _n
	file write myfile `"#########################################################################"' _n
	file close myfile	
	}
qui { // module 5 - clean folder
	files2dta, dir(`c(pwd)')
	gen keep = .
	replace keep = 1 if file == "``function'_summaryqc1_short'-`function'-hw3-rg.log"
	replace keep = 1 if file == "``function'_summaryqc1_short'-`function'.log"
	replace keep = 1 if file == "``function'_summaryqc1_short'-`function'.dta"
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

	
