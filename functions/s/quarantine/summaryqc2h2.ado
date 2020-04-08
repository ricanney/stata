program  summaryqc2h2
local    function summaryqc2h2
local    `function'_version v1
syntax ,  summaryqc(string asis) w_hm3(string asis) ldsc(string asis) [samp_prev(string asis) pop_prev(string asis) munge_sumstats(string asis)]

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
	noi bim_path_to_short, bim(`summaryqc')
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
		noi checkfile, file(`ldsc')
		noi di as text"#########################################################################"
		noi di as text""
		}
	}
qui { // module 3 - convert summaryqc to sumstats
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3: apply munge_sumstats"
	noi di as text"#########################################################################"
	local already_created ${init_root}/summaryqc2sumstats/``function'_bim_short'/``function'_bim_short'-summaryqc2sumstats-hw3.sumstats.gz
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		}
	else {
		summaryqc2sumstats, summaryqc(`summaryqc') w_hm3(`w_hm3') munge_sumstats(`munge_sumstats')
		noi checkfile, file(`already_created')
		}
	noi di as text"#########################################################################"
	noi di as text""
	}
qui { // module 3 - define h2
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3: apply `ldsc'"
	noi di as text"#########################################################################"
	local `function'_sumstat ${init_root}/summaryqc2sumstats/``function'_bim_short'/``function'_bim_short'-summaryqc2sumstats-hw3.sumstats.gz
	noi di as text"# > ..... project name "as result"``function'_bim_short'"
	noi di as text"#########################################################################"
	qui { // observed scale
		noi di as text"# > ...... estimate h2 (observed scale)"
		!python `ldsc' --h2 ``function'_sumstat' --ref-ld `w_hm3'w_hm3/w_hm3 --w-ld `w_hm3'w_hm3/w_hm3 --out ``function'_bim_short'-`function'-hw3-observed-h2
		import delim using  ``function'_bim_short'-`function'-hw3-observed-h2.log, varnames(nonames) clear
		split v1,p("Total Observed scale h2: ")
		split v12,p("("")")
		destring v121, replace force
		sum v121
		local `function'_o_h2 `r(max)'
		destring v122, replace force
		sum v122
		local `function'_o_h2_se `r(max)'
		keep v1
		split v1,p("Lambda GC: ")
		destring v12, replace force
		sum v12
		local `function'_o_lambda `r(max)'
		keep v1
		split v1,p("Mean Chi^2: ")
		destring v12, replace force
		sum v12
		local `function'_o_chi `r(max)'
		keep v1	
		split v1,p("Intercept: ")
		split v12,p("("")")
		destring v121, replace force
		sum v121
		local `function'_o_intercept `r(max)'
		destring v122, replace force
		sum v122
		local `function'_o_intercept_se `r(max)'
		keep v1
		noi di as text"# >  h2 observed-scale "as result`"``function'_o_h2' (``function'_o_h2_se')"'
		noi di as text"# > ........ lambda GC "as result`"``function'_o_lambda'"'
		noi di as text"# > ....... mean chi^2 "as result`"``function'_o_chi'"'
		noi di as text"# > ........ intercept "as result`"``function'_o_intercept' (``function'_o_intercept_se')"'
		noi di as text"#########################################################################"
		}
	qui { // liability scale
		clear
		set obs 1
		gen pop = "`pop_prev'"
		if pop == "" {
			noi di as text"# > ...... estimate h2 (liability scale)"
			noi di as text"# > .................. pop_prev / samp_prev not specified"
			noi di as text"#########################################################################"
			noi di as text""
			}
		else {
			noi di as text"# > ...... estimate h2 (liability scale)"
			!python `ldsc' --h2 `already_created' --ref-ld `w_hm3'w_hm3/w_hm3 --w-ld `w_hm3'w_hm3/w_hm3 --pop-prev `pop_prev' --samp-prev `samp_prev' --out ``function'_bim_short'-`function'-hw3-liability-h2
			import delim using  ``function'_bim_short'-`function'-hw3-liability-h2.log, varnames(nonames) clear
			split v1,p("Total Liability scale h2: ")
			split v12,p("("")")
			destring v121, replace force
			sum v121
			local `function'_l_h2 `r(max)'
			destring v122, replace force
			sum v122
			local `function'_l_h2_se `r(max)'
			keep v1
			split v1,p("Lambda GC: ")
			destring v12, replace force
			sum v12
			local `function'_l_lambda `r(max)'
			keep v1
			split v1,p("Mean Chi^2: ")
			destring v12, replace force
			sum v12
			local `function'_l_chi `r(max)'
			keep v1	
			split v1,p("Intercept: ")
			split v12,p("("")")
			destring v121, replace force
			sum v121
			local `function'_l_intercept `r(max)'
			destring v122, replace force
			sum v122
			local `function'_l_intercept_se `r(max)'
			keep v1
			noi di as text"# > h2 liability-scale "as result`"``function'_l_h2' (``function'_l_h2_se')"'
			noi di as text"# >  based on pop-prev "as result`"`pop_prev'"'
			noi di as text"# > based on samp-prev "as result`"`samp_prev'"'
			noi di as text"# > ........ lambda GC "as result`"``function'_l_lambda'"'
			noi di as text"# > ....... mean chi^2 "as result`"``function'_l_chi'"'
			noi di as text"# > ........ intercept "as result`"``function'_l_intercept' (``function'_l_intercept_se')"'
			noi di as text"#########################################################################"
			noi di as text""
			}
		}
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
	file write myfile `"# >  h2 observed-scale ``function'_o_h2' (``function'_o_h2_se')"' _n
	clear
	set obs 1
	gen pop = "`pop_prev'"
	if pop == "" {
		file write myfile `"# > h2 liability-scale <prevalence not specified>"' _n
		file write myfile `"# > ......... pop_prev NaN"' _n
		file write myfile `"# > ........ samp_prev NaN"' _n
		}
	else {	
		file write myfile `"# > h2 liability-scale ``function'_l_h2' (``function'_l_h2_se')"' _n
		file write myfile `"# > ......... pop_prev `pop_prev'"' _n
		file write myfile `"# > ........ samp_prev `samp_prev'"' _n
		file write myfile `"# > ........ lambda GC ``function'_o_lambda'"' _n
		file write myfile `"# > ....... mean chi^2 ``function'_o_chi'"' _n
		file write myfile `"# > ........ intercept ``function'_o_intercept' (``function'_o_intercept_se')"'	_n
		}
	
	file write myfile `"# > ........ saved log ``function'_bim_short'-`function'-hw3-observed-h2.log"' _n
	file write myfile `"# > ........ saved log ``function'_bim_short'-`function'-hw3-liability-h2.log"' _n
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
	replace keep = 1 if file == "``function'_bim_short'-`function'-hw3-observed-h2.log"
	replace keep = 1 if file == "``function'_bim_short'-`function'-hw3-liability-h2.log"
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

	
