*! 1.0.2 Richard Anney 30apr2019
* clean up-front functions
* change output name to include chr and loc

program  graphqq
local    function graphqq
local    `function'_version v2
syntax , [project(string asis) maxp(real 0) minp(real 20) gws(real 7.3) str(real 6) color(string asis) summaryqc(string asis) png(string asis)]

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
	cd ${init_root}	
	!mkdir `function'
	cd `function'
	gen `function'_project = "`project'"
	if `function'_project == "" {
		gen `function'_summaryqc = "`summaryqc'"
		if `function'_summaryqc == "" {
			noi di as text"# > ..... project name "as result"unknown"
			local project unknown
			!mkdir unknown
			cd unknown
			drop `function'_project `function'_summaryqc
			}
		else {
			ref_path_to_short, ref(`summaryqc')
			local project ${ref_short}
			noi di as text"# > ..... project name "as result"`project'"
			!mkdir `project'
			cd `project'
			use `summaryqc'-summaryqc.dta, clear
			}
		}
	else {
		noi di as text"# > ..... project name "as result"`project'"
		!mkdir `project'
		cd `project'
		drop `function'_project
		}
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 2 - define variables
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 2: define variables"
	noi di as text"#########################################################################"
	keep p
	foreach i in p {
		capture confirm variable `i'
		if !_rc {
			noi di as text"# > .................. "as result"`i' "as text"is present"
			}
		else {
			noi di as text"# > .................. "as error"`i' "as text"is absent"
			exit
			}
		capture confirm variable numeric `i'
		if !_rc {
			noi di as text"# > .................. "as result"`i' "as text"is numeric"
			}
		else {
			noi di as text"# > .......... convert "as result"`i' "as text"to numeric"
			destring `i', replace
			}
		}
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 3 - process variables
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3: process variables"
	noi di as text"#########################################################################"
	qui { // count observations
		drop if p == .
		count
		local `function'_N `r(N)'	
		noi di as text"# > .. SNPs in dataset "as result"``function'_N' "
		sum p
		local `function'_minP `r(min)'	
		noi di as text"# > .. minP in dataset "as result"`: display %10.4e ``function'_minP''"
		}
	qui { // define observed
		noi di as text"# > .................. define observed P"
		gen `function'_observed = -log10(p)
		lab var `function'_observed "observed"
		replace `function'_observed = `minp' if `function'_observed > `minp'	
		}
	qui { // define expected
		noi di as text"# > .................. calculate expected P"
		sort p
		gen `function'_n = _n
		gen `function'_expected = -log10(`function'_n/``function'_N')
		lab var `function'_observed "expected"
		}
	qui { // calculate lambda
		noi di as text"# > .................. calculate lambda (genomic inflation factor)"
		qui { // convert p to chi
			gen `function'_chi2_observed = invchi2tail(1,p)
			}
		sum `function'_chi2_observed, detail
		gen lambda = round(`r(p50)' / 0.456,.001)
		sum lambda
		local `function'_lambda `r(max)'
		noi di as text"# > . lambda inflation "as result"``function'_lambda'"
		}
	qui { // define binomial boundaries
		qui { // prune data bins
			keep `function'_observed `function'_expected `function'_n
			for var `function'_observed `function'_expected: replace X = round(X,.005)
			egen dup = seq(),by(`function'_observed `function'_expected)
			keep if dup == 1
			keep `function'_observed `function'_expected `function'_n
			}
		qui { // calculate binomial boundaries
			noi di as text"# > .................. calculate binomial boundaries for expected P"
			expand 3
			egen obs = seq(),by(`function'_observed `function'_expected)
			sort `function'_n obs
			tostring `function'_n, replace
			gen script = ""
			replace script = "qui cii ``function'_N' " + `function'_n if obs == 1
			replace script = `"qui replace `function'_ub = r(ub) if `function'_n == ""' + `function'_n + `"""' if obs == 2
			replace script = `"qui replace `function'_lb = r(lb) if `function'_n == ""' + `function'_n + `"""' if obs == 3
			outsheet script using temp.do, non noq replace
			keep `function'_observed `function'_expected `function'_n
			duplicates drop
			gen `function'_ub = .
			gen `function'_lb = .
			do temp.do
			erase temp.do
			for var `function'_ub `function'_lb : replace X = -log10(X)
			keep    `function'_ub `function'_lb `function'_observed `function'_expected 
			duplicates drop
			sort `function'_expected
			for var  `function'_observed : drop if X < `maxp'
			}
		}
	noi di as text"#########################################################################"
	noi di as text""
	}
qui { // module 4 - plot graph
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 4: plot graph"
	noi di as text"#########################################################################"
	gen color = "`color'"
	if color == "" {
		colorscheme 8, palette(Blues)
		}
	else {
		colorscheme 8, palette(`color')
		}	
	local color1	"mlc("`r(color6)'") mfc("`r(color6)'")"
	local color2	"mlc("`r(color7)'") mfc("`r(color7)'")"
	local color3	"mlc("`r(color8)'") mfc("`r(color8)'")"
	local tmp_symbol "msymbol(o) msize(small)"
	#delimit;
	tw line `function'_expected `function'_expected , lwidth(vthin) lcolor(black)
	|| line `function'_ub       `function'_expected , lpattern(dash) lwidth(vthin) lcolor(black)
	|| line `function'_lb       `function'_expected        , lpattern(dash) lwidth(vthin) lcolor(black)
	|| scatter `function'_observed `function'_expected if (`function'_observed <  `str')                               , `tmp_symbol' `color1'
	|| scatter `function'_observed `function'_expected if (`function'_observed >= `str' & `function'_observed  < `gws'), `tmp_symbol' `color2'
	|| scatter `function'_observed `function'_expected if (`function'_observed >= `gws')                               , `tmp_symbol' `color3'
	legend(off) 
	title("`project'") caption("lambda: ``function'_lambda'", size(small))
	xtitle(" " "Expected (-log10(P))") 
	ytitle("Observed -log10(p)"" ")
	xlabel(#5) 	ylabel(#5)
	fysize(100) fxsize(100)
	nodraw
	saving(`project'_`function'.gph, replace)
	;
	#delimit cr
	noi di as text"# > .... plot graph to "as result"`project'_`function'.gph"
	clear 
	set obs 1
	gen png = "`png'"
	if png == "yes" {
		graph use `project'_`function'.gph
		noi di as text"# > ..... saved png to "as result"`project'_`function'.png"
		graph export `project'_`function'.png, as(png) height(2000) width(3000) replace
		window manage close graph
		}
	else {
		}
	noi di as text"#########################################################################"
	noi di as text""
	}
qui { // module 5 - write log file
	file open myfile using "`project'_`function'.log", write replace
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
	file write myfile `"# > .. SNPs in dataset ``function'_N'"' _n
	file write myfile `"# > .. minP in dataset `: display %10.4e ``function'_minP''"' _n
	file write myfile `"# > ........... lambda ``function'_lambda'"' _n
	file write myfile `"# > ..... saved graphs `project'_`function'.gph"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > overview of `function'"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > 1. convert p to log p"' _n
	file write myfile `"# > 2. create expected p "' _n
	file write myfile `"# > 3. define lambda"' _n
	file write myfile `"# > 4. define binomial boundaries"' _n
	file write myfile `"#########################################################################"' _n
	file close myfile	
	}
noi di as text"#########################################################################"
noi di as text"# Completed: $S_DATE $S_TIME"
noi di as text"#########################################################################"
qui cd `return'
restore
end;
		
