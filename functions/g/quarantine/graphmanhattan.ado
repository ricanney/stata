*! 1.0.2 Richard Anney 30apr2019
* clean up-front functions
* change output name to include chr and loc

program  graphmanhattan
local    function graphmanhattan
local    `function'_version v2
syntax , [project(string asis) maxp(real 2) minp(string asis) gws(real 7.3) str(real 6) color(string asis) summaryqc(string asis) png(string asis)]

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
	keep chr bp p
	foreach i in chr bp p {
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
	drop if p == .
	gen logp = -log10(p)
	lab var logp "observed"
	count
	global rN `r(N)'	
	noi di as text"# > .................. "as result"`r(N)' "as text"SNPs in dataset"
	sum p
	noi di as text"# > .. minP in dataset "as result"`: display %10.4e r(min)'"
	drop if chr > 23		   // drop chromosomes > X (X- XY and other)
	duplicates drop			   // drop any duplicate observations
	replace logp = round(logp,0.005)
	keep chr bp logp
	duplicates drop			   // drop any duplicate observations
	qui { // apply ceiling
		gen minp = "`minp'"
		if minp == "" {
			sum logp
			drop minp
			gen minp = `r(max)'
			replace minp = minp + 2
			replace minp = round(minp,1)
			sum minp
			local minp `r(max)'
			}
		else {
			replace logp = `minp' if logp > `minp' 
			}
		}
	qui { // apply floor
		drop if logp < `maxp'                  
		}
	qui { // convert bp to new_bp
		gen new_bp = bp 
		replace new_bp = new_bp + 270000000 if chr == 2
		replace new_bp = new_bp + 534000000 if chr == 3
		replace new_bp = new_bp + 754000000 if chr == 4
		replace new_bp = new_bp + 966000000 if chr == 5
		replace new_bp = new_bp + 1169000000 if chr == 6
		replace new_bp = new_bp + 1361000000 if chr == 7
		replace new_bp = new_bp + 1543000000 if chr == 8
		replace new_bp = new_bp + 1710000000 if chr == 9
		replace new_bp = new_bp + 1870000000 if chr == 10
		replace new_bp = new_bp + 2025000000 if chr == 11
		replace new_bp = new_bp + 2182000000 if chr == 12
		replace new_bp = new_bp + 2337000000 if chr == 13
		replace new_bp = new_bp + 2473000000 if chr == 14
		replace new_bp = new_bp + 2602000000 if chr == 15
		replace new_bp = new_bp + 2725000000 if chr == 16
		replace new_bp = new_bp + 2837000000 if chr == 17
		replace new_bp = new_bp + 2942000000 if chr == 18
		replace new_bp = new_bp + 3044000000 if chr == 19
		replace new_bp = new_bp + 3124000000 if chr == 20
		replace new_bp = new_bp + 3210000000 if chr == 21
		replace new_bp = new_bp + 3278000000 if chr == 22
		replace new_bp = new_bp + 3350000000 if chr == 23
		replace new_bp = round(new_bp/100000,1)
		drop bp
		duplicates drop			   // drop any duplicate observations
		}
	qui { // define ticks
		foreach chr of num 1/23 {
			sum new_bp if chr == `chr'
			local mtick`chr' `r(mean)'
			}
		}
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 4 - plot manhattan
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 4: plot manhattan"
	noi di as text"#########################################################################"
	qui { // define color
		gen color = "`color'"
		if color == "" {
			colorscheme 8, palette(Blues)
			}
		else {
			colorscheme 8, palette(`color')
			}	
		local color3	"mlc("`r(color7)'") mfc("`r(color7)'")"
		local color4	"mlc("`r(color8)'") mfc("`r(color8)'")"
		}
	qui { // define labels	
		sum logp
		gen tmpx = `r(max)' + 2.1
		replace tmpx = round(tmpx,2)
		sum tmpx
		local tmpmax `r(max)'
		gen tmpmin = `minp'
		local tmp_symbol "msymbol(o) msize(small)"
		gen range = tmpx - `maxp'
		replace range = range / 4
		replace range = round(range,1)
		sum range
		local tmpgap `r(max)'
		}
	qui { // plot graphs
		sum chr
		if `r(max)' == 23 {
			#delimit;
			tw scatter logp new_bp if chr == 1 ,  `tmp_symbol' `color3'
			|| scatter logp new_bp if chr == 2 ,  `tmp_symbol' `color4'
			|| scatter logp new_bp if chr == 3 ,  `tmp_symbol' `color3'
			|| scatter logp new_bp if chr == 4 ,  `tmp_symbol' `color4'
			|| scatter logp new_bp if chr == 5 ,  `tmp_symbol' `color3'
			|| scatter logp new_bp if chr == 6 ,  `tmp_symbol' `color4'
			|| scatter logp new_bp if chr == 7 ,  `tmp_symbol' `color3'
			|| scatter logp new_bp if chr == 8 ,  `tmp_symbol' `color4'
			|| scatter logp new_bp if chr == 9 ,  `tmp_symbol' `color3'
			|| scatter logp new_bp if chr == 10 , `tmp_symbol' `color4'
			|| scatter logp new_bp if chr == 11 , `tmp_symbol' `color3'
			|| scatter logp new_bp if chr == 12 , `tmp_symbol' `color4'
			|| scatter logp new_bp if chr == 13 , `tmp_symbol' `color3'
			|| scatter logp new_bp if chr == 14 , `tmp_symbol' `color4'
			|| scatter logp new_bp if chr == 15 , `tmp_symbol' `color3'
			|| scatter logp new_bp if chr == 16 , `tmp_symbol' `color4'
			|| scatter logp new_bp if chr == 17 , `tmp_symbol' `color3'
			|| scatter logp new_bp if chr == 18 , `tmp_symbol' `color4'
			|| scatter logp new_bp if chr == 19 , `tmp_symbol' `color3'
			|| scatter logp new_bp if chr == 20 , `tmp_symbol' `color4'
			|| scatter logp new_bp if chr == 21 , `tmp_symbol' `color3'
			|| scatter logp new_bp if chr == 22 , `tmp_symbol' `color4'
			|| scatter logp new_bp if chr == 23 , `tmp_symbol' `color3'
			title("`project'")
			ytitle("-log10(p)"" ")  ylabel(`maxp'(`tmpgap')`tmpmax')
			xtitle(" ""Chromosome")	xlabel(none)
			yline(`gws', lp(dash) lc("203 024 029") lw(thin)) 
			yline(`str', lp(dash) lc("065 171 093") lw(thin)) 
			xmlabel(`mtick1' "1" `mtick2' "2" `mtick3' "3" `mtick4' "4" `mtick5' "5" `mtick6' "6" `mtick7' "7" `mtick8' "8" `mtick9' "9" `mtick10' "10" `mtick11' "11" `mtick12' "12" `mtick13' "13" `mtick14' "14" `mtick15' "15" `mtick16' "16" `mtick17' "17" `mtick18' "18" `mtick19' "19" `mtick20' "20" `mtick21' "21" `mtick22' "22" `mtick23' "23" , nogrid)
			fysize(100) fxsize(500)
			legend(off)
			nodraw saving(`project'_`function'.gph, replace)
			;
			#delimit cr
			}
		else {
			#delimit;
			tw scatter logp new_bp if chr == 1 ,  `tmp_symbol' `color3'
			|| scatter logp new_bp if chr == 2 ,  `tmp_symbol' `color4'
			|| scatter logp new_bp if chr == 3 ,  `tmp_symbol' `color3'
			|| scatter logp new_bp if chr == 4 ,  `tmp_symbol' `color4'
			|| scatter logp new_bp if chr == 5 ,  `tmp_symbol' `color3'
			|| scatter logp new_bp if chr == 6 ,  `tmp_symbol' `color4'
			|| scatter logp new_bp if chr == 7 ,  `tmp_symbol' `color3'
			|| scatter logp new_bp if chr == 8 ,  `tmp_symbol' `color4'
			|| scatter logp new_bp if chr == 9 ,  `tmp_symbol' `color3'
			|| scatter logp new_bp if chr == 10 , `tmp_symbol' `color4'
			|| scatter logp new_bp if chr == 11 , `tmp_symbol' `color3'
			|| scatter logp new_bp if chr == 12 , `tmp_symbol' `color4'
			|| scatter logp new_bp if chr == 13 , `tmp_symbol' `color3'
			|| scatter logp new_bp if chr == 14 , `tmp_symbol' `color4'
			|| scatter logp new_bp if chr == 15 , `tmp_symbol' `color3'
			|| scatter logp new_bp if chr == 16 , `tmp_symbol' `color4'
			|| scatter logp new_bp if chr == 17 , `tmp_symbol' `color3'
			|| scatter logp new_bp if chr == 18 , `tmp_symbol' `color4'
			|| scatter logp new_bp if chr == 19 , `tmp_symbol' `color3'
			|| scatter logp new_bp if chr == 20 , `tmp_symbol' `color4'
			|| scatter logp new_bp if chr == 21 , `tmp_symbol' `color3'
			|| scatter logp new_bp if chr == 22 , `tmp_symbol' `color4'
			title("`project'")
			ytitle("-log10(p)"" ")  ylabel(`maxp'(`tmpgap')`tmpmax')
			xtitle(" ""Chromosome")	xlabel(none)
			yline(`gws', lp(dash) lc("203 024 029") lw(thin)) 
			yline(`str', lp(dash) lc("065 171 093") lw(thin)) 
			xmlabel(`mtick1' "1" `mtick2' "2" `mtick3' "3" `mtick4' "4" `mtick5' "5" `mtick6' "6" `mtick7' "7" `mtick8' "8" `mtick9' "9" `mtick10' "10" `mtick11' "11" `mtick12' "12" `mtick13' "13" `mtick14' "14" `mtick15' "15" `mtick16' "16" `mtick17' "17" `mtick18' "18" `mtick19' "19" `mtick20' "20" `mtick21' "21" `mtick22' "22" , nogrid)
			fysize(100) fxsize(500)
			legend(off)
			nodraw saving(`project'_`function'.gph, replace)
			;
			#delimit cr
			}

		}
	noi di as text"# > .... plot graph to "as result"`project'_`function'.gph"
	clear 
	set obs 1
	gen png = "`png'"
	if png == "yes" {
		graph use `project'_`function'.gph
		noi di as text"# > ..... saved png to "as result"`project'_`function'.png"
		graph export `project'_`function'.png, as(png) height(1000) width(3000) replace
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
	file write myfile `"# > ..... saved graphs `project'_`function'.gph"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > overview of `function'"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > 1. convert p to log p"' _n
	file write myfile `"# > 2. convert bp to continuous "' _n
	file write myfile `"# > 3. plot log p over bp"' _n
	file write myfile `"#########################################################################"' _n
	file close myfile	
	}
noi di as text"#########################################################################"
noi di as text"# Completed: $S_DATE $S_TIME"
noi di as text"#########################################################################"
qui cd `return'
restore
end;
	
 
