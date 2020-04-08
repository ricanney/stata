*! 1.0.1 Richard Anney 07may2019

program  snp2phewas
local    function snp2phewas
local    `function'_version v1
syntax, index(string asis) summaryqc_dir(string asis) summaryqc_label(string asis) [report_p(real .001) summaryqc_limit(string asis) report_full(real 15) png(string asis)]

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
*
qui { // module 1 - prepare folders 
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 1: prepare folders"
	noi di as text"#########################################################################"
	cd ${init_root}	
	!mkdir `function'
	cd `function'
	noi di as text"# > ..... project name "as result"`index'"
	!mkdir `index'
	cd `index'
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 2 - define files
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 2: select gwas from databank"
	noi di as text"#########################################################################"
	noi di as text"# > ......... defining gwas to process"
	local already_created `index'-final-`function'.dta
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		use `already_created', clear
		}
	else {
		clear
		set obs 1
		gen gwas = "dummy"
		save `index'-pre-`function'.dta, replace
		gen summaryqc_limit = "`summaryqc_limit'"
		if summaryqc_limit == "" {
			dir2dta, dir(`summaryqc_dir')
			rename folder gwas
			}
		else {
			use `summaryqc_limit', clear
			}
		sort gwas
		gen of = _n
		expand 10
		egen obs = seq(),by(gwas) 
		sort gwas obs
		count if obs == 1
		gen total = `r(N)'
		for var of total: tostring X, replace
		gen script = ""
		replace script = "qui { " if obs == 1
		replace script = `"noi di as text"# > ....... processing "as result""' + gwas + " (" + of + " of " + total + `")""' if obs == 2
		replace script = `"use if snp == "`index'" using `summaryqc_dir'/"' + gwas + "/" + gwas + "-summaryqc.dta, clear" if obs == 3
		replace script = "set obs 1" if obs == 4
		replace script = `"gen gwas = ""' + gwas + `"""' if obs == 5
		replace script = "append using `index'-pre-`function'.dta" if obs == 6
		replace script = "save `index'-pre-`function'.dta, replace" if obs == 7
		replace script = "}" if obs == 8
		drop if script == ""
		egen x = seq(),by(script)
		drop if script == "qui { " & x != 1
		sum if script == "}" 
		drop if script == "}" & x != `r(max)'
		outsheet script using temp.do, non noq replace
		noi do temp.do
		erase temp.do
		use `index'-pre-`function'.dta, clear
		drop if gwas == "dummy"
		drop if snp  == ""
		qui { // map to common allele
			egen obs = seq(),by(snp)
			gen allele1 = a1 if obs == 1
			gen allele2 = a2 if obs == 1
			replace allele1 = allele1[_n-1] if allele1 == ""
			replace allele2 = allele2[_n-1] if allele2 == ""
			for var beta z : replace X = -X if a1 != allele1
			for var or l95 u95: replace X = 1/X if a1 != allele1
			drop a1 a2
			rename (allele1 allele2) (a1 a2)
			drop obs
			}
		qui { // define upper / lower boundaries	
			gen l = beta -se
			gen u = beta +se
			gsort -beta
			egen obs = seq()
			}
		qui { // define label
			gen label = gwas
			replace label = strupper(label)
			gen str6 beta_label = string(beta,"%05.4f")						
			gen str6 lower_label = string(l,"%05.4f")					
			gen str6 upper_label = string(u,"%05.4f")			
			gen str8 p_label = string(p,"%06.5f")
			replace p_label = "= " + p_label
			replace p_label = "< 0.00001" if p < 1e-5
			gen lnor = "beta = " + beta_label + " [" + lower_label + "," + upper_label + "]; p " + p_label
			order label lnor
			}
		save `index'-final-`function'.dta, replace
		noi checkfile, file(`already_created')
		}
	noi di as text"#########################################################################"
	noi di as text""
	} 
qui { // module 3 - plot pheWAS
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3: plot phewas"
	noi di as text"#########################################################################"
	use `index'-final-`function'.dta, clear
	keep l u  beta lnor label p n
	drop if p > `report_p'
	merge 1:1 label using `summaryqc_label'
	keep if _m == 3
	compress
	drop _m
	gsort -beta
	egen obs = seq()
	sum obs if beta <0
	gen yline = `r(min)' - 0.5
	sum l
	gen l2=`r(min)'
	replace l2 = round(l2 -.4,.02)
	gen obs1 = obs -.3
	gen obs2 = obs
	sum yline
	local yline `r(min)' 
	replace label2 = class + ": " + label2
	qui { // define the 22 colors in the class2
		/*
		local color1 "178 223 138"
		local color2 "239 251 253"
		local color3 "240 002 127"
		local color4 "102 102 102"
		local color5 "177 100 154"
		local color6 "000 068 027"
		local color7 "165 015 021"
		local color8 "255 127 000"
		local color9 "008 048 107"
		local color10 "247 129 191"
		local color11 "251 154 153"
		local color12 "000 000 000"
		local color13 "166 054 003"
		local color14 "031 120 180"
		local color15 "255 255 153"
		local color16 "166 206 227"
		local color17 "227 026 028"
		local color18 "106 061 154"
		local color19 "103 000 013"
		local color20 "253 191 111"
		local color21 "051 160 044"
		local color22 "202 178 214"
		*/
		}
	count
	gen report_n = `r(N)'
	if report_n < `report_full' {
		#delimit;
		tw rcap l u obs, horizontal lw(vthin) xlabel(-.3(.05).3)
		|| scatter obs beta [w=n] ,  msymbol(s)  msize(tiny) mfc(red) mlc(black) mlw(vvthin)
		|| scatter obs2 l2,   msymbol(i) mlabel(label2) mlabpos(3) mlabcolor(black) mlabsize(tiny) 
		|| scatter obs1 l2,   msymbol(i) mlabel(lnor)   mlabpos(3) mlabcolor(black)  mlabsize(tiny)
		xline(0) yline(`yline') 
		ylabel(none) ytitle("") legend(off) xtitle("beta (+/- std.err)") 
		nodraw saving(`index'_`function'.gph, replace)
		title("PheWAS for `index'") 
		;
		#delimit cr
		}
	else {
		#delimit;
		tw rcap l u obs, horizontal lw(vthin) xlabel(-.3(.05).3)
		|| scatter obs beta [w=n] ,  msymbol(s)  msize(tiny) mfc(red) mlc(black) mlw(vvthin)
		|| scatter obs l2,   msymbol(i) mlabel(label2) mlabpos(3) mlabcolor(black) mlabsize(tiny) 
		xline(0) yline(`yline') 
		ylabel(none) ytitle("") legend(off) xtitle("beta (+/- std.err)") 
		nodraw saving(`index'_`function'.gph, replace)
		title("PheWAS for `index'") 
		;
		#delimit cr
		}
	noi di as text"# > .... plot graph to "as result"`index'_`function'.gph"
	clear 
	set obs 1
	gen png = "`png'"
	if png == "yes" {
		graph use `index'_`function'.gph
		noi di as text"# > ..... saved png to "as result"`index'_`function'.png"
		graph export `index'_`function'.png, as(png) height(2000) width(4000) replace
		window manage close graph
		}
	else {
		}
	noi di as text"#########################################################################"
	noi di as text""
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 4 - write log file
	file open myfile using "`index'-`function'.log", write replace
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
	file write myfile `"# > ........ index SNP `index'"' _n
	file write myfile `"# > ....... saved graph `index'-`function'.gph"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > overview of `function'"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > 1. parse index from all or defined set of GWAS files"' _n
	file write myfile `"# > 2. limit to P < `report_p'"' _n
	file write myfile `"# > 3. plot beta / se"' _n
	file write myfile `"#########################################################################"' _n
	file close myfile	
	}
qui { // module 5 - clean folder
	files2dta, dir(`c(pwd)')
	gen keep = .
	replace keep = 1 if file == "`index'-final-`function'.dta"
	replace keep = 1 if file == "`index'-`function'.log"
	replace keep = 1 if file == "`index'_`function'.png"
	replace keep = 1 if file == "`index'_`function'.gph"
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

	
