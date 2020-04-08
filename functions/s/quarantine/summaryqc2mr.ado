program  summaryqc2mr
local    function summaryqc2mr
local    `function'_version v1
syntax , summaryqc1(string asis) summaryqc2(string asis) ldref(string asis) glist_hg19(string asis) [p_max(real 1e-5) first(string asis)]

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
	ref_path_to_short, ref(`ldref')
	local `function'_ldref_short ${ref_short}
	foreach i in 1 2 {
		bim_path_to_short, bim(`summaryqc`i'')
		local `function'_summaryqc`i' ${bim_short}
		}
	local `function'_project ``function'_summaryqc1'-``function'_summaryqc2'
	cd ${init_root}	
	!mkdir `function'
	cd `function'
	noi di as text"# > ..... project name "as result"``function'_project'"
	!mkdir ``function'_project'
	cd ``function'_project'
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 2 - define files
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 2: define files"
	noi di as text"#########################################################################"
	foreach i in 1 2 { 
		noi checkfile, file(`summaryqc`i''-summaryqc.dta)
		}
	checktabbed
	clear
	set obs 1
	gen first = "`first'"
	if first == "yes" {
		noi di as text"# > .................. load MR packages"
		noi net install mrrobust, from(https://raw.github.com/remlapmot/mrrobust/master/)
		noi ssc install addplot
		noi ssc install moremata
		noi ssc install heterogi
		noi ssc install kdens
		noi ssc install metan
		}
	else if first == "update" {
		adoupdate mrrobust, update
		adoupdate addplot, update
		adoupdate moremata, update
		adoupdate heterogi, update
		adoupdate kdens, update
		adoupdate metan, update
		}
	noi di as text"#########################################################################"
	noi di as text""
	}
qui { // module 3 - define common SNPs
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3: define common SNPs"
	noi di as text"#########################################################################"
	local already_created ${init_root}/`function'/``function'_project'/``function'_summaryqc1'-`function'-summaryqc.dta
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(${init_root}/`function'/``function'_project'/``function'_summaryqc1'-`function'-summaryqc.dta)
		noi checkfile, file(${init_root}/`function'/``function'_project'/``function'_summaryqc2'-`function'-summaryqc.dta)
		noi checkfile, file(${init_root}/`function'/``function'_project'/``function'_project'-`function'-common-snps.dta)
		}
	else {
		noi di as text"# >  define SNPs using "as result"``function'_summaryqc1'.dta"
		use snp a1 a2 using `summaryqc1'-summaryqc.dta, clear
		recodegenotype
		drop if recode == "W"
		drop if recode == "S"
		keep snp
		noi di as text"# > ...... merge using "as result"``function'_summaryqc2'.dta"
		merge 1:1 snp using `summaryqc2'-summaryqc.dta
		keep if _m == 3
		drop _m
		save ``function'_summaryqc2'-`function'-summaryqc.dta, replace
		keep snp a1 a2
		rename (a1 a2) (b1 b2)
		save ``function'_project'-`function'-common-snps.dta, replace
		noi di as text"# > ........... import "as result"``function'_summaryqc1'.dta"
		use `summaryqc1'-summaryqc.dta, clear
		merge 1:1 snp using ``function'_project'-`function'-common-snps.dta
		keep if _m == 3
		drop _m
		noi di as text"# >  convert strand to "as result"``function'_summaryqc2'.dta"
		gen flip = 1 if a1 != b1
		replace beta = -beta if flip == 1
		replace z    = - z   if flip == 1 
		for var or l95 u95: replace X = 1/X if flip == 1
		gen x95 = u95
		replace u95 = l95 if flip == 1
		replace l95 = x95 if flip == 1
		drop x95
		keep chr bp snp a1 a2 beta se z or l95 u95 p n
		save ``function'_summaryqc1'-`function'-summaryqc.dta, replace
		noi checkfile, file(${init_root}/`function'/``function'_project'/``function'_summaryqc1'-`function'-summaryqc.dta)
		noi checkfile, file(${init_root}/`function'/``function'_project'/``function'_summaryqc2'-`function'-summaryqc.dta)
		noi checkfile, file(${init_root}/`function'/``function'_project'/``function'_project'-`function'-common-snps.dta)
		}
	noi di as text"#########################################################################"
	noi di as text""
	}
qui { // module 4 - run summaryqc2top
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 4: run summaryqc2top"
	noi di as text"#########################################################################"
	foreach i of num 1 2 {
		noi di as text"# > .. define top SNPs "as result"``function'_summaryqc`i''.dta"
		local already_created ${init_root}/summaryqc2top/``function'_summaryqc`i''-`function'/``function'_summaryqc`i''-`function'-summaryqc2top.dta
		capture confirm file `already_created'
		if !_rc {
			noi checkfile, file(`already_created')
			}
		else {
			noi summaryqc2top, summaryqc(${init_root}/`function'/``function'_project'/``function'_summaryqc`i''-`function') ldref(${summaryqc2top_ref}) glist_hg19(${summaryqc2top_glist})
			}
		}
	noi di as text"#########################################################################"
	noi di as text""
	}
qui { // module 5 - define instrument variables (LD independent SNPs P < `: display %10.4e `p_max'')
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 5: define instrument variables (LD independent SNPs P < `: display %10.4e `p_max'')"
	noi di as text"#########################################################################"
	foreach i of num 1 2 {
		local already_created "${init_root}/`function'/``function'_project'/``function'_summaryqc`i''-`function'-iv.dta"
		capture confirm file `already_created'
		if !_rc {
			noi checkfile, file(`already_created')
			}
		else {		
			use ${init_root}/summaryqc2top/``function'_summaryqc`i''-`function'/``function'_summaryqc`i''-`function'-summaryqc2top.dta, clear
			keep if P < `p_max'
			keep SNP
			renvars, lower
			merge 1:1 snp using ${init_root}/`function'/``function'_project'/``function'_summaryqc`i''-`function'-summaryqc.dta
			keep if _m == 3
			drop _m
			save "${init_root}/`function'/``function'_project'/``function'_summaryqc`i''-`function'-iv.dta", replace
			noi checkfile, file(`already_created')
			}
		}
	noi di as text"#########################################################################"
	noi di as text""
	}
qui { // module 6 - run regressions
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 6: run regressions"
	noi di as text"#########################################################################"
	qui { // run model 1 
		local already_created "``function'_summaryqc1'-`function'-mreggerplot.png"
		capture confirm file `already_created'
		if !_rc {
			noi checkfile, file(`already_created')
			}
		else {
			qui { // merge iv with paired GWAS
				use "``function'_summaryqc1'-`function'-iv.dta", clear
				keep snp beta se 
				rename (beta se) (beta_1 se_1)
				merge 1:1 snp using ${init_root}/`function'/``function'_project'/``function'_summaryqc2'-`function'-summaryqc.dta
				keep if _m == 3
				keep snp beta_1 se_1 beta se
				rename (beta se) (beta_2 se_2)
				* covert all exposure variable to +ve
				gen flip = 1 if beta_1 <0
				replace beta_1 = -beta_1 if flip == 1
				replace beta_2 = -beta_2 if flip == 1
				drop flip
				}
			noi { // calculate regressions
				log using "``function'_summaryqc1'-`function'-mregger.log", replace
				noi di as text"# > ............ model "as result"(1) inverse weighted regression"
				noi mregger beta_2 beta_1 [aw=1/(se_2^2)], ivw
				noi di as text"# > ............ model "as result"(2) inverse weighted regression (fixed-effect)"
				noi mregger beta_2 beta_1 [aw=1/(se_2^2)], ivw fe
				noi di as text"# > ............ model "as result"(3) median regression "
				noi mrmedian beta_2 se_2 beta_1 se_1 
				noi di as text"# > ............ model "as result"(4) median regression (weighted) "
				noi mrmedian beta_2 se_2 beta_1 se_1, weighted
				noi di as text"# > ............ model "as result"(5) median regression (penalised-weighted)"
				noi mrmedian beta_2 se_2 beta_1 se_1, penweighted
				noi di as text"# > ............ model "as result"(6) modal regression"
				noi mrmodal beta_2 se_2 beta_1 se_1, phi(2)
				noi di as text"# > ............ model "as result"(7) modal regression (weighted)"
				noi mrmodal beta_2 se_2 beta_1 se_1, weighted phi(2)
				noi di as text"# > ............ model "as result"(8) modal regression (NOME)"
				noi mrmodal beta_2 se_2 beta_1 se_1, nome phi(2)
				noi di as text"# > ............ model "as result"(9) modal regression (weighted/NOME)"
				noi mrmodal beta_2 se_2 beta_1 se_1, weighted nome phi(2)
				noi di as text"# > ............ model "as result"(10) mr-egger regression reporting I^2_GX statistic and heterogeneity Q-test"
				noi mregger beta_2 beta_1 [aw=1/(se_2^2)], ivw gxse(se_1) heterogi
				noi di as text"# > ............ model "as result"(11) mr-egger simex"
				noi mreggersimex beta_2 beta_1 [aw=1/(se_2^2)], gxse(se_1) seed(12345) nodraw
				log close
				}
			qui { // plot regressions 
				sum beta_1
				global max = r(max)*1.1
				global min = r(min)/1.1
				mreggerplot beta_2 se_2 beta_1 se_1, xtitle(" ""Exposure" "``function'_summaryqc1'") ytitle("Outcome" "``function'_summaryqc2'" " ") 
				mregger beta_2 beta_1 [aw=1/(se_2^2)], ivw 
				addplot : function _b[beta_1]*x  , range($min $max) lw(vthin) lp(dot) lc("198 219 239")
				mregger beta_2 beta_1 [aw=1/(se_2^2)], ivw fe
				addplot : function _b[beta_1]*x  , range($min $max) lw(vthin) lp(dot) lc("008 069 148")
				mrmedian beta_2 se_2 beta_1 se_1 
				addplot : function _b[beta]*x  , range($min $max) lw(vthin) lp(dash) lc("199 233 192")
				mrmedian beta_2 se_2 beta_1 se_1, weighted
				addplot : function _b[beta]*x , range($min $max) lw(vthin) lp(dash) lc("065 171 193")
				mrmedian beta_2 se_2 beta_1 se_1, penweighted
				addplot : function _b[beta]*x , range($min $max) lw(vthin) lp(dash) lc("000 090 050")
				mrmodal beta_2 se_2 beta_1 se_1, phi(2)
				addplot : function _b[beta]*x , range($min $max) lw(vthin) lp(dash_dot) lc("252 187 161")
				mrmodal beta_2 se_2 beta_1 se_1, weighted phi(2)
				addplot : function _b[beta]*x , range($min $max) lw(vthin) lp(dash_dot) lc("252 146 073")
				mrmodal beta_2 se_2 beta_1 se_1, nome phi(2)
				addplot : function _b[beta]*x , range($min $max) lw(vthin) lp(dash_dot) lc("239 059 044")
				mrmodal beta_2 se_2 beta_1 se_1, weighted nome phi(2)
				addplot : function _b[beta]*x , range($min $max) lw(vthin) lp(dash_dot) lc("153 000 013") ///
				legend(order(5 "Instruments" 4 "95% CIs" 3 "MR-Egger" 2 "MR-Egger 95% CI" 6 "IVW" 7 "IVW-FE" 8 "Median" 9 "Median (weighted)" 10 "Median (peanlised-weighted)" 11 "Modal" 12 "Modal (weighted)" 13 "Modal (NOME)" 14 "Modal (NOME/weighted)") ///
				rows(5) size(tiny))
				set scheme plotplainblind , permanently
				graph export "``function'_summaryqc1'-`function'-mreggerplot.png", as(png) height(2000) width(3000) replace
				}
			}
		}
	qui { // run model 2 
		local already_created "``function'_summaryqc2'-`function'-mreggerplot.png"
		capture confirm file `already_created'
		if !_rc {
			noi checkfile, file(`already_created')
			}
		else {
			qui { // merge iv with paired GWAS
				use "``function'_summaryqc2'-`function'-iv.dta", clear
				keep snp beta se 
				rename (beta se) (beta_1 se_1)
				merge 1:1 snp using ${init_root}/`function'/``function'_project'/``function'_summaryqc1'-`function'-summaryqc.dta
				keep if _m == 3
				keep snp beta_1 se_1 beta se
				rename (beta se) (beta_2 se_2)
				* covert all exposure variable to +ve
				gen flip = 1 if beta_1 <0
				replace beta_1 = -beta_1 if flip == 1
				replace beta_2 = -beta_2 if flip == 1
				drop flip
				}
			noi { // calculate regressions
				log using "``function'_summaryqc2'-`function'-mregger.log", replace
				noi di as text"# > ............ model "as result"(1) inverse weighted regression"
				noi mregger beta_2 beta_1 [aw=1/(se_2^2)], ivw
				noi di as text"# > ............ model "as result"(2) inverse weighted regression (fixed-effect)"
				noi mregger beta_2 beta_1 [aw=1/(se_2^2)], ivw fe
				noi di as text"# > ............ model "as result"(3) median regression "
				noi mrmedian beta_2 se_2 beta_1 se_1 
				noi di as text"# > ............ model "as result"(4) median regression (weighted) "
				noi mrmedian beta_2 se_2 beta_1 se_1, weighted
				noi di as text"# > ............ model "as result"(5) median regression (penalised-weighted)"
				noi mrmedian beta_2 se_2 beta_1 se_1, penweighted
				noi di as text"# > ............ model "as result"(6) modal regression"
				noi mrmodal beta_2 se_2 beta_1 se_1, phi(2)
				noi di as text"# > ............ model "as result"(7) modal regression (weighted)"
				noi mrmodal beta_2 se_2 beta_1 se_1, weighted phi(2)
				noi di as text"# > ............ model "as result"(8) modal regression (NOME)"
				noi mrmodal beta_2 se_2 beta_1 se_1, nome phi(2)
				noi di as text"# > ............ model "as result"(9) modal regression (weighted/NOME)"
				noi mrmodal beta_2 se_2 beta_1 se_1, weighted nome phi(2)
				noi di as text"# > ............ model "as result"(10) mr-egger regression reporting I^2_GX statistic and heterogeneity Q-test"
				noi mregger beta_2 beta_1 [aw=1/(se_2^2)], ivw gxse(se_1) heterogi
				noi di as text"# > ............ model "as result"(11) mr-egger simex"
				noi mreggersimex beta_2 beta_1 [aw=1/(se_2^2)], gxse(se_1) seed(12345) nodraw
				log close
				}
			qui { // plot regressions 
				sum beta_1
				global max = r(max)*1.1
				global min = r(min)/1.1
				mreggerplot beta_2 se_2 beta_1 se_1, xtitle(" ""Exposure" "``function'_summaryqc1'") ytitle("Outcome" "``function'_summaryqc2'" " ") 
				mregger beta_2 beta_1 [aw=1/(se_2^2)], ivw 
				addplot : function _b[beta_1]*x  , range($min $max) lw(vthin) lp(dot) lc("198 219 239")
				mregger beta_2 beta_1 [aw=1/(se_2^2)], ivw fe
				addplot : function _b[beta_1]*x  , range($min $max) lw(vthin) lp(dot) lc("008 069 148")
				mrmedian beta_2 se_2 beta_1 se_1 
				addplot : function _b[beta]*x  , range($min $max) lw(vthin) lp(dash) lc("199 233 192")
				mrmedian beta_2 se_2 beta_1 se_1, weighted
				addplot : function _b[beta]*x , range($min $max) lw(vthin) lp(dash) lc("065 171 193")
				mrmedian beta_2 se_2 beta_1 se_1, penweighted
				addplot : function _b[beta]*x , range($min $max) lw(vthin) lp(dash) lc("000 090 050")
				mrmodal beta_2 se_2 beta_1 se_1, phi(2)
				addplot : function _b[beta]*x , range($min $max) lw(vthin) lp(dash_dot) lc("252 187 161")
				mrmodal beta_2 se_2 beta_1 se_1, weighted phi(2)
				addplot : function _b[beta]*x , range($min $max) lw(vthin) lp(dash_dot) lc("252 146 073")
				mrmodal beta_2 se_2 beta_1 se_1, nome phi(2)
				addplot : function _b[beta]*x , range($min $max) lw(vthin) lp(dash_dot) lc("239 059 044")
				mrmodal beta_2 se_2 beta_1 se_1, weighted nome phi(2)
				addplot : function _b[beta]*x , range($min $max) lw(vthin) lp(dash_dot) lc("153 000 013") ///
				legend(order(5 "Instruments" 4 "95% CIs" 3 "MR-Egger" 2 "MR-Egger 95% CI" 6 "IVW" 7 "IVW-FE" 8 "Median" 9 "Median (weighted)" 10 "Median (peanlised-weighted)" 11 "Modal" 12 "Modal (weighted)" 13 "Modal (NOME)" 14 "Modal (NOME/weighted)") ///
				rows(5) size(tiny))
				set scheme plotplainblind , permanently
				graph export "``function'_summaryqc2'-`function'-mreggerplot.png", as(png) height(2000) width(3000) replace
				}
			}
		}
	noi di as text"#########################################################################"
	noi di as text""
	}
qui { // module 7 - extract information from log files
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 7: extract information from log files"
	noi di as text"#########################################################################"	
	qui { // 1>2
		local `function'_exposure ``function'_summaryqc1'
		local `function'_outcome  ``function'_summaryqc2'
		!$tabbed ``function'_exposure'-summaryqc2mr-mregger.log
		import delim using ``function'_exposure'-summaryqc2mr-mregger.log.tabbed, clear
		dropmiss, force
		gen v99 = ""
		egen a = concat(v1-v99),p(" ") 
		dropmiss, force
		aorder
		replace v1 = subinstr(v1,"-","",.)
		replace v1 = subinstr(v1,":","",.)
		drop if v1 == ">"
		drop if v1 == ""
		drop if v1 == "+"
		drop if v1 == "("
		gen obs = _n
		gen n = v5 if v1 == "Number"
		gen replications = v3 if v1 == "Replications"
		for var n replication : replace X = X[_n-1] if X == ""
		gen exposure = "``function'_exposure'"
		gen outcome  = "``function'_outcome'"
		qui { // model
			gen model = ""
			replace model = "ivw-fe"                           if a == "# > ............ model (1) inverse weighted regression"
			replace model = "ivw"                              if a == "# > ............ model (2) inverse weighted regression (fixed-effect)"
			replace model = "median"                           if a == "# > ............ model (3) median regression" 
			replace model = "median (weighted)"                if a == "# > ............ model (4) median regression (weighted)" 
			replace model = "median (penalised-weighted)"      if a == "# > ............ model (5) median regression (penalised-weighted)" 
			replace model = "modal"                            if a == "# > ............ model (6) modal regression" 
			replace model = "modal (weighted)"                 if a == "# > ............ model (7) modal regression (weighted)" 
			replace model = "modal (penalised-weighted)"       if a == "# > ............ model (8) modal regression (NOME)" 
			replace model = "modal (NOME)"                     if a == "# > ............ model (9) modal regression (weighted/NOME)" 
			replace model = "mr-egger (slope)"                 if a == "# > ............ model (10) mr-egger regression reporting I^2_GX statistic and heterogeneity Q-test" 
			replace model = "mr-egger-simex (slope)"           if a == "# > ............ model (11) mr-egger simex" 
			replace model = model[_n-1] if model == ""
			drop if model == ""
			}
		qui { //  coef 
			gen keep = ""
			replace keep = "Y" if v1 == "beta"
			replace keep = "Y" if v1 == "beta_1"
			replace keep = "Y" if v1 == "slope"
			replace keep = "Y" if v1 == "_cons"
			keep if keep == "Y"
			}
		rename (v3 v4 v5 v6 v7 v8) (coef se z p l95 u95)
		for var model coef se z p l95 u95: replace X = X[_n-1] if X == ""
		keep  exposure outcome model n replications coef se z p l95 u95
		order exposure outcome model n replications coef se z p l95 u95
		replace replications = "NA" if replications == ""
		replace replications = "NA" if model == "mr-egger (slope)"
		replace replications = "NA" if model == "mr-egger-simex (slope)"
		egen x = seq(),by(model)
		replace model = "mr-egger-simex (_cons)"           if x == 2
		drop x 
		compress
		replace p = "<0.0001" if p == "0.000"
		save ``function'_exposure'-summaryqc2mr-mregger.dta, replace
		noi checkfile, file(``function'_exposure'-summaryqc2mr-mregger.dta)
		import delim using ``function'_exposure'-summaryqc2mr-mregger.log.tabbed, clear
		gen keep = ""
		foreach keep in H I^2 Q-test {
		replace keep = "Y" if v1 == "`keep'"
		}
		keep if keep == "Y"
		compress
		gen H   		= v3 + " (" + v4 + "-" + v5 + ")" in 1
		gen I_2 		= v3 + " (" + v4 + "-" + v5 + ")" in 2
		gen     Q_test 	= v3 + " (d.f. = " + v6 + "; p = " + v9 + ")" in 3
		gen     Q_test2 	= v3 + " (d.f. = " + v6 + "; p < 0.0001)" in 3
		gen p = v9 in 3
		replace Q_test = Q_test2 if p == "0.0000"
		keep H I_2 Q_test 
		for var H I_2 Q_test: replace X = X[_n-1] if X == ""
		for var H I_2 Q_test: drop if X == ""
		duplicates drop
		sxpose, clear
		gen exposure = "``function'_exposure'"
		gen outcome  = "``function'_outcome'"		
		gen model = ""
		replace model = "H" in 1
		replace model = "I^2" in 2
		replace model = "Q-test" in 3
		rename _var1 HETEROGENEITY
		append using ``function'_exposure'-summaryqc2mr-mregger.dta
		renvars, upper
		order EXPOSURE OUTCOME MODEL HETERO
		save ``function'_exposure'-summaryqc2mr-mregger.dta, replace
		}
	qui { // 2>1
		local `function'_exposure ``function'_summaryqc2'
		local `function'_outcome  ``function'_summaryqc1'
		!$tabbed ``function'_exposure'-summaryqc2mr-mregger.log
		import delim using ``function'_exposure'-summaryqc2mr-mregger.log.tabbed, clear
		dropmiss, force
		gen v99 = ""
		egen a = concat(v1-v99),p(" ") 
		dropmiss, force
		aorder
		replace v1 = subinstr(v1,"-","",.)
		replace v1 = subinstr(v1,":","",.)
		drop if v1 == ">"
		drop if v1 == ""
		drop if v1 == "+"
		drop if v1 == "("
		gen obs = _n
		gen n = v5 if v1 == "Number"
		gen replications = v3 if v1 == "Replications"
		for var n replication : replace X = X[_n-1] if X == ""
		gen exposure = "``function'_exposure'"
		gen outcome  = "``function'_outcome'"
		qui { // model
			gen model = ""
			replace model = "ivw-fe"                           if a == "# > ............ model (1) inverse weighted regression"
			replace model = "ivw"                              if a == "# > ............ model (2) inverse weighted regression (fixed-effect)"
			replace model = "median"                           if a == "# > ............ model (3) median regression" 
			replace model = "median (weighted)"                if a == "# > ............ model (4) median regression (weighted)" 
			replace model = "median (penalised-weighted)"      if a == "# > ............ model (5) median regression (penalised-weighted)" 
			replace model = "modal"                            if a == "# > ............ model (6) modal regression" 
			replace model = "modal (weighted)"                 if a == "# > ............ model (7) modal regression (weighted)" 
			replace model = "modal (penalised-weighted)"       if a == "# > ............ model (8) modal regression (NOME)" 
			replace model = "modal (NOME)"                     if a == "# > ............ model (9) modal regression (weighted/NOME)" 
			replace model = "mr-egger (slope)"                 if a == "# > ............ model (10) mr-egger regression reporting I^2_GX statistic and heterogeneity Q-test" 
			replace model = "mr-egger-simex (slope)"           if a == "# > ............ model (11) mr-egger simex" 
			replace model = model[_n-1] if model == ""
			drop if model == ""
			}
		qui { //  coef 
			gen keep = ""
			replace keep = "Y" if v1 == "beta"
			replace keep = "Y" if v1 == "beta_1"
			replace keep = "Y" if v1 == "slope"
			replace keep = "Y" if v1 == "_cons"
			keep if keep == "Y"
			}
		rename (v3 v4 v5 v6 v7 v8) (coef se z p l95 u95)
		for var model coef se z p l95 u95: replace X = X[_n-1] if X == ""
		keep  exposure outcome model n replications coef se z p l95 u95
		order exposure outcome model n replications coef se z p l95 u95
		replace replications = "NA" if replications == ""
		replace replications = "NA" if model == "mr-egger (slope)"
		replace replications = "NA" if model == "mr-egger-simex (slope)"
		egen x = seq(),by(model)
		replace model = "mr-egger-simex (_cons)"           if x == 2
		drop x 
		compress
		replace p = "<0.0001" if p == "0.000"
		save ``function'_exposure'-summaryqc2mr-mregger.dta, replace
		noi checkfile, file(``function'_exposure'-summaryqc2mr-mregger.dta)
		import delim using ``function'_exposure'-summaryqc2mr-mregger.log.tabbed, clear
		gen keep = ""
		foreach keep in H I^2 Q-test {
		replace keep = "Y" if v1 == "`keep'"
		}
		keep if keep == "Y"
		compress
		gen H   		= v3 + " (" + v4 + "-" + v5 + ")" in 1
		gen I_2 		= v3 + " (" + v4 + "-" + v5 + ")" in 2
		gen     Q_test 	= v3 + " (d.f. = " + v6 + "; p = " + v9 + ")" in 3
		gen     Q_test2 	= v3 + " (d.f. = " + v6 + "; p < 0.0001)" in 3
		gen p = v9 in 3
		replace Q_test = Q_test2 if p == "0.0000"
		keep H I_2 Q_test 
		for var H I_2 Q_test: replace X = X[_n-1] if X == ""
		for var H I_2 Q_test: drop if X == ""
		duplicates drop
		sxpose, clear
		gen exposure = "``function'_exposure'"
		gen outcome  = "``function'_outcome'"		
		gen model = ""
		replace model = "H" in 1
		replace model = "I^2" in 2
		replace model = "Q-test" in 3
		rename _var1 HETEROGENEITY
		append using ``function'_exposure'-summaryqc2mr-mregger.dta
		renvars, upper
		order EXPOSURE OUTCOME MODEL HETERO
		save ``function'_exposure'-summaryqc2mr-mregger.dta, replace		
		}
	append using ``function'_outcome'-summaryqc2mr-mregger.dta
	save ``function'_project'-summaryqc2mr-mregger.dta, replace
	noi checkfile, file(``function'_project'-summaryqc2mr-mregger.dta)
	noi di as text"#########################################################################"
	noi di as text""
	}
qui { // module y - write log file
	file open myfile using "``function'_project'-`function'.log", write replace
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
	file write myfile `"# > ...... summaryqc1 ``function'_summaryqc1'"' _n
	file write myfile `"# > ...... summaryqc2 ``function'_summaryqc2'"' _n
	file write myfile `"# > ....... reference ``function'_ldref_short'"' _n
	file write myfile `"# > ....... gene-list ``function'_glist_hg19_short'"' _n
	file write myfile `"# > ................. reporting thresholds"' _n
	file write myfile `"# > ........ iv p-max `: display %10.4e `p_top''"' _n
	file write myfile`"#########################################################################"' _n
	file write myfile `"# > ...... saved logs ``function'_summaryqc1'-summaryqc2mr-mregger.log"' _n
	file write myfile `"# > .... saved graphs ``function'_summaryqc1'-summaryqc2mr-mreggerplot.png"' _n
	file write myfile `"# > .... saved tables ``function'_summaryqc1'-summaryqc2mr-mregger.dta"' _n
	file write myfile`"#########################################################################"' _n
	file write myfile `"# > ...... saved logs ``function'_summaryqc2'-summaryqc2mr-mregger.log"' _n
	file write myfile `"# > .... saved graphs ``function'_summaryqc2'-summaryqc2mr-mreggerplot.png"' _n
	file write myfile `"# > .... saved tables ``function'_summaryqc2'-summaryqc2mr-mregger.dta"' _n	
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > overview of `function'"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > 1. define common snps across gwas"' _n
	file write myfile `"# > 2. create top ld-independent snps per parsed gwas (summaryqc2top)"' _n
	file write myfile `"# > 3. define instrument variables per GWAS at P< `: display %10.4e `p_max''"' _n
	file write myfile `"# > 4. join with outcome"' _n
	file write myfile `"# > 5. perform MR regressions (x11) "' _n
	file write myfile `"# > 6. plot regressions"' _n
	file write myfile `"# > 7. create table"' _n
	file write myfile `"#########################################################################"' _n
	file close myfile	
	}
qui { // module z - clean folder
	files2dta, dir(`c(pwd)')
	gen keep = .
	replace keep = 1 if file == "``function'_project'-`function'.log"
	replace keep = 1 if file == "``function'_project'-`function'-mregger.dta"
	replace keep = 1 if file == "``function'_summaryqc1'-`function'-mregger.log"
	replace keep = 1 if file == "``function'_summaryqc2'-`function'-mregger.log"
	replace keep = 1 if file == "``function'_summaryqc1'-`function'-mreggerplot.png"
	replace keep = 1 if file == "``function'_summaryqc2'-`function'-mreggerplot.png"
	replace keep = 1 if file == "``function'_summaryqc1'-`function'-iv.dta"
	replace keep = 1 if file == "``function'_summaryqc2'-`function'-iv.dta"
	replace keep = 1 if file == "``function'_summaryqc1'-`function'-summaryqc.dta"
	replace keep = 1 if file == "``function'_summaryqc2'-`function'-summaryqc.dta"
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






	


