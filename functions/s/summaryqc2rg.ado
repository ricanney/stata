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
		bim_path_to_short, bim(`=_var1[`i']')
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
		local already_created ${init_root}/summaryqc2sumstats/``function'_summaryqc`i'_short'/``function'_summaryqc`i'_short'-summaryqc2sumstats-hw3.sumstats
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
	gen     sumstat = "${init_root}/summaryqc2sumstats/" + summaryqc + "/" + summaryqc + "-summaryqc2sumstats-hw3.sumstats,"
	replace sumstat = "${init_root}/summaryqc2sumstats/" + summaryqc + "/" + summaryqc + "-summaryqc2sumstats-hw3.sumstats" in ``function'_summaryqc_N'
	sxpose, clear
	drop in 1
	gen _var999 = ""
	gen _var = ""
	for var _var1-_var999: replace _var = _var + X
	!python `ldsc' --rg `=_var[1]' --ref-ld `w_hm3'w_hm3/w_hm3 --w-ld `w_hm3'w_hm3/w_hm3 --out ``function'_summaryqc1_short'-`function'-hw3-rg
	import delim using  ``function'_summaryqc1_short'-`function'-hw3-rg.log, varnames(nonames) clear		
	split v1,p("${init_root}""/summaryqc2sumstats/")
	keep if v13 !=""
	keep v1 
	split v1,p("--rg ")
	drop if v11 == ""
	keep v1
	split v1,p("${init_root}""/summaryqc2sumstats/")
	keep v13 v15
	drop if v15 == ""
	split v13,p("/")
	split v15,p("/")
	foreach i of num 1/10 {
		replace v152 = subinstr(v152,"  "," ",.)
		}
	split v152,p(" ")
	rename (v131 v151 v1522-v15211) (p1 p2 rg se z p h2_obs h2_obs_se h2_int h2_int_se gcov_int gcov_int_se)
	keep p1 p2 rg se z p h2_obs h2_obs_se h2_int h2_int_se gcov_int gcov_int_se
	for var rg se z p h2_obs h2_obs_se h2_int h2_int_se gcov_int gcov_int_se: destring X, replace force
	noi di as text"# > ...... index study "as result`"``function'_summaryqc1_short'"'
	noi di as text"#########################################################################"
	count
	foreach i of num 1 / `r(N)' {
			noi di as text"# > . comparison study "as result"`=p2[`i']'"
			noi di as text"# > ............... rg "as result"`: display %05.3f `=rg[`i']'' (`: display %05.3f `=se[`i']'')"
			noi di as text"# > ................ z "as result"`: display %05.3f `=z[`i']''"
			noi di as text"# > ... gcov intercept "as result"`: display %05.3f `=gcov_int[`i']'' (`: display %05.3f `=gcov_int_se[`i']'')"
			noi di as text"# > ................ p "as result"`: display %10.4e `=p[`i']''"
			noi di as text"#########################################################################"
			}
	lab var p1 					"Primary GWAS"
	lab var p2 					"Comparison GWAS"
	lab var rg 					"Genetic Correlation"
	lab var se 					"Standard Error (rg)"
	lab var z  					"Z-score"
	lab var p  					"P-value"
	lab var h2_obs 			"SNP based heritability (Comparison GWAS)"
	lab var h2_obs_se 	"Standard Error (h2)"
	lab var h2_int 			"SNP based hertiability intercept (Comparison GWAS)"
	lab var h2_int_se 	"Standard Error (h2_int)"
	lab var gcov_int 		"Genetic Covariance intercept"
	lab var gcov_int_se	"Standard Error (gcov_int)"
	save ``function'_summaryqc1_short'-`function'.dta, replace
	noi di as text""
	}
qui { // module 5 - plot correlation graph
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 5 - plot correlation graph"
	noi di as text"#########################################################################"
	gen sig = .
	replace sig = 1 if p < .05
	gen upper = rg + se
	gen lower = rg - se
	replace rg = 1 if rg > 1
	replace rg = -1 if rg < -1
	replace upper = 1 if upper > 1
	replace lower = -1 if lower < -1
	split p2,p("-2005-""-2006-""-2007-""-2008-""-2009-""-2010-""-2011-""-2012-""-2013-""-2014-""-2015-""-2016-""-2017-""-2018-""-2019-""-2020-")
	drop if p1 == p2
	for var rg se : drop if X == .
	encode p22, gen(trait)
	count
	gen max = `r(N)' + .5
	sum max
	twoway 	rspike lower upper trait, lc(gray) horizontal || ///
	scatter trait rg, msymbol(O) mcolor(gray) mlcolor(black) msize(medium) || ///
	scatter trait rg if sig == 1, msymbol(O) mcolor(red) mlcolor(black) msize(large) ///
	legend(off) ///
	ytitle(" ") ylabel(0.5 " " 1(1)`r(N)' `r(mean)' " ",valuelabel nogrid angle(0) labs(small)) ///
	xtitle("Genetic Correlation (rg)") xlabel(-1.1 " " -1(.2)1, nogrid) xline(0) ///
	title("``function'_summaryqc1_short'") nodraw saving(``function'_summaryqc1_short'-`function'.gph, replace)
	noi checkfile, file(``function'_summaryqc1_short'-`function'.gph)
	noi di as text""
	}
qui { // module 5 - write log file
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
		file write myfile `"# > ............... rg `: display %05.3f `=rg[`i']'' (`: display %05.3f `=se[`i']'')"' _n
		file write myfile `"# > ................ z `: display %05.3f `=z[`i']''"' _n
		file write myfile `"# > ................ p `: display %10.4e `=p[`i']''"' _n
		file write myfile `"# > ... gcov intercept `: display %05.3f `=gcov_int[`i']'' (`: display %05.3f `=gcov_int_se[`i']'')"' _n
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
qui { // module 6 - clean folder
	files2dta, dir(`c(pwd)')
	gen keep = .
	replace keep = 1 if file == "``function'_summaryqc1_short'-`function'-hw3-rg.log"
	replace keep = 1 if file == "``function'_summaryqc1_short'-`function'.log"
	replace keep = 1 if file == "``function'_summaryqc1_short'-`function'.dta"
	replace keep = 1 if file == "``function'_summaryqc1_short'-`function'.gph"
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

	
