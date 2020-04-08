*! 1.0.2 Richard Anney 17apr2019
* clean up-front functions
*! 1.0.1 Richard Anney 24jan2019
* adaptation from graphplinkkin0
* - update to start from bim file
* - add folder structure 

program  bim2graphkin0
local    function bim2graphkin0
local    `function'_version v2
syntax , bim(string asis) [d(real 0.3540) f(real 0.1770) s(real 0.0884) t(real 0.0442)]

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
qui { // module 3 - calculate kinship matrix	
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3: calculate kinship matrix"
	noi di as text"#########################################################################"
	local already_created ${init_root}/bim2kin0/``function'_bim_short'/``function'_bim_short'_kin0.dta
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		}
	else {
		noi di as text"# > .................. create "as result "``function'_bim_short'_kin0.dta"
		noi bim2kin0, bim(`bim')
		noi checkfile, file(`already_created')
		}
	use `already_created', clear
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 4 - plotting kinship matrix
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 4: plotting kinship data"
	noi di as text"#########################################################################"	
	count
	if `r(N)' > 0 {
		local format msiz(medlarge) msymbol(O) mfc(red) mlc(black) mlabsize(small) mlw(vvthin)
		local xlabel 0(0.1).5
		qui { 	
			tw scatter ibs kin, `format'       ///
			title("data: ``function'_bim_short'") ///
			ytitle("Proportion of Zero IBS") ///
			xlabel("`xlabel'")          ///
			xtitle("Estimated Kinship Coefficient") ///
			xline(`d', lpattern(dash) lwidth(vthin) lcolor(red))  ///
			xline(`f', lpattern(dash) lwidth(vthin) lcolor(red))  ///
			xline(`s', lpattern(dash) lwidth(vthin) lcolor(red)) ///
			xline(`t', lpattern(dash) lwidth(vthin) lcolor(red)) ///
			nodraw saving(``function'_bim_short'_`function'-scatter.gph, replace)
			noi di as text"# > . saving graph as "as result"``function'_bim_short'_`function'-scatter.gph"
			}
		gen rel = ""
		replace rel = "3" if kinship > `t'
		replace rel = "2" if kinship > `s'
		replace rel = "1" if kinship > `f'
		replace rel = "d" if kinship > `d'
		replace rel =  "" if kinship == .
		foreach rel in 3 2 1 d { 
			count if rel == "`rel'"
			local `function'_N`rel' `r(N)'
			}
		qui { 
			tw hist kinship , width(0.005) freq                     ///
			xline(`t', lpattern(dash) lwidth(vthin) lcolor(red)) ///
			xline(`s', lpattern(dash) lwidth(vthin) lcolor(red)) ///
			xline(`f', lpattern(dash) lwidth(vthin) lcolor(red)) ///
			xline(`d', lpattern(dash) lwidth(vthin) lcolor(red)) ///
			xlabel("`xlabel'") legend(off)                             ///
			title("data: ``function'_bim_short'") ///
			caption("Twin/Duplicate Pairs; N = ``function'_Nd'"           ///
				   "1st Degree Relative Pairs ; N = ``function'_N1'"     ///
				   "2nd Degree Relative Pairs ; N = ``function'_N2'"     ///
				   "3rd Degree Relative Pairs ; N = ``function'_N3'")    ///
			nodraw saving(``function'_bim_short'_`function'-hist.gph, replace)
			noi di as text"# > . saving graph as "as result"``function'_bim_short'_`function'-hist.gph"
			}
		noi di as text"# > . saving graph as "as result"``function'_bim_short'_`function'-hist.gph"
		noi di as text"# > ...... duplicates "as result "``function'_Nd'"
		noi di as text"# > .... first degree "as result "``function'_N1'"
		noi di as text"# > ... second degree "as result "``function'_N2'"
		noi di as text"# > .... third degree "as result "``function'_N3'"
		noi di as text"# > .......... export "as result "``function'_bim_short'_`function'.relpairs"
		outsheet if rel != "" using ``function'_bim_short'_`function'.relpairs, noq replace 
		}
	else {
		foreach i in scatter hist {
			noi di as text"# > ................. nothing to plot - generating blank plot"
			tw scatteri 1 1, msymbol(i) ylab("") xlab("") ytitle("") xtitle("") yscale(off) xscale(off) plotregion(lpattern(blank)) nodraw saving(``function'_bim_short'_`function'-`i'.gph, replace)
			noi di as text"# > . saving graph as "as result"``function'_bim_short'_`function'-`i'.gph"
			}
		}
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
	file write myfile `"# > .... indiv. in bim ${bim2count_ind}"' _n
	file write myfile `"# > ....... duplicates ``function'_Nd'"' _n
	file write myfile `"# > ..... first degree ``function'_N1'"' _n
	file write myfile `"# > .... second degree ``function'_N2'"' _n
	file write myfile `"# > ..... third degree ``function'_N3'"' _n
	file write myfile `"# > ...... saved graph ``function'_bim_short'_`function'-hist.gph"' _n
	file write myfile `"# > ...... saved graph ``function'_bim_short'_`function'-scatter.gph"' _n
	file write myfile `"# > ...... saved graph ``function'_bim_short'_`function'.relpairs"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > overview of `function'"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > 1. create/use and _kin0.dta files from bim "' _n
	file write myfile `"# > 2. define relatedness according to thresholds"' _n
	file write myfile `"# > 3. plot ibs x kin"' _n
	file write myfile `"# > 4. plot kin"' _n
	file write myfile `"# > 5. report related individuals"' _n
	file write myfile `"#########################################################################"' _n
	file close myfile	
	}
qui { // module 6 - clean folder
	files2dta, dir(`c(pwd)')
	gen keep = .
	replace keep = 1 if file == "``function'_bim_short'_`function'-scatter.gph"
	replace keep = 1 if file == "``function'_bim_short'_`function'-hist.gph"
	replace keep = 1 if file == "``function'_bim_short'_`function'.relpairs"
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
	
