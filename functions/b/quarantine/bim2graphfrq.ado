*! 1.0.2 Richard Anney 17apr2019
* clean up-front functions
*! 1.0.1 Richard Anney 24jan2019
* adaptation from graphplinkfrq
* - update to start from bim file
* - add folder structure 

program  bim2graphfrq
local    function bim2graphfrq
local    `function'_version v2
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
	foreach file in bim bed fam {
		noi checkfile, file(`bim'.`file')
		}
	noi checkfile, file(${plink})
	checktabbed
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 3 - calculate frq by chromosome	
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3: calculate allele frequency"
	noi di as text"#########################################################################"
	local already_created ${init_root}/bim2frq/``function'_bim_short'/``function'_bim_short'_frq.dta
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		}
	else {
		noi di as text"# > .................. create "as result "``function'_bim_short'_frq.dta"
		noi bim2frq, bim(`bim')
		noi checkfile, file(`already_created')
		}
	use `already_created', clear
	for var maf : lab var X "minor allele frequency"
	count
	local `function'_snpN `r(N)'
	noi di as text"# > .... SNPs in file "as result"``function'_snpN'"
	count if maf < .01
	local `function'_lowN `r(N)'
	noi di as text"# > ...... SNPs < .1% "as result"``function'_lowN'"
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 4 - plot frq by chromosome	
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 4: plot allele frequency"
	noi di as text"#########################################################################"
	sum maf
	if `r(min)' != `r(max)' {
		tw hist maf,  width(0.002) start(0) percent ///
		xlabel(0(.05)0.5) ///
		xline(.01 , lpattern(dash) lwidth(vthin) lcolor(red) ) ///
		legend(off) title("data: ``function'_bim_short'") ///
		caption("SNPs in dataset; N = ``function'_snpN'" ///
		        "SNPs with maf < .01 ; N = ``function'_lowN'") ///
		nodraw saving(``function'_bim_short'_`function'.gph, replace)
		noi di as text"# > . saving graph as "as result"``function'_bim_short'_`function'.gph"
		}
	else {
		noi di as text"# > ................. nothing to plot - generating blank plot"
		tw scatteri 1 1, msymbol(i) ylab("") xlab("") ytitle("") xtitle("") yscale(off) xscale(off) plotregion(lpattern(blank)) nodraw saving(``function'_bim_short'_`function'.gph, replace)
		noi di as text"# > . saving graph as "as result"``function'_bim_short'_`function'.gph"
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
	file write myfile `"# > ...... SNPs in bim ``function'_snpN'"' _n
	file write myfile `"# > ....... SNPs < .1% ``function'_lowN'"' _n
	file write myfile `"# > ...... saved graph ``function'_bim_short'_`function'.gph"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > overview of `function'"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > 1. create/use and _frq.dta files from bim "' _n
	file write myfile `"# > 2. plot frequencies"' _n
	file write myfile `"#########################################################################"' _n
	file close myfile	
	}
qui { // module 5 - clean folder
	files2dta, dir(`c(pwd)')
	gen keep = .
	replace keep = 1 if file == "``function'_bim_short'_`function'.gph"
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
	
