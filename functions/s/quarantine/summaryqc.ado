*! 1.0.1 Richard Anney 19dec2018

program  summaryqc
local    function summaryqc
local    `function'_version v2
syntax ,  ref(string asis) input(string asis) project(string asis) [png(string asis)]

noi di as text""
noi di as text"#########################################################################"
noi di as text"# > .................. `function'"
noi di as text"# > ........... author "as result"Richard Anney"
noi di as text"# > ... report bugs to "as result"anneyr@cardiff.ac.uk"
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
	ref_path_to_short, ref(`ref')
	local `function'_ref_short ${ref_short}
	ref_path_to_short, ref(`input')
	local `function'_input_short ${ref_short}
	cd ${init_root}	
	!mkdir `function'
	cd `function'
	noi di as text"# > ..... project name "as result"`project'"
	noi di as text"# > ..... derived from "as result"``function'_input_short'"
	!mkdir `project'
	cd `project'
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 2 - define files
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 2: define files"
	noi di as text"#########################################################################"
	noi di as text"# > ......... locating files to process"
	noi checkfile, file(`ref')
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 3 - process file
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3: process ``function'_input_short'"
	noi di as text"#########################################################################"
	qui { // SNPs in
		count
		local summaryqc_Nin `r(N)'
		noi di as text"# > ..... SNPs in file "as result "`summaryqc_Nin'"
		noi di as text"#########################################################################"
		}
	qui { // P out of bounds
		noi di as text"# > .................. drop SNPs with P out-of bounds"
		drop if p > 1
		drop if p < 0
		drop if p == .
		count
		local count `r(N)'
		gen count =  `count' - `summaryqc_Nin'
		sum count
		local summaryqc_oobSNP `r(max)'
		drop count
		count
		local running `r(N)'
		noi di as text"# > ..... SNPs dropped "as result "`summaryqc_oobSNP'"
		noi di as text"# > ..... SNPs in file "as result "`running'"
		noi di as text"#########################################################################"
		}
	qui { // duplicates
		noi di as text"# > .................. drop duplicates"
		duplicates drop 
		duplicates drop snp, force
		count
		local count `r(N)'
		gen count = `count' - `running'
		sum count
		local summaryqc_dupSNP `r(max)'
		drop count
		count
		local running `r(N)'
		noi di as text"# > ..... SNPs dropped "as result "`summaryqc_dupSNP'"
		noi di as text"# > ..... SNPs in file "as result "`running'"
		noi di as text"#########################################################################"
		}
	qui { // info
		capture confirm numeric var info
		if !_rc {
			noi di as text"# > .................. drop SNPs with INFO out-of-bounds"
			drop if info == .
			drop if info > 2 & info ! = .
			drop if info < .8
			drop info
			count
			local count `r(N)'
			gen count3 = `count' - `summaryqc_dupSNP'
			sum count3
			local summaryqc_infoSNP `r(max)'
			count
			local running `r(N)'
			drop count
			noi di as text"# > ..... SNPs dropped "as result "`summaryqc_infoSNP'"
			noi di as text"# > ..... SNPs in file "as result "`running'"
			noi di as text"#########################################################################"
			}
		else {
			noi di as text"# > .................. drop SNPs with INFO out-of-bounds"
			count
			local running `r(N)'
			local summaryqc_infoSNP NA
			noi di as text"# > ..... SNPs dropped "as result "0" as text " - info not present"
			noi di as text"# > ..... SNPs in file "as result "`running'"
			noi di as text"#########################################################################"
			}
		}
	qui { // direction
		capture confirm string var direction
		if !_rc {
			noi di as text"# > .................. drop SNPs with missingness within meta-analysis"
			replace direction = subinstr(direction, "-", "",.)
			replace direction = subinstr(direction, "+", "",.)
			gen count4 = length(direction)
			drop if count4 > 2
			drop count
			count
			local count `r(N)'
			gen count5 = `count' - `summaryqc_dupSNP'
			sum count5
			local summaryqc_directionSNP `r(max)'
			count
			local running `r(N)'
			drop count
			noi di as text"# > ..... SNPs dropped "as result "`summaryqc_directionSNP'"
			noi di as text"# > ..... SNPs in file "as result "`running'"
			noi di as text"#########################################################################"
			}
		else {
			noi di as text"# > .................. drop SNPs with missingness within meta-analysis"
			count
			local running `r(N)'
			local summaryqc_directionSNP NA

			noi di as text"# > ..... SNPs dropped "as result "0" as text " - direction not present"
			noi di as text"# > ..... SNPs in file "as result "`running'"
			noi di as text"#########################################################################"
			}		
		}
	}
qui { // module 4 - saving file
	noi di as text""
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 4: saving file"
	noi di as text"#########################################################################"	
	count
	local summaryqc_Nout `r(N)'
	noi di as text"# > ..... SNPs in file "as result "`summaryqc_Nout'"
	order chr bp snp a1 a2 beta se z or l95 u95 p n
	keep  chr bp snp a1 a2 beta se z or l95 u95 p n
	sort  chr bp
	for var chr bp: tostring X, replace	
	compress
	noi di as text"# > ..... save data as "as result"`project'-`function'.dta"
	save `project'-`function'.dta, replace
	noi di as text"#########################################################################"	
	}
qui { // module 4 - plot updated manhattan
	noi di as text""
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 5: plot updated manhattan"
	noi di as text"#########################################################################"		
	graphmanhattan, summaryqc(${init_root}/summaryqc/`project'/`project') png(`png')
	checkfile, file(${init_root}/graphmanhattan/`project'/`project'_graphmanhattan.gph)
	noi di as text""
	}
qui { // module 6 - write log file
	file open myfile using "`project'-`function'.log", write replace
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
	file write myfile`"# > ..... project name `project'"' _n
	file write myfile`"# > ....... input GWAS ``function'_input_short'"' _n
	file write myfile`"# > .... SNP reference ``function'_ref_short'"' _n
	file write myfile`"# > .... SNPs in input `summaryqc_Nin'"' _n 
	file write myfile`"# > .. P out-of-bounds `summaryqc_oobSNP'"' _n 
	file write myfile`"# > ....... duplicates `summaryqc_dupSNP'"' _n
	file write myfile`"# > info out-of-bounds `summaryqc_infoSNP'"' _n
	file write myfile`"# >  direction missing `summaryqc_directionSNP'"' _n
	file write myfile`"# > ... SNPs in output `summaryqc_Nin'"' _n 
	file write myfile`"# > ....... saved data `project'-`function'.dta"' _n
	file write myfile`"#########################################################################"' _n
	file write myfile`"# > overview of `function'"' _n
	file write myfile`"#########################################################################"' _n
	file write myfile`"# > 1. using pre-processed summary files (refid aligned)"' _n
	file write myfile`"# > 2. remove SNPs with P-values out-of-bounds"' _n
	file write myfile`"# > 3. remove duplicate identifiers"' _n
	file write myfile`"# > 4. remove SNPs if info out-of-bounds (<.8)"' _n
	file write myfile`"# > 5. remove SNPs if absent in 2+ studies (metal)"' _n
	file write myfile`"# > 7. plot against reference populations"' _n
	file write myfile`"#########################################################################"' _n
	file close myfile	
	}
qui { // module 8 - clean folder
	files2dta, dir(`c(pwd)')
	gen keep = .
	replace keep = 1 if file == "`project'-`function'.dta"
	replace keep = 1 if file == "`project'-`function'.log"
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

	
