program  summaryqc2metal
local    function summaryqc2metal
local    `function'_version v1
syntax , summaryqc1(string asis) summaryqc2(string asis) ref(string asis)

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
	ref_path_to_short, ref(`ref')
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
	noi checkfile, file(${metal})
	checktabbed
	noi di as text"#########################################################################"
	noi di as text""
	}
qui { // module 3 - run meta
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3a: create chr files"
	noi di as text"#########################################################################"
	noi di as text"# > .......... process "as result"``function'_summaryqc1'"
	foreach j of num 1/22 {
		local already_created ${init_root}/`function'/``function'_project'/``function'_summaryqc1'.chr`j'
		capture confirm file `already_created'
		if !_rc {
			checkfile, file(`already_created')
			}
		else {
			use chr snp a1 a2 beta se n using `summaryqc1'-summaryqc.dta, clear
			recodegenotype
			drop if recode == "W"
			drop if recode == "S"
			renvars, upper
			foreach chr of num 1/22 {
				outsheet SNP A1 A2 BETA SE N if CHR == "`chr'"  using ``function'_summaryqc1'.chr`chr', noq replace
				checkfile, file(${init_root}/`function'/``function'_project'/``function'_summaryqc1'.chr`chr')
				}
			}
		
		}
	noi di as text"# > .......... process "as result"``function'_summaryqc2'"
	foreach j of num 1/22 {
		local already_created ${init_root}/`function'/``function'_project'/``function'_summaryqc2'.chr`j'
		capture confirm file `already_created'
		if !_rc {
			checkfile, file(`already_created')
			}
		else {		
			use chr snp a1 a2 beta se n using `summaryqc2'-summaryqc.dta, clear
			recodegenotype
			drop if recode == "W"
			drop if recode == "S"
			renvars, upper
			foreach chr of num 1/22 {
				outsheet SNP A1 A2 BETA SE N if CHR == "`chr'"  using ``function'_summaryqc2'.chr`chr', noq replace
				checkfile, file(${init_root}/`function'/``function'_project'/``function'_summaryqc2'.chr`chr')
				}
			}
		}
	noi di as text"#########################################################################"
	noi di as text""
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3b: create meta files"
	noi di as text"#########################################################################"
	foreach chr of num 1/22 {
		local already_created ${init_root}/`function'/``function'_project'/``function'_project'-chr`chr'.meta
		capture confirm file `already_created'
		if !_rc {
			checkfile, file(`already_created')
			}
		else {
			file open myfile using "${init_root}/`function'/``function'_project'/``function'_project'-chr`chr'.meta", write replace
			file write myfile "SCHEME STDERR" _n
			file write myfile "VERBOSE OFF" _n
			file write myfile "MARKER SNP" _n
			file write myfile "ALLELE A1 A2" _n
			file write myfile "EFFECT BETA" _n
			file write myfile "STDERR SE" _n
			file write myfile "CUSTOMVARIABLE TotalSampleSize" _n
			file write myfile "LABEL TotalSampleSize as N" _n
			file write myfile `"PROCESS ${init_root}/`function'/``function'_project'/``function'_summaryqc1'.chr`chr'"' _n
			file write myfile `"PROCESS ${init_root}/`function'/``function'_project'/``function'_summaryqc2'.chr`chr'"' _n
			file write myfile `"OUTFILE ${init_root}/`function'/``function'_project'/``function'_project'-chr`chr'- .tbl"' _n
			file write myfile "ANALYZE" _n
			file close myfile	
			}
		}
	noi di as text"# > .................. "as result"completed"
	noi di as text"#########################################################################"
	noi di as text""
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3c: run meta analysis"
	noi di as text"#########################################################################"
	foreach chr of num 1/22 {
		local already_created ${init_root}/`function'/``function'_project'/``function'_project'-chr`chr'-1.tbl
		capture confirm file `already_created'
		if !_rc {
			checkfile, file(`already_created')
			}
		else {
			!${metal} ${init_root}/`function'/``function'_project'/``function'_project'-chr`chr'.meta     > ${init_root}/`function'/``function'_project'/``function'_project'-chr`chr'.log
			checkfile, file(`already_created')
			}
		}
	noi di as text"# > .................. "as result"completed"
	noi di as text"#########################################################################"
	noi di as text""
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3d: define chr bp in study"
	noi di as text"#########################################################################"
	local already_created ${init_root}/`function'/``function'_project'/``function'_project'-`function'-chrbp.dta
	capture confirm file `already_created'
	if !_rc {
		checkfile, file(`already_created')
		}
	else {
		use chr bp snp using `summaryqc1'-summaryqc.dta, clear
		merge 1:1 snp  using `summaryqc2'-summaryqc.dta
		keep if _m == 3
		drop _m
		keep snp chr bp
		compress
		save ``function'_project'-`function'-chrbp.dta, replace
		}
	noi di as text"# > .................. "as result"completed"
	noi di as text"#########################################################################"
	noi di as text""
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3e: merge metal files"
	noi di as text"#########################################################################"
	local already_created ${init_root}/summaryqc/``function'_project'/``function'_project'-summaryqc.dta
	capture confirm file `already_created'
	if !_rc {
		checkfile, file(`already_created')
		}
	else {
		foreach chr of num 1/22 {
			import delim using ${init_root}/`function'/``function'_project'/``function'_project'-chr`chr'-1.tbl, clear
			rename (mark pv) (snp p_meta)
			save ${init_root}/`function'/``function'_project'/``function'_project'-chr`chr'.dta, replace
			}
		foreach chr of num 1/21 {
			append using ${init_root}/`function'/``function'_project'/``function'_project'-chr`chr'.dta
			}
		merge 1:1 snp using ``function'_project'-`function'-chrbp.dta
		keep if _m == 3
		rename (allele1 allele2 effect stderr p_meta totalsamplesize) (a1 a2 beta se p n)
		for var a1 a2: replace X = strupper(X)
		drop _m
		sum p if p ! = 0
		replace p = `r(min)' if p == 0
		gen z    = beta / se
		gen or   = exp(beta)
		gen l95  = or - (1.96 * se)
		gen u95  = or + (1.96 * se)
		noi summaryqc, ref(`ref') input(``function'_project'-`function') project(``function'_project') png(yes)
		}
	noi di as text"# > .................. "as result"completed"
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
	file write myfile `"# > ...... saved data ``function'_project'-`function'-summaryqc.dta"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > overview of `function'"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > 1. create chr files from summaryqc"' _n
	file write myfile `"# > 2. create meta files from summaryqc"' _n
	file write myfile `"# > 3. perform ivw meta analysis of 2 datasets using metal"' _n
	file write myfile `"# > 4. combine tbl output files "' _n
	file write myfile `"# > 5. apply summaryqc to data"' _n
	file write myfile `"#########################################################################"' _n
	file close myfile	
	}
qui { // module z - clean folder
	files2dta, dir(`c(pwd)')
	gen keep = .
	replace keep = 1 if file == "``function'_project'-`function'.log"
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
