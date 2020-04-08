*! 1.0.2 Richard Anney 21jan2019

program  bim2array
local    function bim2array
local    `function'_version v2
syntax , bim(string asis) dir(string asis) 

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
qui { // module 2 - check files
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 2: check files"
	noi di as text"#########################################################################"
	noi checkfile, file(`bim'.bim)
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 3 - import SNPs
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3: import SNPs"
	noi di as text"#########################################################################"
	qui { // count snps in array
		import delim using `bim'.bim, clear
		bim2count, bim(`bim')
		noi di as text"# > ..... SNPs in file "as result "${bim2count_snp}"
		}
	qui { // import snp list
		noi di as text"# > ........ importing "as result"``function'_bim_short'.bim"
		capture confirm variable v2
		if !_rc { 
		}
		else {
		rename v1 v
		split v,p(" ")
		drop v
		}
		keep v2
		tostring v2, replace
		rename v2 snp
		save ``function'_bim_short'-`function'.dta, replace
		}
	qui { // create `bim'.array
		clear
		set obs 1
		gen a = "array jaccard-index"
		outsheet a using ``function'_bim_short'.array, non noq replace
		}
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 4 - populate ``function'_bim_short'.array
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 4: populate ``function'_bim_short'.array"
	noi di as text"#########################################################################"		
	qui { // define arrays by number
		files2dta, dir(`dir')
		erase _files2dta.dta
		drop if files == "_files2dta.dta"
		drop if files == "``function'_bim_short'-`function'.dta"
		split files, p(".dt")
		keep if files2 == "a"
		keep files1
		duplicates drop
		gen n = _n
		tostring n, replace
		gen a = "global `function'" + n + " " + files1
		rename a script
		keep script
		do_temp_script
		count
		local `function'_Nvs `r(N)'

		noi di as text"# > .................. this can take time - comparing to "as result"``function'_Nvs'" as text " arrays/ datasets"
		}
	qui { // populate bim to array
		foreach num of num 1 / `r(N)' { 
			use ``function'_bim_short'-`function'.dta, replace 
			duplicates drop
			merge 1:1 snp using `dir'/${`function'`num'}.dta
			gen array = "${`function'`num'}"
			sum _m 
			gen all = `r(N)'
			count if _m == 3
			gen ab = `r(N)'
			gen jaccard = ab/(all)
			keep array jaccard
			sum jaccard
			qui di as text"# > .. jaccard index = " as result trim("`: display %5.4f r(min)'") as text " for array " as result "${`function'`num'}" 
			filei + "${`function'`num'} `r(min)'" ``function'_bim_short'.array
			}
		erase ``function'_bim_short'-`function'.dta
		}
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 5 - determine most likely array
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 5: determine most likely array"
	noi di as text"#########################################################################"
	import delim using "``function'_bim_short'.array", clear delim(" ") varnames(1) case(preserve)
	duplicates drop
	gsort -j
	gen a1 = "global `function'_array " + array in 1
	gen str6 jaccard2 = string(jaccard, "%5.4f") 
	gen a2 = "global `function'_jaccard " + jaccard2 in 1 
	keep a1 a2
	gen n = 1
	keep in 1
	reshape long a, j(x) i(n)
	keep a
	outsheet a using `function'_array_tmp.do, non noq replace
	do `function'_array_tmp.do
	erase `function'_array_tmp.do
	noi di as text"# > ............. from "as result "``function'_bim_short'.array" 
	noi di as text"# >  most likely array "as result "${`function'_array}" 
	noi di as text"# > .... jaccard index "as result "${`function'_jaccard}" 
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 6 - save files
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 6: save files and clean up temporary folder"
	import delim using "``function'_bim_short'.array", clear delim(" ") varnames(1) case(preserve)
	duplicates drop
	gsort -j	
	outsheet using ``function'_bim_short'.array, noq replace
	outsheet using `bim'.array, noq replace
	}
qui { // module 7 - write log file
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
	file write myfile `"# > .... SNPs in array ${bim2count_snp}"' _n
	file write myfile `"# > .......... map vs. ``function'_Nvs' arrays"' _n
	file write myfile `"# > ...... most likely ${`function'_array}"' _n
	file write myfile `"# > .... jaccard index ${`function'_jaccard}"' _n
	file write myfile `"# > ....... saved data ``function'_bim_short'.array"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > overview of `function'"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > 1. import markers"' _n
	file write myfile `"# > 2. merge against known SNPs"' _n
	file write myfile `"# > 3. calculate the jaccard overlap of named markers"' _n
	file write myfile `"# > 4. define array as that with highest jaccard coefficient"' _n
	file write myfile `"#########################################################################"' _n
	file close myfile	
	}
qui { // module 8 - clean folder
	files2dta, dir(`c(pwd)')
	gen keep = .
	replace keep = 1 if file == "``function'_bim_short'.array"
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

