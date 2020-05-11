*! 1.0.2 Richard Anney 5may2020

program  summaryqc2cadd
local    function summaryqc2cadd
local    `function'_version v2

syntax , summaryqc(string asis) ldref(string asis) caddref(string asis) [minp(real 1e-5) ]

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
	qui { // define short names
		ref_path_to_short, ref(`summaryqc')
		local `function'_summaryqc_short ${ref_short}
		ref_path_to_short, ref(`ldref')
		local `function'_ldref_short ${ref_short}
		ref_path_to_short, ref(`caddref')
		local `function'_clist_short ${ref_short}	
		}
	cd ${init_root}	
	!mkdir `function'
	cd `function'
	noi di as text"# > ..... project name "as result"``function'_summaryqc_short'"
	!mkdir ``function'_summaryqc_short'
	cd ``function'_summaryqc_short'
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 2 - define files
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 2: define files"
	noi di as text"#########################################################################"
	noi checkfile, file(${plink})
	checktabbed
	noi di as text"# > ......... locating files to process"
	noi checkfile, file(`summaryqc'-summaryqc.dta)
	foreach file in bim bed fam {
		noi checkfile, file(`ldref'.`file')
		}
	noi checkfile, file(`caddref')
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 3 - define ld pairs
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3: define ld pairs"
	noi di as text"#########################################################################"
	noi di as text"# > ..... loading from "as result"``function'_summaryqc_short'-summaryqc2top.dta"
	noi di as text"# > ..... limit to P < "as result"`: display %10.4e `minp''"
	local already_created ${init_root}/summaryqc2top/``function'_summaryqc_short'/``function'_summaryqc_short'-summaryqc2top.dta
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		}
	else {	
		noi summaryqc2top, summaryqc(`summaryqc') ldref(`ldref') 
		noi checkfile, file(`already_created')
		}
	local already_created ${init_root}/snp2ldpairs/``function'_summaryqc_short'-summaryqc2top/``function'_summaryqc_short'-summaryqc2top-snp2ldpairs.dta
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		}
	else {	
		use ${init_root}/summaryqc2top/``function'_summaryqc_short'/``function'_summaryqc_short'-summaryqc2top.dta, clear
		keep if p < `minp'
		outsheet snp using ``function'_summaryqc_short'-summaryqc2top.snplist, non noq replace
		snp2ldpairs , ldref(${ldref}) snplist(${init_root}/`function'/``function'_summaryqc_short'/``function'_summaryqc_short'-summaryqc2top)		
		noi checkfile, file(`already_created')
		}
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 4 - map to cadd
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 4: map to cadd"
	noi di as text"#########################################################################"
	qui { // define cadd mapped to proxy
		use `already_created', clear
		keep snp_proxy
		rename snp_proxy snp
		duplicates drop
		noi di as text"# > .. merge with data "as result"``function'_clist_short'"
		merge 1:m snp using `caddref'
		drop if _m == 2
		count if _m == 3
		if `r(N)' == 0 {
			noi di as text"# > .................. "as result"no CADD annotations observed"
			use ${init_root}/summaryqc2top/``function'_summaryqc_short'/``function'_summaryqc_short'-summaryqc2top.dta, clear
			gen cadd = "none"
			gen cadd_r2 = "NA"
			lab var cadd "cadd annotation (max)"
			lab var cadd_r2 "Linkage Disequilibrium between index and reported annotation"
			keep if p < `minp'
			compress
			save ``function'_summaryqc_short'-`function'.dta,replace	
			}
		else {	
			drop _m
			rename snp snp_proxy
			duplicates drop
			qui { // define most significant cadd per snp
				for var consscore phred: destring X, replace
				for var polyp sift: replace X = "" if X == "NA"
				gsort -conss phred
				egen x = seq(),by(snp)
				}			
			tostring conss, replace
			gen a = string(phred,"%03.1f")
			gen cadd = ""
			replace cadd = snp_proxy + "-" + a1 + " " + consequence + " " + consscore + " (cadd_phred = " + a + ")"
			replace cadd = snp_proxy + "-" + a1 + " " + consequence + " " + consscore + " (cadd_phred = " + a + "; " + polyp + "/" + sift + ")" if sift != ""
			keep if x == 1
			keep snp_proxy cadd
			merge m:m snp_proxy using `already_created'
			keep if _m == 3
			keep snp snp_proxy cadd r2
			merge m:1 snp using ${init_root}/summaryqc2top/``function'_summaryqc_short'/``function'_summaryqc_short'-summaryqc2top.dta
			keep  snp position a1 a2 beta se z or p cadd r2
			order snp position a1 a2 beta se z or p cadd r2
			split_cadd
			gsort -cadd_consscore -cadd_phred
			egen x = seq(), by(snp )
			keep if x == 1
			sort p
			keep  snp position a1 a2 beta se z or p cadd  r2
			lab var cadd "cadd annotation (max)"
			replace cadd = "none" if cadd == ""
			gen cadd_r2 = string(r2,"%04.2f")
			replace cadd_r2 = "NA" if cadd_r2 == ""
			lab var cadd_r2 "Linkage Disequilibrium between index and reported annotation"
			drop r2
			keep if p < `minp'
			compress
			save ``function'_summaryqc_short'-`function'.dta,replace
			}
		}
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 4 - write log file
	file open myfile using "``function'_summaryqc_short'-`function'.log", write replace
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
	file write myfile `"# > ....... summaryqc ``function'_summaryqc_short'"' _n
	file write myfile `"# > .... LD reference ``function'_ldref_short'"' _n
	file write myfile `"# > .. CADD reference ``function'_clist_short'"' _n
	file write myfile `"# > .. min-P reported `: display %10.4e `minp''"' _n
	file write myfile `"# > ...... saved data ``function'_summaryqc_short'-`function'.dta"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > overview of `function'"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > 1. identify top N SNPs"' _n
	file write myfile `"# > 2. extend LD to capture proxy_SNPs"' _n
	file write myfile `"# > 3. map eQTL pairs (abridged to *.dta)"' _n
	file write myfile `"# > 4. map eQTL pairs (detailed to detailed.dta"' _n
	file write myfile `"#########################################################################"' _n
	file close myfile	
	}
qui { // module 5 - clean folder
	files2dta, dir(`c(pwd)')
	gen keep = .
	replace keep = 1 if file == "``function'_summaryqc_short'-`function'.log"
	replace keep = 1 if file == "``function'_summaryqc_short'-`function'.dta"
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

   
