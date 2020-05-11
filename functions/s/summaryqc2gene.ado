*! 1.0.1 Richard Anney 22jan2020
*! 1.0.2 Richard Anney 5may2020

program  summaryqc2gene
local    function summaryqc2gene
local    `function'_version v2

syntax , summaryqc(string asis) ldref(string asis) generef(string asis) [minp(real 1e-5)]

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
		local `function'_glist_short Homo_sapiens.GRCh37.87.gtf
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
	noi checkfile, file(`summaryqc'-summaryqc.dta)
	noi checkfile, file(`generef'/Homo_sapiens.GRCh37.87.gtf_exon.dta)
	local already_created `generef'/Homo_sapiens.GRCh37.87.gtf.hg-list
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		}
	else {	
		use `generef'/Homo_sapiens.GRCh37.87.gtf_exon.dta, clear
		keep chr start end symbol 
		for var chr start end: tostring X, replace
		duplicates drop
		outsheet chr start end symbol using `generef'/Homo_sapiens.GRCh37.87.gtf.hg-list, non noq delim(" ") replace
		noi checkfile, file(`already_created')
		}
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 3 - identify genes overlapping ranges
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3: map by regions"
	noi di as text"#########################################################################"
	qui { // load summaryqc2top
		noi di as text"# > ..... loading from "as result"``function'_summaryqc_short'-summaryqc2top.dta"
		local already_created ${init_root}/summaryqc2top/``function'_summaryqc_short'/``function'_summaryqc_short'-summaryqc2top.dta
		capture confirm file `already_created'
		if !_rc {
			noi checkfile, file(`already_created')
			}
		else {	
			noi summaryqc2top, summaryqc(`summaryqc') ldref(`ldref')
			noi checkfile, file(`already_created')
			}
		}
	qui { // create summaryqc2gene.dta
		local already_created ${init_root}/`function'/``function'_summaryqc_short'/``function'_summaryqc_short'-`function'.dta
		capture confirm file `already_created'
		if !_rc {
			noi checkfile, file(`already_created')
			}
		else {
			use ${init_root}/summaryqc2top/``function'_summaryqc_short'/``function'_summaryqc_short'-summaryqc2top.dta, clear
			keep ld_range
			rename ld intersect_genes_range 
			split intersect_genes_range,p("chr"":""..")
			for var intersect_genes_range3 intersect_genes_range4: destring X, replace
			replace intersect_genes_range4 = intersect_genes_range3 + 1 if intersect_genes_range4 == .
			for var intersect_genes_range3 intersect_genes_range4: tostring X, replace
			gen new_range = "chr" + intersect_genes_range2 + ":" + intersect_genes_range3 + ".." + intersect_genes_range4
			replace intersect_genes_range = new_range if new_range != intersect_genes_range
			keep intersect_genes_range
			save ${init_root}/`function'/``function'_summaryqc_short'/``function'_summaryqc_short'-range-list.dta, replace
			intersect_genes, range_list(${init_root}/`function'/``function'_summaryqc_short'/``function'_summaryqc_short'-range-list.dta) generef(`generef')
			use ${init_root}/intersect_genes/Homo_sapiens.GRCh37.87.gtf_exon/``function'_summaryqc_short'-range-list-intersect_genes.dta, clear
			keep intersect symbol biotype
			duplicates drop
			sort biotype symbol
			egen x = seq(),by(inter bio)
			reshape wide symbol, i(intersect biotype) j(x)
			gen symbol999 = ""
			gen symbol = ""
			aorder
			for var symbol1 - symbol999: replace symbol = symbol + " " + X
			keep intersect symbol biotype
			rename intersect ld_range
			merge m:1 ld_range using ${init_root}/summaryqc2top/``function'_summaryqc_short'/``function'_summaryqc_short'-summaryqc2top.dta
			drop if _m == 1
			drop _m
			lab var symbol "Gene Symbol"
			order snp-p
			sort p bio
			local fmt: format symbol
			local fmt: subinstr local fmt "%" "%-"
			format symbol `fmt'
			replace symbol = "x!x" + symbol
			replace symbol = subinstr(symbol,"x!x ","",.)
			keep if p < `minp'
			compress
			save ``function'_summaryqc_short'-`function'.dta, replace
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
	file write myfile `"# > .. gene reference ``function'_glist_short'"' _n
	file write myfile `"# > .. min-P reported `minp'"' _n
	file write myfile `"# > ...... saved data ``function'_summaryqc_short'-`function'.dta"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > overview of `function'"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > 1. map symbols to ranges"' _n
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

   
