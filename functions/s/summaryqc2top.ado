*! 1.0.1 Richard Anney 28jan2019
* - convert summaryqc to top hits based on ld-clumping routine
* - assume snp chr bp p
* - update to include summaryqc
*! 1.0.2 Richard Anney 289jan2020
* - drop GENES - replace with summaryqc2top_annot

program  summaryqc2top
local    function summaryqc2top
local    `function'_version v2

syntax , summaryqc(string asis) ldref(string asis) [minp(real 1e-5)]

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
	ref_path_to_short, ref(`summaryqc')
	local `function'_summaryqc_short ${ref_short}
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
	noi checkfile, file(${plink2})
	checktabbed
	noi di as text"# > ......... locating files to process"
	noi checkfile, file(`summaryqc'-summaryqc.dta)
	foreach file in bim bed fam {
		noi checkfile, file(`ldref'.`file')
		}
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 3 - perform data clumping
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3: perform summaryqc clumping"
	noi di as text"#########################################################################"
	noi di as text"# >  loading snps from "as result"``function'_summaryqc_short'"
	noi di as text"# > ..... limit to P < "as result"`: display %10.4e `minp''"
	local already_created ``function'_summaryqc_short'.top
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		}
	else {	
		use snp p if p < = .05 using `summaryqc'-summaryqc.dta, clear
		keep snp p
		renvars, upper
		outsheet using  ``function'_summaryqc_short'.top, noq replace
		noi checkfile, file(`already_created')
		}
	noi di as text"# > ...... clump using "as result"``function'_ldref_short'"
	local already_created ``function'_summaryqc_short'.clumped
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		}
	else {		
		!${plink} --bfile  `ldref' --clump ``function'_summaryqc_short'.top --clump-allow-overlap --clump-p1 `minp' --clump-p2 0.05 --clump-kb 1000 --clump-r2 0.1 --out ``function'_summaryqc_short'
		!${tabbed} ``function'_summaryqc_short'.clumped
		noi checkfile, file(`already_created')
		}
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 4 - define LD Range
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 4: define LD Range"
	noi di as text"#########################################################################"
	local already_created ${init_root}/bim2dta/``function'_ldref_short'/``function'_ldref_short'_bim.dta
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		}
	else {	
		bim2dta, bim(`ldref')
		noi checkfile, file(`already_created')
		}
	noi di as text"# > .............. add "as result"LD Ranges"
	local already_created ``function'_summaryqc_short'-`function'.dta
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		}
	else {
		import delim using ``function'_summaryqc_short'.clumped.tabbed,clear varnames(1) case(upper)
		keep SNP SP2
		replace SP2 = subinstr(SP2,"(1)","",.)
		split SP2,p(",")
		drop SP2
		rename SNP INDEX
		reshape long SP2,i(INDEX) 
		drop if SP2 == ""
		drop _j
		rename SP2 snp
		merge m:1 snp using ${init_root}/bim2dta/``function'_ldref_short'/``function'_ldref_short'_bim.dta
		keep if _m == 3
		keep INDEX snp chr bp
		destring bp, replace
		sort INDEX bp
		egen FROM = seq(), by(INDEX)
		gsort INDEX -bp
		egen TO = seq(), by(INDEX)
		gen keep = ""
		replace keep = "Y" if FROM == 1
		replace keep = "Y" if   TO == 1
		keep if keep == "Y"
		sort INDEX FROM
		egen x = seq(),by(INDEX)
		tostring bp, replace
		replace bp = "chr" + chr + ":" + bp if x == 1
		keep INDEX bp x
		reshape wide bp, i(INDEX) j(x)
		gen LD_RANGE = bp1 + ".." + bp2
		keep INDEX LD_RANGE
		rename INDEX snp	
		noi di as text"# > ............ merge "as result"``function'_summaryqc_short'-summaryqc"
		merge 1:1 snp using `summaryqc'-summaryqc.dta
		keep if _m == 3
		renvars, upper
		for var CHR BP OR SE P BETA Z: destring X, replace force
		compress
		sort SNP 
		tostring CHR, replace
		tostring BP,replace
		gen POSITION = "chr" + CHR + ":" + BP
		drop CHR BP
		gen lnOR = ln(OR)
		gen lnLB = lnOR - (1.96*SE)
		gen lnUB = lnOR + (1.96*SE)
		gen LB = exp(lnLB) 
		gen UB = exp(lnUB) 
		gen str4 X= string(OR,"%05.4f")
		gen str4 Y= string(LB,"%05.4f")
		gen str4 Z_= string(UB,"%05.4f") 
		drop lnOR lnLB lnUB LB UB OR 
		gen OR = X + " (" + Y + "-" + Z_ + ")"
		drop X Y Z_
		keep  SNP POSITION A1 A2 BETA SE Z OR P LD_RANGE
		order SNP POSITION A1 A2 BETA SE Z OR P LD_RANGE  
		lab var SNP 			"Marker Name"
		lab var POSITION 	"Chromosome:Physical Location"
		lab var OR        "Odds Ratio ((5% Confidence Intervals)"
		lab var LD_RANGE  "Clumped Region"
		sort P
		renvars, lower
		compress
		save ``function'_summaryqc_short'-`function'.dta, replace
		noi checkfile, file(`already_created')
		}
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 5 - write log file
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
	file write myfile `"# > ....... reference ``function'_ldref_short'"' _n
	file write myfile `"# > .. min-P reported `minp'"' _n
	file write myfile `"# > ...... saved data ``function'_summaryqc_short'-`function'.dta"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > overview of `function'"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > 1. import top N SNPs"' _n
	file write myfile `"# > 2. clump versus reference"' _n
	file write myfile `"# > 3. prepare data for publication"' _n
	file write myfile `"#########################################################################"' _n
	file close myfile	
	}
qui { // module 6 - clean folder
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

   
