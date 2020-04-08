program  summaryqc2pairwise
local    function summaryqc2pairwise
local    `function'_version v2
syntax , summaryqc1(string asis) summaryqc2(string asis) ldref(string asis) glist_hg19(string asis) [p_top(real 1e-8) p_max(real 1e-5) redo_models(string asis)]

* v2 - remove metal - include meta-analysis calculation in script

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
	noi di as text"#########################################################################"
	noi di as text""
	}
qui { // module 3 - run meta analysis
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3: run meta analysis"
	noi di as text"#########################################################################"
	local already_created ${init_root}/`function'/``function'_project'/``function'_project'-`function'-meta.dta
	capture confirm file `already_created'
	if !_rc {
		checkfile, file(`already_created')
		}
	else {
		noi di as text"# > ........... import "as result"``function'_summaryqc1'"
		use snp a1 a2 beta z se n p using `summaryqc1'-summaryqc.dta, clear
		noi di as text"# > .................. drop ambiguous genotypes"
		drop if a1 == "A" & a2 == "T"
		drop if a1 == "T" & a2 == "A"
		drop if a1 == "C" & a2 == "G"
		drop if a1 == "G" & a2 == "C"
		noi di as text"# > .................. drop monomorphic genotypes"
		drop if a1 == a2
		for var a1 a2 z beta se n p : rename X X_1
		noi di as text"# > ....... merge with "as result"``function'_summaryqc2'"
		merge 1:1 snp using `summaryqc2'-summaryqc.dta
		keep if _m == 3
		noi di as text"# > .................. drop ambiguous genotypes"
		drop if a1 == "A" & a2 == "T"
		drop if a1 == "T" & a2 == "A"
		drop if a1 == "C" & a2 == "G"
		drop if a1 == "G" & a2 == "C"
		noi di as text"# > .................. drop monomorphic genotypes"
		drop if a1 == a2
		for var a1 a2 z beta se n p : rename X X_2
		keep  chr bp snp a1_1 a2_1 beta_1 z_1 se_1 n_1 p_1 a1_2 a2_2 beta_2 z_2 se_2 n_2 p_2
		order chr bp snp a1_1 a2_1 beta_1 z_1 se_1 n_1 p_1 a1_2 a2_2 beta_2 z_2 se_2 n_2 p_2
		noi di as text"# > .................. drop incompatible genotypes"
		rename (a1_1 a2_1) (a1 a2)
		recodegenotype
		rename (a1 a2 recodegenotype_gt) (a1_1 a2_1 gt_1)
		rename (a1_2 a2_2) (a1 a2)
		recodegenotype
		rename (a1 a2 recodegenotype_gt) (a1_2 a2_2 gt_2)
		drop if gt_1 == "K" & gt_2 == "R"
		drop if gt_1 == "K" & gt_2 == "Y"
		drop if gt_1 == "M" & gt_2 == "R"
		drop if gt_1 == "M" & gt_2 == "Y"
		drop if gt_1 == "R" & gt_2 == "K"
		drop if gt_1 == "R" & gt_2 == "M"
		drop if gt_1 == "Y" & gt_2 == "K"
		drop if gt_1 == "Y" & gt_2 == "M"
		noi di as text"# > .................. calculate inverse variance z-score for a1 matched to b1"
		gen w_1   = 1/(se_1^2)
		gen w_2   = 1/(se_2^2)
		gen se_3   = sqrt(1/(w_1 + w_2))
		gen beta_3a = ((beta_1 * w_1 + beta_2 * w_2) / (w_1 + w_2)) 
		gen z_3a    = beta_3a/se_3
		gen p_3a    = 2*(1-normal(abs(-z_3a)))
		noi di as text"# > .................. calculate inverse variance z-score for a1 matched to b2"
		gen beta_2b = -beta_2
		gen beta_3b = ((beta_1 * w_1 + beta_2b * w_2) / (w_1 + w_2)) 
		gen z_3b    = beta_3b/se_3
		gen p_3b    = 2*(1-normal(abs(-z_3b)))
		gen p_3     = p_3a
		noi di as text"# > .................. select strongest p-value for model 3"
		replace p_3  = p_3b if p_3b < p_3a
		qui { // <alternative (skipped) code> calculate from N
			/*
			gen w_1 = sqrt(n_1)
			gen w_2 = sqrt(n_2)
			gen z_2a = z_2
			gen z_3a = (z_1 * w_1 + z_2a * w_2)/sqrt(w_1^2 + w_2^2)
			gen p_3a = (1-normal(abs(-z_3a)))
			gen z_2b = -z_2
			gen z_3b = (z_1 * w_1 + z_2b * w_2)/sqrt(w_1^2 + w_2^2)
			gen p_3b = (1-normal(abs(-z_3b)))
			gen p_3     = p_3a
			replace p_3  = p_3b if p_3b < p_3a		
			*/
			}
		gen beta = .
		gen z    = .
		gen se   = se_3
		gen n    = n_1 + n_2
		foreach i in a b {
			replace beta = beta_3`i' if p_3 == p_3`i'
			replace z    = z_3`i'    if p_3 == p_3`i'
			}
		gen or   = exp(beta)
		gen l95  = or - (1.96 * se)
		gen u95  = or + (1.96 * se)
		gen p = p_3
		tostring chr, replace
		gen model = ""
		replace model = "model-1" if (p_1 <= p_2) & (p_1 <= p_3)
		replace model = "model-2" if (p_2 <= p_1) & (p_2 <= p_3)
		replace model = "model-3" if (p_3 <= p_1) & (p_3 <= p_2)
		save ``function'_project'-`function'-meta.dta, replace
		checkfile, file(`already_created')
		}
	noi di as text"# > .................. "as result"completed"
	noi di as text"#########################################################################"
	noi di as text""
	}
qui { // module 4 - create summaryqc
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 4: create summaryqc"
	noi di as text"#########################################################################"
	local already_created ${init_root}/summaryqc/``function'_project'-model3/``function'_project'-model3-summaryqc.dta
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		}
	else {	
		use ``function'_project'-`function'-meta.dta, clear
		gen a1 = a1_1
		gen a2 = a2_1
		keep  chr bp snp a1 a2 beta z se or l95 u95 n p
		!mkdir ${init_root}/summaryqc/``function'_project'-model3
		save ${init_root}/summaryqc/``function'_project'-model3/``function'_project'-model3-summaryqc.dta, replace
		noi checkfile, file(`already_created')
		}
	noi di as text"# > .................. "as result"completed"
	noi di as text"#########################################################################"
	noi di as text""
	}
qui { // module 5 - select model 3 SNPs and clump
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 4: select model 3 SNPs and clump"
	noi di as text"#########################################################################"
	noi di as text"# > .................. "as result"updating model selection parameters"
	noi di as text"# > . model 3 SNPs P < "as result"`: display %10.4e `p_top''"
	noi di as text"# > .. max_P in models "as result"`: display %10.4e `p_max''"
	clear 
	set obs 1
	gen redo_models = "`redo_models'"
	if redo_models == "yes" {
		files2dta, dir(`c(pwd)')
		gen drop = .
		replace drop = 1 if file == "``function'_project'-`function'.dta"
		keep if drop == 1
		gen script = "erase " + file
		outsheet script using temp.do, non noq replace 
		do temp.do
		erase temp.do
		clear
		}
	local already_created ${init_root}/`function'/``function'_project'/``function'_project'-`function'.dta
	capture confirm file `already_created'
	if !_rc {
		checkfile, file(`already_created')
		}
	else {
		use ${init_root}/`function'/``function'_project'/``function'_project'-`function'-meta.dta, clear
		for var p_1 p_2 p: drop if X > `p_max'
		keep if model == "model-3"
		keep snp chr bp a1_1 a2_1 p
		rename (a1_1 a2_1) (a1 a2)
		checkloc_name
		keep snp p loc
		noi di as text"# > .................. "as result"updating model selection parameters"
		noi di as text"# > .................. "as result"limit to ld-ref SNPs"
		merge 1:1 snp using ${init_root}/bim2dta/``function'_ldref_short'/``function'_ldref_short'_bim.dta
		keep if _m == 3
		renvars, upper
		outsheet using  ``function'_project'-`function'.top, noq replace
		!plink2 --bfile  `ldref' --clump ``function'_project'-`function'.top --clump-allow-overlap --clump-range `glist_hg19' --clump-p1 `p_top' --clump-p2 `p_max' --clump-kb 1000 --clump-r2 0.2 --out ``function'_project'-`function'
		!${tabbed} ``function'_project'-`function'.clumped
		import delim using ${init_root}/`function'/``function'_project'/``function'_project'-`function'.clumped.tabbed,clear varnames(1) case(upper)
		dropmiss, force
		keep CHR SNP BP P
		for var CHR BP P: destring X, replace
		sort SNP
		save ${init_root}/`function'/``function'_project'/``function'_project'-`function'.clumped.dta,replace
		!${tabbed} ${init_root}/`function'/``function'_project'/``function'_project'-`function'.clumped.ranges
		import delim using ${init_root}/`function'/``function'_project'/``function'_project'-`function'.clumped.ranges.tabbed,clear varnames(1) case(upper)
		dropmiss, force
		keep SNP POS KB RANGES
		replace RANGES = subinstr(RANGES , "[", "",.) 
		replace RANGES = subinstr(RANGES , "]", "",.) 
		split RANGES, p(",")
		drop RANGES
		reshape long RANGES, i(SNP POS KB) j(OBS)
		drop OBS
		duplicates drop
		gsort -RANGES
		egen OBS = seq(),by(SNP POS KB)
		replace RANGES = "NONE" if RANGES == "" & OBS == 1
		drop if RANGES == ""
		drop OBS
		sort RANGES
		egen OBS = seq(),by(SNP POS KB)
		reshape wide RANGES , i(SNP POS KB) j(OBS)
		gen RANGES999 = ""
		egen RANGES = concat(RANGES1-RANGES999),p(" ")
		compress
		keep SNP POS KB RANGES
		sort SNP
		merge m:m SNP using ${init_root}/`function'/``function'_project'/``function'_project'-`function'.clumped.dta
		drop _merge
		duplicates drop
		sort SNP
		save ${init_root}/`function'/``function'_project'/``function'_project'-`function'.clumped.ranges.dta,replace
		use ${init_root}/`function'/``function'_project'/``function'_project'-`function'-meta.dta, clear
		keep snp p_1 p_2 p_3 model
		renvars, upper
		compress
		merge m:m SNP using ${init_root}/`function'/``function'_project'/``function'_project'-`function'.clumped.ranges.dta
		ta _merge
		keep if _merge == 3 
		drop _merge
		sort SNP
		tostring CHR, replace
		tostring BP,replace
		gen POSITION = "chr" + CHR + ":" + BP
		drop CHR BP
		gen GENES = RANGES
		replace GENES = subinstr(GENES , "[", "",.) 
		replace GENES = subinstr(GENES , "]", "",.) 
		replace GENES = subinstr(GENES , ",", " ",.) 
		rename POS LD_RANGE
		rename (P_1 P_2 P_3) (P_MODEL1 P_MODEL2 P_MODEL3)
		keep  SNP POSITION P_MODEL1 P_MODEL2 P_MODEL3 LD_RANGE GENES 
		order SNP POSITION P_MODEL1 P_MODEL2 P_MODEL3 LD_RANGE GENES 
		lab var SNP "MARKER"
		lab var POSITION "MARKER LOCATION"
		lab var LD_RANGE "ASSOCIATED REGION (CLUMP)"
		lab var GENES "GENES IN REGION"
		lab var P_MODEL1 "P (``function'_summaryqc1')"
		lab var P_MODEL2 "P (``function'_summaryqc2')"
		lab var P_MODEL3 "P (combined)"
		sort P_MODEL3
		save ${init_root}/`function'/``function'_project'/``function'_project'-`function'.dta, replace
		}
	use ${init_root}/`function'/``function'_project'/``function'_project'-`function'.dta, clear
	count
	local `function'_regions `r(N)'
	noi di as text"# > .. regions flagged "as result"``function'_regions'"
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module x - create scripts for graphmiami
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 7: create supporting scripts"
	noi di as text"#########################################################################"
	noi di as text"# > .................. "as result"create graphmiami script"
	use ${init_root}/`function'/``function'_project'/``function'_project'-`function'.dta, clear
	split LD_RANGE,p("chr"":""..")
	for var LD_RANGE3 LD_RANGE4: destring X, replace
	gen u = .
	gen d = .
	replace u = round(LD_RANGE3 - 1000000,500000)
	replace d = round(LD_RANGE4 + 1000000,500000)
	replace u = u - 500000 if GENES == "NONE"
	replace d = d + 500000 if GENES == "NONE"
	for var u d: tostring X, replace
	gen script = ""
	replace script = "graphmiami, summaryqc1(`summaryqc1') summaryqc2(`summaryqc2') chr(" + LD_RANGE2 + ") from(" + u + ") to(" + d + ") generef(\${graphgene_ref}) ldref(\${graphlocus_ld}) png(yes)"
	noi di as text"# > .......... created "as result"``function'_project'-`function'-graphmiami.do"
	outsheet script using ``function'_project'-`function'-graphmiami.do, non noq replace
	noi di as text"# > .................. "as result"create locuszoom script"
	use ${init_root}/`function'/``function'_project'/``function'_project'-`function'.dta, clear
	split LD_RANGE,p("chr"":""..")
	for var LD_RANGE3 LD_RANGE4: destring X, replace
	gen u = .
	gen d = .
	replace u = round(LD_RANGE3 - 1000000,500000)
	replace d = round(LD_RANGE4 + 1000000,500000)
	replace u = u - 500000 if GENES == "NONE"
	replace d = d + 500000 if GENES == "NONE"
	for var u d: tostring X, replace
	gen script = ""
	replace script = "graphlocuszoom, index(" + SNP + ") generef(\${graphgene_ref}) ldref(\${graphlocus_ld}) recombref(\${graphlocus_recomb}) summaryqc(${init_root}/summaryqc/``function'_project'-model3/``function'_project'-model3) png(yes) range(" + u + "-" + d + ")"
	outsheet script using ``function'_project'-`function'-graphlocuszoom.do, non noq replace
	noi di as text"# > .......... created "as result"``function'_project'-`function'-graphlocuszoom.do"
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
	file write myfile `"# > ...... p-combined `: display %10.4e `p_top''"' _n
	file write myfile `"# > max P in GWAS 1/2 `: display %10.4e `p_max''"' _n
	file write myfile `"# > . regions flagged ``function'_regions'"' _n
	file write myfile `"# > ...... saved data ``function'_project'-`function'.dta"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > overview of `function'"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > 1. merge 2x summaryqc files"' _n
	file write myfile `"# > 2. calculate same and opposite allele meta analysis of 2 datasets "' _n
	file write myfile `"# > 3. create summaryqc file of p-max "' _n
	file write myfile `"# > 4. define models and select on thresholds"' _n
	file write myfile `"# > 5. clump markers using reference genotypes"' _n
	file write myfile `"# > 6. create supporting scripts for downstream plots"' _n
	file write myfile `"#########################################################################"' _n
	file close myfile	
	}
qui { // module z - clean folder
	files2dta, dir(`c(pwd)')
	gen keep = .
	replace keep = 1 if file == "``function'_project'-`function'.log"
	replace keep = 1 if file == "``function'_project'-`function'.dta"
	replace keep = 1 if file == "``function'_project'-`function'-meta.dta"
	replace keep = 1 if file == "``function'_project'-`function'-graphmiami.do"
	replace keep = 1 if file == "``function'_project'-`function'-graphlocuszoom.do"
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
