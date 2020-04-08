*! 1.0.2 Richard Anney 18dec2018
* - add manhattan plot from original (pre-cleaned dataset)
* - add ld_ref to highlight reference is for determining LD in the dataset
* - add routine to perform "minus" regions of high LD


program  summaryqc2score
local    function summaryqc2score
local    `function'_version v2
syntax , summaryqc(string asis) genotypeqc(string asis) project(string asis) ld_ref(string asis) [min(real 3) png(string asis)]

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
	bim_path_to_short, bim(`summaryqc')
	local `function'_summaryqc_short ${bim_short}
	cd ${init_root}	
	!mkdir `function'
	cd `function'
	noi di as text"# > ..... project name "as result"`project'-genotypes-``function'_summaryqc_short'-gwas"
	!mkdir `project'-genotypes
	cd     `project'-genotypes
	!mkdir ``function'_summaryqc_short'-gwas
	cd     ``function'_summaryqc_short'-gwas	
	qui { // define genotypeqc datasets to process
		clear
		set obs 1
		gen genotypeqc = "`genotypeqc'"
		split genotypeqc,p(","" ")
		drop genotypeqc
		sxpose, clear
		count
		local `function'_bim_count `r(N)'
		gen count = _n
		tostring count, replace
		egen obs = seq(),by(_var1)
		sort _var1 obs
		gen script = ""
		replace script = "global `function'_bim" + count + " " + _var1
		outsheet script using `function'script.do, non noq replace
		do `function'script.do	
		forvalues i = 1(1)``function'_bim_count' {
			bim_path_to_short, bim(${`function'_bim`i'})
			local `function'_bim`i'_short ${bim_short}
			}
		erase `function'script.do
		}
	ref_path_to_short, ref(`ld_ref')
	local `function'_ld_ref_short ${ref_short}	
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 2 - define files
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 2: define files"
	noi di as text"#########################################################################"
	noi di as text"# > ......... locating summaryqc file"
	noi checkfile, file(`summaryqc'-summaryqc.dta)
	noi di as text"# > ......... locating PLINK binary files"
	forvalues i = 1(1)``function'_bim_count' {
			foreach file in bim bed fam {
				noi checkfile, file(${`function'_bim`i'}.`file')
				}
			}
	noi di as text"# > ......... locating PLINK files"
	noi checkfile, file(${plink})
	noi checkfile, file(${plink2})
	checktabbed
	noi di as text"# > ......... locating ld reference files"
	foreach file in bim bed fam {
		noi checkfile, file(`ld_ref'.`file')
		}
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 3 - identify commonSNPs list across reference / summaryqc and genotypeqc
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3: make commonSNP list"
	noi di as text"#########################################################################"
	qui { // PLINK binaries
		noi di as text"# > ......... defining overlapping SNPS file from PLINK binaries"
		forvalues i = 1(1)``function'_bim_count' {
			bim2count, bim(${`function'_bim`i'})
			noi di as text"# > ......... binaries "as result"``function'_bim`i'_short'"
			noi di as text"# > ..... SNPs in file "as result "${bim2count_snp}"
			noi di as text"# > .... indiv in file "as result "${bim2count_ind}"
			local `function'_bim`i'N ${bim2count_snp}
			}
		local already_created `project'-genotypes-``function'_summaryqc_short'-gwas-snps.dta
		capture confirm file `already_created'
		if !_rc {
			noi checkfile, file(`already_created')
			use `already_created', clear 
			count 
			local `function'_merge0N `r(N)'
			}
		else {
			noi di as text"# > ......... locating bim.dta files"
			forvalues i = 1(1)``function'_bim_count' {
				local already_created ${init_root}/bim2dta/``function'_bim`i'_short'/``function'_bim`i'_short'_bim.dta
				capture confirm file `already_created'
				if !_rc {
					noi checkfile, file(`already_created')
					}
				else {	
					bim2dta, bim(${`function'_bim`i'})
					noi checkfile, file(`already_created')
					}
				}
			clear
			set obs ``function'_bim_count'
			gen bim_dta = ""
			forvalues i = 1(1)``function'_bim_count' {
				replace bim_dta = "${init_root}/bim2dta/``function'_bim`i'_short'/``function'_bim`i'_short'_bim.dta" in `i'
				}
			gen n = _n
			expand 6
			egen x = seq(),by(bim)
			sort n x
			gen script = ""
			replace script = "use snp gt using " + bim_dta + ", clear" if n == 1 & x == 1
			replace script = `"drop if gt == "W""'  if n == 1 & x == 2
			replace script = `"drop if gt == "S""'  if n == 1 & x == 3
			replace script = `"drop if gt == "ID""' if n == 1 & x == 4
			replace script = "merge 1:1 snp using " + bim_dta  if n != 1 & x == 2
			replace script = "keep if _m == 3" if n != 1 & x == 3
			replace script = "keep snp"        if n != 1 & x == 4
			replace script = "save `project'-genotypes-``function'_summaryqc_short'-gwas-snps.dta,replace" if n == ``function'_bim_count' & x == 5
			drop if script == ""
			outsheet script using `function'script.do, non noq replace
			noi di as text"# > ......... creating overlapping SNPS file from genotypes"
			do `function'script.do
			erase `function'script.do
			noi checkfile, file(`project'-genotypes-``function'_summaryqc_short'-gwas-snps.dta)
			count 
			local `function'_merge0N `r(N)'
			}
			noi di as text"# > . SNPs after merge "as result "``function'_merge0N'"
		}
	qui { // summaryqc
	noi di as text"#########################################################################"
		noi di as text"# > ......... defining overlapping SNPS file from summaryqc"
		noi di as text"# > ........ summaryqc "as result"``function'_summaryqc_short'"
		use `summaryqc'-summaryqc.dta, clear
		count 
		local `function'_sumN `r(N)'
		noi di as text"# > ..... SNPs in file "as result "``function'_sumN'"
		local already_created `project'-genotypes-``function'_summaryqc_short'-gwas-summaryqc.dta
		capture confirm file `already_created'
		if !_rc {
			noi checkfile, file(`already_created')
			use `already_created', clear 
			count 
			local `function'_merge1N `r(N)'
			}
		else {	
			merge 1:1 snp using `project'-genotypes-``function'_summaryqc_short'-gwas-snps.dta
			keep if _m == 3
			count 
			local `function'_merge1N `r(N)'
			drop _m
			save `project'-genotypes-``function'_summaryqc_short'-gwas-summaryqc.dta, replace
			noi checkfile, file(`already_created')
			}
		noi di as text"# > . SNPs after merge "as result "``function'_merge1N'"
		}
	qui { // ld-reference
	noi di as text"#########################################################################"
		noi di as text"# > ......... defining overlapping SNP list in LD reference files"
		bim2count, bim(`ld_ref')
		noi di as text"# > ......... binaries "as result"``function'_ld_ref_short'"
		noi di as text"# > ..... SNPs in file "as result "${bim2count_snp}"
		noi di as text"# > .... indiv in file "as result "${bim2count_ind}"
		local `function'_ld_refN ${bim2count_snp}	
		local already_created `project'-genotypes-``function'_summaryqc_short'-gwas-all.dta
		capture confirm file `already_created'
		if !_rc {	
			noi checkfile, file(`already_created')
			use `already_created', clear 
			count 
			local `function'_merge2N `r(N)'
			}
		else {
			local already_created ${init_root}/bim2dta/``function'_ld_ref_short'/``function'_ld_ref_short'_bim.dta
			capture confirm file `already_created'
			if !_rc {
				noi checkfile, file(`already_created')
				}
			else {	
				bim2dta, bim(`ld_ref')
				noi checkfile, file(`already_created')
				}
			use snp using `already_created', clear
			merge 1:1 snp using `project'-genotypes-``function'_summaryqc_short'-gwas-summaryqc.dta
			keep if _m == 3
			count 
			local `function'_merge2N `r(N)'
			drop _m
			save `project'-genotypes-``function'_summaryqc_short'-gwas-all.dta, replace
			noi checkfile, file(`project'-genotypes-``function'_summaryqc_short'-gwas-all.dta)
			}
		noi di as text"# > . SNPs after merge "as result "``function'_merge2N'"
		}
	noi di as text"#########################################################################"
	noi di as text""
	}
qui { // module 4 - map to common strand
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 4: parse and map to common strand"
	noi di as text"#########################################################################"
	noi di as text"# > .................. load common SNP list"
	local already_created `project'-genotypes-``function'_summaryqc_short'-gwas-strand.dta
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		}
	else {
		use snp a1 a2 using `project'-genotypes-``function'_summaryqc_short'-gwas-all.dta, clear
		noi di as text"# > ......... define genotype code"
		recodegenotype
		rename recodegenotype_gt strand_gt 
		keep snp strand_gt
		save `project'-genotypes-``function'_summaryqc_short'-gwas-strand.dta, replace
		noi checkfile, file(`already_created')
		}
	noi di as text"# > .................. define extract list"
	local already_created `project'-genotypes-``function'_summaryqc_short'-gwas.extract
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		}
	else {
		use snp using `project'-genotypes-``function'_summaryqc_short'-gwas-strand.dta, clear
		outsheet snp using `project'-genotypes-``function'_summaryqc_short'-gwas.extract, non noq replace
		noi checkfile, file(`already_created')
		}
	noi di as text"# > .................. define flip list"
	forvalues i = 1(1)``function'_bim_count' {
		local already_created ``function'_bim`i'_short'-genotypes-``function'_summaryqc_short'-gwas.flip
		capture confirm file `already_created'
		if !_rc {
			noi checkfile, file(`already_created')
			}
		else {
			use snp gt using  ${init_root}/bim2dta/``function'_bim`i'_short'/``function'_bim`i'_short'_bim.dta, clear
			merge 1:1 snp using `project'-genotypes-``function'_summaryqc_short'-gwas-strand.dta
			keep if _m == 3
			outsheet snp if gt != strand_gt using ``function'_bim`i'_short'-genotypes-``function'_summaryqc_short'-gwas.flip, non noq replace
			noi checkfile, file(`already_created')
			}
		}
	noi di as text"# > .................. make binaries limited to common snps"
	forvalues i = 1(1)``function'_bim_count' {
		local already_created ``function'_bim`i'_short'-genotypes-``function'_summaryqc_short'-gwas.bed
		capture confirm file `already_created'
		if !_rc {
			noi checkfile, file(`already_created')
			}
		else {
			create_temp_dir
			!$plink 	--bfile 	  ${`function'_bim`i'} ///
					--extract   ${init_root}/summaryqc2score/`project'-genotypes/``function'_summaryqc_short'-gwas/`project'-genotypes-``function'_summaryqc_short'-gwas.extract ///
					--make-bed  ///
					--out       ``function'_bim`i'_short'
			!$plink 	--bfile     ``function'_bim`i'_short' ///
					--flip      ${init_root}/summaryqc2score/`project'-genotypes/``function'_summaryqc_short'-gwas/``function'_bim`i'_short'-genotypes-``function'_summaryqc_short'-gwas.flip ///
					--make-bed  ///
					--out       ${init_root}/summaryqc2score/`project'-genotypes/``function'_summaryqc_short'-gwas/``function'_bim`i'_short'-genotypes-``function'_summaryqc_short'-gwas
			delete_temp_dir
			noi checkfile, file(`already_created')
			}
		}
	local already_created ``function'_ld_ref_short'-genotypes-``function'_summaryqc_short'-gwas.bed
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		}
	else {
		create_temp_dir
		!$plink 	--bfile 	`ld_ref' ///
				--extract ${init_root}/summaryqc2score/`project'-genotypes/``function'_summaryqc_short'-gwas/`project'-genotypes-``function'_summaryqc_short'-gwas.extract ///
				--make-bed ///
				--out ${init_root}/summaryqc2score/`project'-genotypes/``function'_summaryqc_short'-gwas/``function'_ld_ref_short'-genotypes-``function'_summaryqc_short'-gwas
		delete_temp_dir
		noi checkfile, file(`already_created')
		}
	noi di as text"#########################################################################"
	noi di as text""
	}
qui { // module 5 - define score files (all SNPs)
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 5: define score files"
	noi di as text"#########################################################################"
	qui { // define weights
		noi di as text"# > .................. define weights"
		noi checkfile, file(`project'-genotypes-``function'_summaryqc_short'-gwas-all.dta)
		local already_created `project'-genotypes-``function'_summaryqc_short'-gwas-weights.dta
		capture confirm file `already_created'
		if !_rc {
			noi checkfile, file(`already_created')
			use `already_created', clear
			}
		else {
			noi di as text`"# > .................. define "risk" allele"'
			use `project'-genotypes-``function'_summaryqc_short'-gwas-all.dta, clear		
			gen flip = .
			replace flip = 1 if or < 1
			gen risk_or  = or
			gen risk     = a1
			gen alt      = a2
			replace risk_or    = 1/or   if flip == 1
			replace risk       = a2     if flip == 1
			replace alt        = a1     if flip == 1
			gen weight = ln(risk_or)
			keep snp risk alt weight p 
			renvars, upper
			save `project'-genotypes-``function'_summaryqc_short'-gwas-weights.dta, replace
			capture confirm file `already_created'
			}
		}
	qui { // define pval threshold
		noi di as text"# > .................. define pval thresholds"
		sum P
		local `function'_minP `: display %10.4e r(min)'
		noi di as text"# > ...... min pval is "as result "``function'_minP'"
		if `r(min)' < 1E-8 {
			local `function'_threshold 5E-1 1E-1 5E-2 1E-2 1E-3 1E-4 1E-5 1E-6 1E-7 1E-8
			}
		else if `r(min)' < 1E-7 {
			local `function'_threshold 5E-1 1E-1 5E-2 1E-2 1E-3 1E-4 1E-5 1E-6 1E-7 
			}
		else if `r(min)' < 1E-6 {
			local `function'_threshold 5E-1 1E-1 5E-2 1E-2 1E-3 1E-4 1E-5 1E-6 
			}
		else if `r(min)' < 1E-5 {
			local `function'_threshold 5E-1 1E-1 5E-2 1E-2 1E-3 1E-4 1E-5 
			}
		else if `r(min)' < 1E-4 {
			local `function'_threshold 5E-1 1E-1 5E-2 1E-2 1E-3 1E-4 
			}
		else if `r(min)' < 1E-3 {
			local `function'_threshold 5E-1 1E-1 5E-2 1E-2 1E-3 
			}
		else if `r(min)' < 1E-2 {
			local `function'_threshold 5E-1 1E-1 5E-2 1E-2 
			}
		else if `r(min)' < 1E-1 {
			local `function'_threshold 5E-1 1E-1 
			}
		noi di as text"# > thresholds include "as result "``function'_threshold'"	
		noi di as text"#########################################################################"
		}
	qui { // define clumps
		noi di as text"# > .................. define ld independent snps (clump)"
		noi checkfile, file(``function'_ld_ref_short'-genotypes-``function'_summaryqc_short'-gwas.bed)
		local `function'_prune   --clump-p1 1 --clump-p2 1 --clump-r2 0.2 --clump-kb 1000
		foreach i in ``function'_threshold' {
			local already_created ``function'_ld_ref_short'-genotypes-``function'_summaryqc_short'-gwas-`i'.clumped
			capture confirm file `already_created'
			if !_rc {
				noi di as text"# > .................. clumps defined"
				}
			else {	
				outsheet SNP P if P < `i' using ``function'_ld_ref_short'-genotypes-``function'_summaryqc_short'-gwas-`i'.input-clump, noq replace
				noi di as text"# > .. clumping at p < "as result "`i'"
				!${plink}	--bfile ``function'_ld_ref_short'-genotypes-``function'_summaryqc_short'-gwas 	///
						--clump ``function'_ld_ref_short'-genotypes-``function'_summaryqc_short'-gwas-`i'.input-clump 			///
						``function'_prune' 												///
						--out ``function'_ld_ref_short'-genotypes-``function'_summaryqc_short'-gwas-`i'
				}
			}
		noi di as text"#########################################################################"
		}
	qui { // define score q-score
		noi di as text"# > .................. define score, q-score-file, q-score-ranges files"
		foreach i in ``function'_threshold' {
			local already_created ``function'_ld_ref_short'-genotypes-``function'_summaryqc_short'-gwas-`i'.q-score-range
			capture confirm file `already_created'
			if !_rc {
				import delim using ``function'_ld_ref_short'-genotypes-``function'_summaryqc_short'-gwas-`i'.score, clear varnames(nonames)
				count
				noi di as text"# >  thresholds = `i' "as result"`r(N)' "as text"LD independent SNPs in model"
				if `r(N)' > `min' & `i' == 1E-8 {
					local `function'_threshold 5E-1 1E-1 5E-2 1E-2 1E-3 1E-4 1E-5 1E-6 1E-7 1E-8
					}
				else if `r(N)' > `min' & `i' == 1E-7 {
					local `function'_threshold 5E-1 1E-1 5E-2 1E-2 1E-3 1E-4 1E-5 1E-6 1E-7 
					}
				else if `r(N)' > `min' & `i' == 1E-6 {
					local `function'_threshold 5E-1 1E-1 5E-2 1E-2 1E-3 1E-4 1E-5 1E-6 
					}				
				else if `r(N)' > `min' & `i' == 1E-5 {
					local `function'_threshold 5E-1 1E-1 5E-2 1E-2 1E-3 1E-4 1E-5 
					}				
				else if `r(N)' > `min' & `i' == 1E-4 {
					local `function'_threshold 5E-1 1E-1 5E-2 1E-2 1E-3 1E-4 
					}				
				else if `r(N)' > `min' & `i' == 1E-3 {
					local `function'_threshold 5E-1 1E-1 5E-2 1E-2 1E-3  
					}			
				else if `r(N)' > `min' & `i' == 1E-2 {
					local `function'_threshold 5E-1 1E-1 5E-2 1E-2 
					}			
				else if `r(N)' > `min' & `i' == 5E-2 {
					local `function'_threshold 5E-1 1E-1 5E-2 
					}
				else if `r(N)' > `min' & `i' == 1E-1 {
					local `function'_threshold 5E-1 1E-1 
					}			
				else if `r(N)' > `min' & `i' == 5E-1 {
					local `function'_threshold 5E-1 
					}			
				}
			else {	
				!${tabbed}         ``function'_ld_ref_short'-genotypes-``function'_summaryqc_short'-gwas-`i'.clumped
				import delim using ``function'_ld_ref_short'-genotypes-``function'_summaryqc_short'-gwas-`i'.clumped.tabbed, varnames(1) clear 
				renvars, upper
				merge 1:1 SNP using `project'-genotypes-``function'_summaryqc_short'-gwas-weights.dta
				keep if _m == 3
				outsheet SNP RISK WEIGHT using ``function'_ld_ref_short'-genotypes-``function'_summaryqc_short'-gwas-`i'.score, non noq replace
				outsheet SNP P           using ``function'_ld_ref_short'-genotypes-``function'_summaryqc_short'-gwas-`i'.q-score-file, non noq replace
				count
				noi di as text"# >  thresholds = `i' "as result"`r(N)' "as text"LD independent SNPs in model"
				if `r(N)' > `min' & `i' == 1E-8 {
					local `function'_threshold 5E-1 1E-1 5E-2 1E-2 1E-3 1E-4 1E-5 1E-6 1E-7 1E-8
					}
				else if `r(N)' > `min' & `i' == 1E-7 {
					local `function'_threshold 5E-1 1E-1 5E-2 1E-2 1E-3 1E-4 1E-5 1E-6 1E-7 
					}
				else if `r(N)' > `min' & `i' == 1E-6 {
					local `function'_threshold 5E-1 1E-1 5E-2 1E-2 1E-3 1E-4 1E-5 1E-6 
					}				
				else if `r(N)' > `min' & `i' == 1E-5 {
					local `function'_threshold 5E-1 1E-1 5E-2 1E-2 1E-3 1E-4 1E-5 
					}				
				else if `r(N)' > `min' & `i' == 1E-4 {
					local `function'_threshold 5E-1 1E-1 5E-2 1E-2 1E-3 1E-4 
					}				
				else if `r(N)' > `min' & `i' == 1E-3 {
					local `function'_threshold 5E-1 1E-1 5E-2 1E-2 1E-3  
					}			
				else if `r(N)' > `min' & `i' == 1E-2 {
					local `function'_threshold 5E-1 1E-1 5E-2 1E-2 
					}			
				else if `r(N)' > `min' & `i' == 5E-2 {
					local `function'_threshold 5E-1 1E-1 5E-2 
					}
				else if `r(N)' > `min' & `i' == 1E-1 {
					local `function'_threshold 5E-1 1E-1 
					}			
				else if `r(N)' > `min' & `i' == 5E-1 {
					local `function'_threshold 5E-1 
					}	
				clear
				set obs 1
				gen a = "P`i'	0	`i'"
				outsheet a using ``function'_ld_ref_short'-genotypes-``function'_summaryqc_short'-gwas-`i'.q-score-range, non noq replace
				}
			}
		noi di as text"# >  min SNPs in model "as result "`min'"
		noi di as text"# > updated thresholds "as result "``function'_threshold'"	
		noi di as text"#########################################################################"
		}
	qui { // calculate profile score
		noi di as text"# > .................. calculate profile scores"
		forvalues i = 1(1)``function'_bim_count' {
			noi di as text"# > ......... binaries "as result"``function'_bim`i'_short'"
			local already_created ``function'_bim`i'_short'-genotypes-``function'_summaryqc_short'-gwas-profiles.dta
			capture confirm file `already_created'
			if !_rc {
				noi checkfile, file(`already_created')
				use `already_created', clear
				}
			else {
				noi di as text"# > .................. calculate profile scores using PLINK across thresholds"
				foreach j in ``function'_threshold' {
					!${plink} --bfile 			``function'_bim`i'_short'-genotypes-``function'_summaryqc_short'-gwas					///
							--score           	``function'_ld_ref_short'-genotypes-``function'_summaryqc_short'-gwas-`j'.score        	///
							--q-score-file    	``function'_ld_ref_short'-genotypes-``function'_summaryqc_short'-gwas-`j'.q-score-file  	///
							--q-score-range   	``function'_ld_ref_short'-genotypes-``function'_summaryqc_short'-gwas-`j'.q-score-range 	///
							--out            	``function'_bim`i'_short'-genotypes-``function'_summaryqc_short'-gwas
					clear 
					set obs 1
					gen p = `j'
					if p == 1E-8 {
						local `function'_tag  p1e_8
						}
					else if p == 1E-7 {
						local `function'_tag  p1e_7
						}
					else if p == 1E-6 {
						local `function'_tag  p1e_6
						}
					else if p == 1E-5 {
						local `function'_tag  p1e_5
						}
					else if p == 1E-4 {
						local `function'_tag  p1e_4
						}
					else if p == 1E-3 {
						local `function'_tag  p1e_3
						}
					else if p == 1E-2 {
						local `function'_tag  p1e_2
						}
					else if p == 5E-2 {
						local `function'_tag  p5e_2
						}
					else if p == 1E-1 {
						local `function'_tag  p1e_1
						}
					else if p == 5E-1 {
						local `function'_tag  p5e_1
						}
					di "``function'_tag'"
					!$tabbed           ``function'_bim`i'_short'-genotypes-``function'_summaryqc_short'-gwas.P`j'.profile
					import delim using ``function'_bim`i'_short'-genotypes-``function'_summaryqc_short'-gwas.P`j'.profile.tabbed, case(lower) clear
					keep fid iid cnt cnt2 score
					for var fid iid: tostring X, replace
					for var cnt cnt2 score: rename X ``function'_tag'_X
					save ``function'_bim`i'_short'-genotypes-``function'_summaryqc_short'-gwas.P`j'.profiles.dta,replace
					}
				noi di as text"# > .................. combine profile scores"
				local already_created ${init_root}/bim2eigenvec/``function'_bim`i'_short'/``function'_bim`i'_short'_eigenvec.dta
				capture confirm file `already_created'
				if !_rc {
					noi checkfile, file(`already_created')
					use `already_created', clear
					}
				else {	
					bim2eigenvec, bim(${`function'_bim`i'})
					noi checkfile, file(`already_created')
					use `already_created', clear
					}
				local already_created ${init_root}/fam2dta/``function'_bim`i'_short'/``function'_bim`i'_short'_fam.dta
				capture confirm file `already_created'
				if !_rc {
					noi checkfile, file(`already_created')
					use `already_created', clear
					}
				else {	
					fam2dta, fam(${`function'_bim`i'})
					noi checkfile, file(`already_created')
					use `already_created', clear
					}	
				use ${init_root}/bim2eigenvec/``function'_bim`i'_short'/``function'_bim`i'_short'_eigenvec.dta, clear
				merge 1:1 fid iid using  ${init_root}/fam2dta/``function'_bim`i'_short'/``function'_bim`i'_short'_fam.dta
				keep if _m == 3
				keep fid iid sex pc1-pc10
				foreach j in ``function'_threshold' {
					merge 1:1 fid iid using 	``function'_bim`i'_short'-genotypes-``function'_summaryqc_short'-gwas.P`j'.profiles.dta
					keep if _m == 3
					drop _m
					}
				aorder
				order fid iid sex pc1-pc10
				save	``function'_bim`i'_short'-genotypes-``function'_summaryqc_short'-gwas-profiles.dta,replace
				outsheet using ``function'_bim`i'_short'-genotypes-``function'_summaryqc_short'-gwas-profiles.csv, noq replace
				noi checkfile, file(``function'_bim`i'_short'-genotypes-``function'_summaryqc_short'-gwas-profiles.dta)
				noi di as text"#########################################################################"
				}
			}
		}
	noi di as text"#########################################################################"
	noi di as text""
	}
qui { // module 6 - define score files (minus regions of high LD (e.g. HLA)
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 6: define score files (minus HLA)"
	noi di as text"#########################################################################"
	qui { // define regions to drop
		local already_created ${init_root}/bim2ldexclude/``function'_ld_ref_short'-genotypes-``function'_summaryqc_short'-gwas/``function'_ld_ref_short'-genotypes-``function'_summaryqc_short'-gwas-bim2ldexclude.exclude
		capture confirm file `already_created'
		if !_rc {
			noi checkfile, file(`already_created')
			}
		else {
			bim2ldexclude, bim(${init_root}/summaryqc2score/`project'-genotypes/``function'_summaryqc_short'-gwas/``function'_ld_ref_short'-genotypes-``function'_summaryqc_short'-gwas)
			}
		local `function'_prune --clump-p1 1 --clump-p2 1 --clump-r2 0.2 --clump-kb 1000 --exclude `already_created'
		}
	qui { // define clumps
		noi di as text"# > .................. define ld independent snps (clump)"
		noi checkfile, file(``function'_ld_ref_short'-genotypes-``function'_summaryqc_short'-gwas.bed)
		foreach i in ``function'_threshold' {
			local already_created ``function'_ld_ref_short'-genotypes-``function'_summaryqc_short'-gwas-`i'-nohla.clumped
			capture confirm file `already_created'
			if !_rc {
				noi di as text"# > .................. clumps defined"
				}
			else {	
				noi di as text"# > .. clumping at p < "as result "`i'"
				!${plink}	--bfile ``function'_ld_ref_short'-genotypes-``function'_summaryqc_short'-gwas 	///
						--clump ``function'_ld_ref_short'-genotypes-``function'_summaryqc_short'-gwas-`i'.input-clump 			///
						``function'_prune' 												///
						--out ``function'_ld_ref_short'-genotypes-``function'_summaryqc_short'-gwas-`i'-nohla
				}
			}
		noi di as text"#########################################################################"
		}
	qui { // define score q-score
		noi di as text"# > .................. define score, q-score-file, q-score-ranges files"
		foreach i in ``function'_threshold' {
			local already_created ``function'_ld_ref_short'-genotypes-``function'_summaryqc_short'-gwas-`i'-nohla.q-score-range
			capture confirm file `already_created'
			if !_rc {
				import delim using ``function'_ld_ref_short'-genotypes-``function'_summaryqc_short'-gwas-`i'-nohla.score, clear varnames(nonames)
				count
				noi di as text"# >  thresholds = `i' "as result"`r(N)' "as text"LD independent SNPs in model"
				if `r(N)' > `min' & `i' == 1E-8 {
					local `function'_threshold2 5E-1 1E-1 5E-2 1E-2 1E-3 1E-4 1E-5 1E-6 1E-7 1E-8
					}
				else if `r(N)' > `min' & `i' == 1E-7 {
					local `function'_threshold2 5E-1 1E-1 5E-2 1E-2 1E-3 1E-4 1E-5 1E-6 1E-7 
					}
				else if `r(N)' > `min' & `i' == 1E-6 {
					local `function'_threshold2 5E-1 1E-1 5E-2 1E-2 1E-3 1E-4 1E-5 1E-6 
					}				
				else if `r(N)' > `min' & `i' == 1E-5 {
					local `function'_threshold2 5E-1 1E-1 5E-2 1E-2 1E-3 1E-4 1E-5 
					}				
				else if `r(N)' > `min' & `i' == 1E-4 {
					local `function'_threshold2 5E-1 1E-1 5E-2 1E-2 1E-3 1E-4 
					}				
				else if `r(N)' > `min' & `i' == 1E-3 {
					local `function'_threshold2 5E-1 1E-1 5E-2 1E-2 1E-3  
					}			
				else if `r(N)' > `min' & `i' == 1E-2 {
					local `function'_threshold2 5E-1 1E-1 5E-2 1E-2 
					}			
				else if `r(N)' > `min' & `i' == 5E-2 {
					local `function'_threshold2 5E-1 1E-1 5E-2 
					}
				else if `r(N)' > `min' & `i' == 1E-1 {
					local `function'_threshold2 5E-1 1E-1 
					}			
				else if `r(N)' > `min' & `i' == 5E-1 {
					local `function'_threshold2 5E-1 
					}			
				}
			else {	
				!${tabbed}         ``function'_ld_ref_short'-genotypes-``function'_summaryqc_short'-gwas-`i'-nohla.clumped
				import delim using ``function'_ld_ref_short'-genotypes-``function'_summaryqc_short'-gwas-`i'-nohla.clumped.tabbed, varnames(1) clear 
				renvars, upper
				merge 1:1 SNP using `project'-genotypes-``function'_summaryqc_short'-gwas-weights.dta
				keep if _m == 3
				outsheet SNP RISK WEIGHT using ``function'_ld_ref_short'-genotypes-``function'_summaryqc_short'-gwas-`i'-nohla.score, non noq replace
				outsheet SNP P           using ``function'_ld_ref_short'-genotypes-``function'_summaryqc_short'-gwas-`i'-nohla.q-score-file, non noq replace
				count
				noi di as text"# >  thresholds = `i' "as result"`r(N)' "as text"LD independent SNPs in model"
				if `r(N)' > `min' & `i' == 1E-8 {
					local `function'_threshold2 5E-1 1E-1 5E-2 1E-2 1E-3 1E-4 1E-5 1E-6 1E-7 1E-8
					}
				else if `r(N)' > `min' & `i' == 1E-7 {
					local `function'_threshold2 5E-1 1E-1 5E-2 1E-2 1E-3 1E-4 1E-5 1E-6 1E-7 
					}
				else if `r(N)' > `min' & `i' == 1E-6 {
					local `function'_threshold2 5E-1 1E-1 5E-2 1E-2 1E-3 1E-4 1E-5 1E-6 
					}				
				else if `r(N)' > `min' & `i' == 1E-5 {
					local `function'_threshold2 5E-1 1E-1 5E-2 1E-2 1E-3 1E-4 1E-5 
					}				
				else if `r(N)' > `min' & `i' == 1E-4 {
					local `function'_threshold2 5E-1 1E-1 5E-2 1E-2 1E-3 1E-4 
					}				
				else if `r(N)' > `min' & `i' == 1E-3 {
					local `function'_threshold2 5E-1 1E-1 5E-2 1E-2 1E-3  
					}			
				else if `r(N)' > `min' & `i' == 1E-2 {
					local `function'_threshold2 5E-1 1E-1 5E-2 1E-2 
					}			
				else if `r(N)' > `min' & `i' == 5E-2 {
					local `function'_threshold2 5E-1 1E-1 5E-2 
					}
				else if `r(N)' > `min' & `i' == 1E-1 {
					local `function'_threshold2 5E-1 1E-1 
					}			
				else if `r(N)' > `min' & `i' == 5E-1 {
					local `function'_threshold2 5E-1 
					}	
				clear
				set obs 1
				gen a = "P`i'	0	`i'"
				outsheet a using ``function'_ld_ref_short'-genotypes-``function'_summaryqc_short'-gwas-`i'-nohla.q-score-range, non noq replace
				}
			}
		noi di as text"# >  min SNPs in model "as result "`min'"
		noi di as text"# > updated thresholds "as result "``function'_threshold'"	
		noi di as text"#########################################################################"
		}
	qui { // calculate profile score
		noi di as text"# > .................. calculate profile scores"
		forvalues i = 1(1)``function'_bim_count' {
			noi di as text"# > ......... binaries "as result"``function'_bim`i'_short'"
			local already_created ``function'_bim`i'_short'-genotypes-``function'_summaryqc_short'-gwas-nohla-profiles.dta
			capture confirm file `already_created'
			if !_rc {
				noi checkfile, file(`already_created')
				use `already_created', clear
				}
			else {
				noi di as text"# > .................. calculate profile scores using PLINK across thresholds"
				foreach j in ``function'_threshold' {
					!${plink} --bfile 			``function'_bim`i'_short'-genotypes-``function'_summaryqc_short'-gwas					///
							--score           	``function'_ld_ref_short'-genotypes-``function'_summaryqc_short'-gwas-`j'-nohla.score        	///
							--q-score-file    	``function'_ld_ref_short'-genotypes-``function'_summaryqc_short'-gwas-`j'-nohla.q-score-file  	///
							--q-score-range   	``function'_ld_ref_short'-genotypes-``function'_summaryqc_short'-gwas-`j'-nohla.q-score-range 	///
							--out            	``function'_bim`i'_short'-genotypes-``function'_summaryqc_short'-gwas-nohla
					clear 
					set obs 1
					gen p = `j'
					if p == 1E-8 {
						local `function'_tag  p1e_8
						}
					else if p == 1E-7 {
						local `function'_tag  p1e_7
						}
					else if p == 1E-6 {
						local `function'_tag  p1e_6
						}
					else if p == 1E-5 {
						local `function'_tag  p1e_5
						}
					else if p == 1E-4 {
						local `function'_tag  p1e_4
						}
					else if p == 1E-3 {
						local `function'_tag  p1e_3
						}
					else if p == 1E-2 {
						local `function'_tag  p1e_2
						}
					else if p == 5E-2 {
						local `function'_tag  p5e_2
						}
					else if p == 1E-1 {
						local `function'_tag  p1e_1
						}
					else if p == 5E-1 {
						local `function'_tag  p5e_1
						}
					di "``function'_tag'"
					!$tabbed           ``function'_bim`i'_short'-genotypes-``function'_summaryqc_short'-gwas-nohla.P`j'.profile
					import delim using ``function'_bim`i'_short'-genotypes-``function'_summaryqc_short'-gwas-nohla.P`j'.profile.tabbed, case(lower) clear
					keep fid iid cnt cnt2 score
					for var fid iid: tostring X, replace
					for var cnt cnt2 score: rename X ``function'_tag'_X
					save ``function'_bim`i'_short'-genotypes-``function'_summaryqc_short'-gwas-nohla.P`j'.profiles.dta,replace
					}
				noi di as text"# > .................. combine profile scores"
				use ${init_root}/bim2eigenvec/``function'_bim`i'_short'/``function'_bim`i'_short'_eigenvec.dta, clear
				merge 1:1 fid iid using  ${init_root}/fam2dta/``function'_bim`i'_short'/``function'_bim`i'_short'_fam.dta
				keep if _m == 3
				keep fid iid sex pc1-pc10
				foreach j in ``function'_threshold' {
					merge 1:1 fid iid using 	``function'_bim`i'_short'-genotypes-``function'_summaryqc_short'-gwas-nohla.P`j'.profiles.dta
					keep if _m == 3
					drop _m
					}
				aorder
				order fid iid sex pc1-pc10
				save	``function'_bim`i'_short'-genotypes-``function'_summaryqc_short'-gwas-nohla-profiles.dta,replace
				outsheet using ``function'_bim`i'_short'-genotypes-``function'_summaryqc_short'-gwas-nohla-profiles.csv, noq replace
				noi checkfile, file(``function'_bim`i'_short'-genotypes-``function'_summaryqc_short'-gwas-nohla-profiles.dta)
				noi di as text"#########################################################################"
				}
			}
		}
	noi di as text"#########################################################################"
	noi di as text""
	}
qui { // module 7 - manhattan of intersect
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 7: plotting manhattan "
	noi di as text"#########################################################################"
	qui { // graph ``function'_summaryqc_short'_`function'-graphmanhattan
		local already_created `project'-genotypes-``function'_summaryqc_short'-gwas-`function'_graphmanhattan.gph
		capture confirm file `already_created'
		if !_rc {
			noi checkfile, file(`already_created')
			}	
		else {	
			local already_created ${init_root}/graphmanhattan/``function'_summaryqc_short'/``function'_summaryqc_short'_graphmanhattan.gph
			capture confirm file `already_created'
			if !_rc {
				noi checkfile, file(`already_created')
				}
			else {
				use `summaryqc'-summaryqc.dta, clear
				graphmanhattan, project(`summaryqc') png(`png')
				noi checkfile, file(`already_created')
				}
			local already_created ${init_root}/graphmanhattan/`project'-genotypes-``function'_summaryqc_short'-gwas/`project'-genotypes-``function'_summaryqc_short'-gwas_graphmanhattan.gph
			capture confirm file `already_created'
			if !_rc {
				noi checkfile, file(`already_created')
				}
			else {
				use `project'-genotypes-``function'_summaryqc_short'-gwas-all.dta, clear
				graphmanhattan, project(`project'-genotypes-``function'_summaryqc_short'-gwas) png(`png')
				noi checkfile, file(`already_created')
				}
			graph combine ///
			${init_root}/graphmanhattan/``function'_summaryqc_short'/``function'_summaryqc_short'_graphmanhattan.gph ///
			${init_root}/graphmanhattan/`project'-genotypes-``function'_summaryqc_short'-gwas/`project'-genotypes-``function'_summaryqc_short'-gwas_graphmanhattan.gph, ///
			col(1) nodraw saving(`project'-genotypes-``function'_summaryqc_short'-gwas-`function'_graphmanhattan.gph, replace)
			}
		}
	qui { // graph ``function'_summaryqc_short'-nohla_`function'-graphmanhattan
		local already_created `project'-genotypes-``function'_summaryqc_short'-gwas-nohla-`function'_graphmanhattan.gph
		capture confirm file `already_created'
		if !_rc {
			noi checkfile, file(`already_created')
			local already_created ${init_root}/bim2ldexclude/``function'_ld_ref_short'-genotypes-``function'_summaryqc_short'-gwas/``function'_ld_ref_short'-genotypes-``function'_summaryqc_short'-gwas-bim2ldexclude.exclude
			import delim using `already_created', clear varnames(nonames)
			rename v1 snp
			merge 1:1 snp using `project'-genotypes-``function'_summaryqc_short'-gwas-all.dta
			drop if _m == 3
			sum p
			local `function'_minP2 `: display %10.4e r(min)'
			}	
		else {
			local already_created ${init_root}/graphmanhattan/`project'-genotypes-``function'_summaryqc_short'-gwas-nohla/`project'-genotypes-``function'_summaryqc_short'-gwas-nohla_graphmanhattan.gph
			capture confirm file `already_created'
			if !_rc {
				noi checkfile, file(`already_created')
				local already_created ${init_root}/bim2ldexclude/``function'_ld_ref_short'-genotypes-``function'_summaryqc_short'-gwas/``function'_ld_ref_short'-genotypes-``function'_summaryqc_short'-gwas-bim2ldexclude.exclude
				import delim using `already_created', clear varnames(nonames)
				rename v1 snp
				merge 1:1 snp using `project'-genotypes-``function'_summaryqc_short'-gwas-all.dta
				drop if _m == 3
				sum p
				local `function'_minP2 `: display %10.4e r(min)'
				}
			else {
				local already_created ${init_root}/bim2ldexclude/``function'_ld_ref_short'-genotypes-``function'_summaryqc_short'-gwas/``function'_ld_ref_short'-genotypes-``function'_summaryqc_short'-gwas-bim2ldexclude.exclude
				capture confirm file `already_created'
				if !_rc {
					noi checkfile, file(`already_created')
					}
				else {
					bim2ldexclude, bim(${init_root}/summaryqc2score/`project'-genotypes/``function'_summaryqc_short'-gwas/${init_root}/bim2ldexclude/``function'_ld_ref_short'-genotypes-``function'_summaryqc_short'-gwas/``function'_ld_ref_short'-genotypes-``function'_summaryqc_short'-gwas)
					}
				import delim using `already_created', clear varnames(nonames)
				rename v1 snp
				merge 1:1 snp using `project'-genotypes-``function'_summaryqc_short'-gwas-all.dta
				drop if _m == 3
				sum p
				local `function'_minP2 `: display %10.4e r(min)'
				graphmanhattan, project(`project'-genotypes-``function'_summaryqc_short'-gwas-nohla) png(`png')
				noi checkfile, file(`already_created')
				}		
			graph combine ///
			${init_root}/graphmanhattan/``function'_summaryqc_short'/``function'_summaryqc_short'_graphmanhattan.gph ///
			${init_root}/graphmanhattan/`project'-genotypes-``function'_summaryqc_short'-gwas-nohla/`project'-genotypes-``function'_summaryqc_short'-gwas-nohla_graphmanhattan.gph, ///
			col(1) nodraw saving(`project'-genotypes-``function'_summaryqc_short'-gwas-nohla-`function'_graphmanhattan.gph, replace)
			}
		}
	clear 
	set obs 1
	gen png = "`png'"
	if png == "yes" {
		local already_created `project'-genotypes-``function'_summaryqc_short'-gwas-`function'_graphmanhattan.png
		capture confirm file `already_created'
		if !_rc {
			noi checkfile, file(`already_created')
			}	
		else {
			graph use    `project'-genotypes-``function'_summaryqc_short'-gwas-`function'_graphmanhattan.gph
			graph export `project'-genotypes-``function'_summaryqc_short'-gwas-`function'_graphmanhattan.png, as(png) height(2500) width(4000) replace
			window manage close graph
			}
		local already_created `project'-genotypes-``function'_summaryqc_short'-gwas-nohla-`function'_graphmanhattan.png
		capture confirm file `already_created'
		if !_rc {
			noi checkfile, file(`already_created')
			}
		else {
			graph use    `project'-genotypes-``function'_summaryqc_short'-gwas-nohla-`function'_graphmanhattan.gph
			graph export `project'-genotypes-``function'_summaryqc_short'-gwas-nohla-`function'_graphmanhattan.png, as(png) height(2500) width(4000) replace
			window manage close graph
			}
		}
	else {
		}
	noi di as text"#########################################################################"
	noi di as text""
	}
qui { // module 8 - write log file
	file open myfile using "`project'-genotypes-``function'_summaryqc_short'-gwas-`function'.log", write replace
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
	file write myfile `"# > ..... project name `project'-genotypes-``function'_summaryqc_short'-gwas"' _n
	file write myfile `"# > .... plink version ${plink}"' _n
	file write myfile `"# > .... plink version ${plink2}"' _n
	file write myfile `"#########################################################################"' _n
	forvalues i = 1(1)``function'_bim_count' {
	file write myfile `"# > ............. bim`i' ``function'_bim`i'_short'"' _n
	}
	file write myfile `"# > ..... ld reference ``function'_ld_ref_short'"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > models"' _n
	file write myfile `"#########################################################################"' _n
	forvalues i = 1(1)``function'_bim_count' {
	file write myfile `"# > ..... SNPs in bim`i' ``function'_bim`i'N'"' _n
	}
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > ........ summaryqc ``function'_summaryqc_short'"' _n
	file write myfile `"# > .. SNPs in dataset ``function'_sumN'"' _n
	file write myfile `"# > .. SNPs post merge ``function'_merge2N'"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > .................. `project'-genotypes-``function'_summaryqc_short'-gwas"' _n
	file write myfile `"# > ............. minP ``function'_minP'"' _n
	file write myfile `"# > P thresholds model ``function'_threshold'"' _n
	foreach i in ``function'_threshold' {
	import delim using ``function'_ld_ref_short'-genotypes-``function'_summaryqc_short'-gwas-`i'.score, clear varnames(nonames)
	count
	file write myfile `"# > .... SNPs in model `r(N)' for P < `i'"' _n
	}		
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > .................. `project'-genotypes-``function'_summaryqc_short'-gwas-nohla"' _n
	file write myfile `"# > ............. minP ``function'_minP2'"' _n
	file write myfile `"# > P thresholds model ``function'_threshold2'"' _n
	foreach i in ``function'_threshold' {
	import delim using ``function'_ld_ref_short'-genotypes-``function'_summaryqc_short'-gwas-`i'-nohla.score, clear varnames(nonames)
	count
	file write myfile `"# > .... SNPs in model `r(N)' for P < `i'"' _n
	}		
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > datasets created (.dta / .csv)"' _n
	file write myfile `"#########################################################################"' _n	
	forvalues i = 1(1)``function'_bim_count' {
	file write myfile `"# > ....... saved data ``function'_bim`i'_short'-genotypes-``function'_summaryqc_short'-gwas-profiles.dta"' _n
	}
	file write myfile `"# > ..... saved graphs `project'-genotypes-``function'_summaryqc_short'-gwas-`function'_graphmanhattan.gph"' _n
	forvalues i = 1(1)``function'_bim_count' {
	file write myfile `"# > ....... saved data ``function'_bim`i'_short'-genotypes-``function'_summaryqc_short'-gwas-nohla-profiles.dta"' _n
	}
	file write myfile `"# > ..... saved graphs `project'-genotypes-``function'_summaryqc_short'-gwas-nohla-`function'_graphmanhattan.gph"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > overview of `function'"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > 1. define bim files"' _n
	file write myfile `"# > 2. identify common SNP list across summaryqc, genotypeqc and ld-reference files"' _n
	file write myfile `"# > 3. parse datasets to common SNPs"' _n
	file write myfile `"# > 4. map genotypeqc to common strand"' _n
	file write myfile `"# > 5. define pval weights (ln(OR))"' _n
	file write myfile `"# > 6. define pval threshold in common summaryqc"' _n
	file write myfile `"# > 7. define ld independent SNPsc"' _n
	file write myfile `"# > 8. define score, q-score-file, q-score-ranges files"' _n
	file write myfile `"# > 9. calculate profiles"' _n
	file write myfile `"# > 10. combine profiles (merge with eigenvec)"' _n
	file write myfile `"# > 11. repeat 6-10 with high LD regions removed"' _n
	file write myfile `"# > 12. plot manhattan for parsed dataset"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile "# *.profiles.* file information"                                         _n
	file write myfile "#########################################################################" _n
	file write myfile "# fid .......... (family identifier) ............................ string"      _n
	file write myfile "# iid .......... (individual identifier) ........................ string"      _n
	file write myfile "# sex .......... (sex) .......................................... 1 = male; 2= female" _n
	file write myfile "# *threshold*_cnt .... (number of alleles present in the model) . numeric"     _n
	file write myfile "# *threshold*_cnt2 ... (total number of named alleles observed) . numeric"     _n
	file write myfile "# *threshold*_score .. (weighted score) ......................... numeric"     _n
	file write myfile "#########################################################################" _n
	file write myfile "# Scores were calculated using PLINK. " _n
	file write myfile "# Scores were calculated using mean imputation for missing SNPs " _n
	file write myfile "# A minimum of `min' LD-indep. SNPs are required per threshold to generate a " _n
	file write myfile "# score.                             " _n
	file write myfile "#########################################################################" _n		
	file close myfile	
	}
qui { // module 9 - clean folder
	cd ${init_root}/`function'/`project'-genotypes/``function'_summaryqc_short'-gwas/
	files2dta, dir(${init_root}/`function'/`project'-genotypes/``function'_summaryqc_short'-gwas/)
	!tar -czvf `project'-genotypes-``function'_summaryqc_short'-gwas-raw.tar.gz *.profile *.score *.q-score-range *.q-score-file *.bed *.bim *.fam
	split files, p("-`function'")
	keep if files2 == ""
	drop files1 files2
	gen keep = .
	forvalues i = 1(1)``function'_bim_count' {
		replace keep = 1 if files == "``function'_bim`i'_short'-genotypes-``function'_summaryqc_short'-gwas-profiles.dta"
		replace keep = 1 if files == "``function'_bim`i'_short'-genotypes-``function'_summaryqc_short'-gwas-profiles.csv"	
		replace keep = 1 if files == "``function'_bim`i'_short'-genotypes-``function'_summaryqc_short'-gwas-nohla-profiles.dta"
		replace keep = 1 if files == "``function'_bim`i'_short'-genotypes-``function'_summaryqc_short'-gwas-nohla-profiles.csv"
		}
	split files,p(".")
	replace keep = 1 if files == "`project'-genotypes-``function'_summaryqc_short'-gwas-`function'.log"
	drop if keep == 1
	gen script = "erase " + files
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
	
