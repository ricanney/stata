*! 1.0.2 Richard Anney 17apr2019
* clean up-front functions
*! 1.0.1 Richard Anney 21dec2018
* updated to deal with non ACGT variants
* update to expect a1 a2

program  recodegenotype
local    function recodegenotype
local    `function'_version v1
syntax 

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
qui { // module 1 - check allele variables
	noi di as text""	
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 1: check allele variables"
	noi di as text"#########################################################################"
	noi di as text"# > .......... checking that alleles are defined correctly "
	foreach i in a1 a2 { 
		capture confirm var `i'
		if _rc==0 {
			noi di as text"# > ....... `i' "as result"present"
			}
		else {
			noi di as text"# > ....... `i' "as error"absent"
			exit
			}
		}
	noi di as text"#########################################################################"
	}
qui { // module 2 - converting to genotype codes
	noi di as text""	
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 2: converting to genotype codes"
	noi di as text"#########################################################################"
	count
	noi di as text"# > .......... SNPs in dataset (N) "as result `r(N)'
	qui { // convert to indel
		noi di as text"# > .......... "as result"convert alleles to indel (I/D)"
		gen recodegenotype_counta1 = length(a1)
		gen recodegenotype_counta2 = length(a2)
		replace a1 = "I" if recodegenotype_counta1 > recodegenotype_counta2
		replace a2 = "I" if recodegenotype_counta2 > recodegenotype_counta1
		replace a1 = "D" if a2 == "I"
		replace a2 = "D" if a1 == "I"
		compress
		}
	qui { // generate recodegenotype_gt
		gen recodegenotype_gt = a1 + "|" + a2
		replace recodegenotype_gt = "M"  if recodegenotype_gt == "A|C"
		replace recodegenotype_gt = "R"  if recodegenotype_gt == "A|G"
		replace recodegenotype_gt = "W"  if recodegenotype_gt == "A|T"
		replace recodegenotype_gt = "M"  if recodegenotype_gt == "C|A"
		replace recodegenotype_gt = "S"  if recodegenotype_gt == "C|G"
		replace recodegenotype_gt = "Y"  if recodegenotype_gt == "C|T"
		replace recodegenotype_gt = "R"  if recodegenotype_gt == "G|A"
		replace recodegenotype_gt = "S"  if recodegenotype_gt == "G|C"
		replace recodegenotype_gt = "K"  if recodegenotype_gt == "G|T"
		replace recodegenotype_gt = "W"  if recodegenotype_gt == "T|A"
		replace recodegenotype_gt = "Y"  if recodegenotype_gt == "T|C"
		replace recodegenotype_gt = "K"  if recodegenotype_gt == "T|G"
		replace recodegenotype_gt = "ID" if recodegenotype_gt == "I|D"
		replace recodegenotype_gt = "ID" if recodegenotype_gt == "D|I"
		lab var recodegenotype_gt "Genotype (RYMKSW(ID))"
		}
	qui { // identify non-standard genotypes
		gen recodegenotype_error = 1
		replace recodegenotype_error = . if recodegenotype_gt == "R"
		replace recodegenotype_error = . if recodegenotype_gt == "Y"
		replace recodegenotype_error = . if recodegenotype_gt == "M"
		replace recodegenotype_error = . if recodegenotype_gt == "K"
		replace recodegenotype_error = . if recodegenotype_gt == "W"
		replace recodegenotype_error = . if recodegenotype_gt == "S"
		replace recodegenotype_error = . if recodegenotype_gt == "ID"
		count if recodegenotype_error == 1
		if `r(N)' != 0 {
			count if recodegenotype_error == 1
			noi di as error"#########################################################################"
			noi di as error"# > `function' ...................................... WARNING"
			noi di as error"#########################################################################"
			noi di as error"# > .......... non standard SNPs in dataset (N) `r(N)'"
			noi di as error"# > .......... exclude problematic snps using recodegenotype_error.exclude"
			noi di as error"# > .......... update-alleles using recodegenotype_error.update-alleles"
			outsheet snp if recodegenotype_error == . using recodegenotype_error.exclude, non noq replace
			keep chr snp bp a1 a2 recodegenotype_gt recodegenotype_error
			keep if recodegenotype_error == 1
			gen recodegenotype_count = length(a1)
			gen recodegenotype_b1 = ""
			gen recodegenotype_b2 = ""
			sum recodegenotype_count
			forval i = 1 / `r(max)' { 
				gen recodegenotype_a1_`i' = substr(a1, `i', 1) 
				gen recodegenotype_a2_`i' = substr(a2, `i', 1) 
				replace recodegenotype_b1 = recodegenotype_b1 + recodegenotype_a1_`i' if   recodegenotype_a1_`i' != recodegenotype_a2_`i'
				replace recodegenotype_b2 = recodegenotype_b2 + recodegenotype_a2_`i' if   recodegenotype_a1_`i' != recodegenotype_a2_`i'
				drop recodegenotype_a1_`i' recodegenotype_a2_`i'
				} 
			outsheet snp a1 a2 recodegenotype_b1 recodegenotype_b2 using recodegenotype_error.update-alleles, non noq replace
			noi di as error"# > .......... exiting"
			noi di as error"#########################################################################"
			exit
			}		
		else {
			drop recodegenotype_counta1 recodegenotype_counta2 recodegenotype_error
			}
		}
	}
noi di as text"#########################################################################"
noi di as text"# Completed: $S_DATE $S_TIME"
noi di as text"#########################################################################"
end;
