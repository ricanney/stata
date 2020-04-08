*! 1.0.1 Richard Anney 21dec2018
* update to expect a1 a2 b1 b2

program  recodestrand
global   recodestrand_version v1
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
	noi di as text"# > recodestrand ........................................ checking that alleles are defined correctly"
	foreach i in a1 a2 b1 b2 { 
		capture confirm var `i'
		if _rc==0 {
			noi di as text"# > recodestrand ...................................... `i' "as result"present"
			}
		else {
			noi di as text"# > recodestrand ...................................... `i' "as error"absent"
			exit
			}
		}
	noi di as text"#########################################################################"
	}
qui { // module 2 - define genotypes
	noi di as text""	
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 2: define genotypes"
	noi di as text"#########################################################################"	
	noi di as text"# > recodestrand ........................................ create genotype variables"
	rename (a1 a2 b1 b2)  (_a1 _a2 _b1 _b2)
	gen a1 = _a1 
	gen a2 = _a2
	recodegenotype
	rename recodegenotype_gt a_gt
	drop a1 a2
	gen a1 = _b1 
	gen a2 = _b2	
	recodegenotype
	rename recodegenotype_gt b_gt	
	drop a1 a2
	rename (_a1 _a2 _b1 _b2) (a1 a2 b1 b2)
	noi di as text"#########################################################################"	
	}
qui { // module 3 - remove incompatible genotypes
	noi di as text""	
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3: remove incompatible genotypes"
	noi di as text"#########################################################################"	
	noi di as text"# > recodestrand ................. SNPs prior to cleaning "as result  `r(N)'
	gen keep = .
	foreach i in a b {
		foreach j in R Y M K ID {
			replace keep = 1 if `i'_gt == "`j'"
			}
		}
	replace keep = . if a_gt == "R" & b_gt == "M"
	replace keep = . if a_gt == "R" & b_gt == "K"
	replace keep = . if a_gt == "Y" & b_gt == "M"
	replace keep = . if a_gt == "Y" & b_gt == "K"
	replace keep = . if a_gt == "M" & b_gt == "R"
	replace keep = . if a_gt == "M" & b_gt == "Y"
	replace keep = . if a_gt == "K" & b_gt == "R"
	replace keep = . if a_gt == "K" & b_gt == "Y"
	replace keep = . if a_gt == "ID" & b_gt != "ID"
	replace keep = . if b_gt == "ID" & a_gt != "ID"
	keep if keep == 1
	drop keep
	noi di as text"# > recodestrand .................. SNPs post to cleaning "as result  `r(N)'
	noi di as text"#########################################################################"	

	}
qui { // module 4 - identify alleles to flip
	noi di as text""	
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 4: identify alleles to flip"
	noi di as text"#########################################################################"	
	gen flip = .
	replace flip = 1 if a_gt != b_gt
	count if flip == 1
	noi di as text"# > recodestrand ........................... SNPs to flip "as result  `r(N)'

	}
noi di as text"#########################################################################"
noi di as text"# Completed: $S_DATE $S_TIME"
noi di as text"#########################################################################"
end;

