*! 1.0.3 Richard Anney 17apr2019
* clean up-front functions
*! 1.0.2 Richard Anney 21dec2018
* updated to deal with non ACGT variants
* now crashes out if non-standard genotypes exist

program  bim2dta
local    function bim2dta
local    `function'_version v3
syntax , bim(string asis) 

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
qui { // module 2 - define files
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 2: define files"
	noi di as text"#########################################################################"
	noi checkfile, file(`bim'.bim)
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 3 - import file
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3: import file"
	noi di as text"#########################################################################"	
	bim2count, bim(`bim')
	local `function'_Nin ${bim2count_snp}
	noi di as text"# > ..... SNPs in file "as result "${bim2count_snp}"
	noi di as text"# > ........... import "as result "``function'_bim_short'"
	qui { // importing bim file (checking for delimit format)
		import delim  using `bim'.bim, clear
		capture confirm variable v2
		if !_rc { 
			}
		else {
			rename v1 v
			split v,p(" ")
			drop v
			}
		rename (v1 v2 v4 v5 v6) (chr snp bp a1 a2)
		noi di as text"#########################################################################"
		}
	}
qui { // module 4 - recoding genotype
	noi di as text""	
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 4: creating gt variable using "as result"recodegenotype"
	noi di as text"#########################################################################"
	keep chr bp snp a1 a2
	qui recodegenotype
	capture confirm var recodegenotype_error
	if _rc==0 {
		exit
		}
	rename recodegenotype_gt gt
	}
qui { // module 5 - rename chromosomes
	noi di as text""	
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 5: update chromosome names"
	noi di as text"#########################################################################"
	for var chr bp: tostring X, replace force
	replace chr = subinstr(chr,"chr","",.)
	replace chr = strupper(chr)	
	replace chr = "23" if chr == "X"
	replace chr = "24" if chr == "Y"
	replace chr = "25" if chr == "XY"
	replace chr = "26" if chr == "MT"
	}
qui { // module 6 - create loc_name 
	noi di as text""	
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 6: creating loc_name variable using "as result"checkloc_name"
	noi di as text"#########################################################################"
	qui checkloc_name
	}
qui { // module 7 - saving file
	noi di as text""	
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 7: saving file"
	noi di as text"#########################################################################"	
	order chr snp bp a1 a2 gt loc_name
	keep  chr snp bp a1 a2 gt loc_name
	lab var chr 			"Chromosome"
	lab var snp 			"Marker Name"
	lab var bp 				"Physical Location"
	lab var a1 				"Allele 1 (ACGT(D/I))"
	lab var a2 				"Allele 2 (ACGT(D/I))"
	lab var gt 				"Genotype (RYMKSW(ID))"
	lab var loc_name 	"Marker loc_name (chr<#>:<pos>-<gt>)"
	compress
	noi di as text"# > ..... save data as "as result"``function'_bim_short'_bim.dta"
	save ``function'_bim_short'_bim.dta, replace
	}
qui { // module 8 - write log file
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
	file write myfile `"# > ..... SNPs in file ${bim2count_snp}"' _n
	file write myfile `"# > ...... saved graph ``function'_bim_short'_bim.dta"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > overview of `function'"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > 1. import bim file "' _n
	file write myfile `"# > 2. create gt variable"' _n
	file write myfile `"# > 3. convert indels to I/D"' _n
	file write myfile `"# > 4. convert chromosomes to numeric codes"' _n
	file write myfile `"# > 5. create loc_name variable"' _n
	file write myfile `"#########################################################################"' _n
	file close myfile	
	}
qui { // module 9 - clean folder
	files2dta, dir(`c(pwd)')
	gen keep = .
	replace keep = 1 if file == "``function'_bim_short'_bim.dta"
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

