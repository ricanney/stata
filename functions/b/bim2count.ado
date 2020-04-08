*! 1.0.2 Richard Anney 16apr2019
* clean up-front functions
* allow absence of fam or bim file in routine

program  bim2count
local    function bim2count
local    `function'_version v2
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
	noi di as text"# SECTION - 2: check files"
	noi di as text"#########################################################################"
	noi checkfile, file(`bim'.bim)
	noi checkfile, file(`bim'.fam)
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 3 - count
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3: count"
	noi di as text"#########################################################################"
	qui { // snps
		capture confirm file `bim'.bim
		if !_rc {
			clear
			set obs 1
			gen os = "`c(os)'"
			if os == "Unix" { 
				!wc -l `bim'.bim  > bim.count
				import delim using bim.count, clear varnames(nonames)
				erase bim.count
				split v1,p(" ")
				destring v11, replace
				sum v11
				global bim2count_snp `r(max)'
				}
			else if os == "Windows" { 	
				import delim using  `bim'.bim, clear varnames(nonames) colrange(1:)
				count
				clear
				set obs 1
				global bim2count_snp `r(N)'
				}
			noi di as text"# > ..... SNPs in file "as result "${bim2count_snp}"
			}
		else {
			global bim2count_snp NaN
			noi di as text"# > ..... SNPs in file "as result "${bim2count_snp}"
			}
		}
	qui { // individuals
		capture confirm file `bim'.fam
		if !_rc {
			clear
			set obs 1
			gen os = "`c(os)'"
			if os == "Unix" { 
				!wc -l `bim'.fam  > fam.count
				import delim using fam.count, clear varnames(nonames)
				erase fam.count
				split v1,p(" ")
				destring v11, replace
				sum v11
				global bim2count_ind `r(max)'
				}
			else if os == "Windows" { 	
				import delim using  `bim'.fam, clear varnames(nonames) colrange(1:)
				count
				clear
				set obs 1
				global bim2count_ind `r(N)'
				}
			noi di as text"# > .... indiv in file "as result "${bim2count_ind}"
			}
		else {
			global bim2count_ind NaN
			noi di as text"# > .... indiv in file "as result "${bim2count_ind}"
			}
		}
	}
qui { // module 4 - write log file
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
	file write myfile `"# > .... indiv in file ${bim2count_ind}"' _n
	file write myfile `"#########################################################################"' _n
	file close myfile	
	}	
restore
qui cd `return'
noi di as text"#########################################################################"
noi di as text"# Completed: $S_DATE $S_TIME"
noi di as text"#########################################################################"
end;	
