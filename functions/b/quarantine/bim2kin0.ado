*! 1.0.1 Richard Anney 17apr2019

program  bim2kin0
local    function bim2kin0
local    `function'_version v1
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
	foreach file in bim bed fam {
		noi checkfile, file(`bim'.`file')
		}
	noi checkfile, file(${plink})
	noi checkfile, file(${plink2})
	checktabbed
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 3 - define kinship matrix
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3: define kinship matrix"
	noi di as text"#########################################################################"		
	bim2count, bim(`bim')
	noi di as text"# > ..... SNPs in file "as result "${bim2count_snp}"
	noi di as text"# > .... indiv in file "as result "${bim2count_ind}"
	noi di as text"# > .................. define ld pruned subset"
	local already_created ${init_root}/bim2ld_subset/``function'_bim_short'/``function'_bim_short'-bim2ld_subset.extract
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		}
	else {
		noi di as text"# > .................. create "as result "``function'_bim_short'-bim2ld_subset.extract"
		noi bim2ld_subset, bim(`bim')
		noi checkfile, file(`already_created')
		}
	noi di as text"# > .................. create kinship matrix "as result "${bim2count_ind}"as text" x "as result "${bim2count_ind}"
	noi di as text"# > .................. kinship matrix filter "as result ".0221"
	!$plink2 --bfile `bim' --extract `already_created' --make-king-table --king-table-filter .0221 --out ``function'_bim_short'
	import delim using ``function'_bim_short'.kin0, clear case(lower)
	for var fid1-id2      : tostring X, replace 
	for var hethet-kinship: destring X, replace force
	replace kin = 0 if kin <0
	dropmiss, force
	noi di as text"# > .................. saving file as "as result"``function'_bim_short'_kin0.dta"
	save ``function'_bim_short'_kin0.dta, replace		
	noi di as text"#########################################################################"	
	noi di as text""	
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
	file write myfile `"# > .... plink version ${plink}"' _n
	file write myfile `"# > .... plink version ${plink2}"' _n
	file write myfile `"# > .... indiv. in bim ${bim2count_ind}"' _n
	file write myfile `"# > ....... saved data ``function'_bim_short'_kin0.dta"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > overview of `function'"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > 1. define ld pruned subset "' _n
	file write myfile `"# > 2. create relatedness matrix (report filter at .0221)"' _n
	file write myfile `"#########################################################################"' _n
	file close myfile	
	}
qui { // module 5 - clean folder
	files2dta, dir(`c(pwd)')
	gen keep = .
	replace keep = 1 if file == "``function'_bim_short'_kin0.dta"
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
