*! 1.0.2 Richard Anney 24apr2019
*! 1.0.1 Richard Anney 18jan2019
* - update to include bim name in output
* - update to save files in bim2eigenvec file and source

program  bim2eigenvec
local    function bim2eigenvec
local    `function'_version v2
syntax , bim(string asis) [pc(real 10)]

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
	noi bim_path_to_short, bim(`bim')
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
	qui { // check files
		foreach file in bim bed fam { 
			noi checkfile, file(`bim'.`file')
			}
		checkfile, file(${plink})
		checkfile, file(${plink2})
		checktabbed
		noi di as text"#########################################################################"
		}
	}
qui { // module 3 - ld-prune
	noi di as text""	
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3: ld-prune input binaries with additional" as result" bim2ldexclude"
	noi di as text"#########################################################################"
	local already_created ${init_root}/bim2ldexclude/``function'_bim_short'/``function'_bim_short'-bim2ldexclude.exclude
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		}
	else {
		noi di as text"# > .................. create "as result "``function'_bim_short'-bim2ldexclude.exclude"
		noi bim2ldexclude, bim(`bim')
		noi checkfile, file(`already_created')
		}
	noi di as text"# > .................. ld-prune using PLINK"
	!$plink --bfile `bim' --exclude `already_created' --indep-pairwise 1000 5 0.2 --out ``function'_bim_short'
	!$plink --bfile `bim' --extract ``function'_bim_short'.prune.in --make-bed    --out ``function'_bim_short'-`function'
	noi di as text"#########################################################################"	
	}
qui { // module 4 - calculate PCA
	noi di as text""	
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 4: defining principle component using " as result" PLINK2"
	noi di as text"#########################################################################"	
	bim2count, bim(`bim')
	noi di as text"# > .................. define PC using PLINK2"
	!$plink2 --bfile ``function'_bim_short'-`function' --pca `pc' --out ``function'_bim_short'-`function'
	noi di as text"# > .................. processing eigenvec files"
	!$tabbed ``function'_bim_short'-`function'.eigenvec
	import delim using ``function'_bim_short'-`function'.eigenvec.tabbed, clear
	dropmiss, force
	for var fid iid: tostring X,replace
	count
	local `function'_eigenN `r(N)'
	noi di as text"# > .................. eigenvec saved to "as result"``function'_bim_short'_eigenvec.dta"
		lab var fid 			"Family Identifier"
		lab var iid				"Individual Identifier"
		foreach num of num 1/`pc' {
			lab var pc`num'		"Ancestry Principal Component `num'"
			}
	save ``function'_bim_short'_eigenvec.dta,replace
	noi di as text"# > .................. processing eigenval files"
	!$tabbed ``function'_bim_short'-`function'.eigenval
	import delim using ``function'_bim_short'-`function'.eigenval.tabbed, clear
	gen pc = _n
	ren v1 eigenval
	dropmiss, force
	order pc
	noi di as text"# > .................. eigenval saved to "as result"``function'_bim_short'_eigenval.dta"
	lab var pc 				"Principal Component"
	lab var eigenval 	"Eigenvalue"
	save ``function'_bim_short'_eigenval.dta,replace
	}
qui { // module 5 - write log file
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
	file write myfile `"# > ............ N PCs `pc'"' _n
	file write myfile `"# > .... indiv. in bim ${bim2count_ind}"' _n
	file write myfile `"# > indiv. in eigenvec ``function'_eigenN'"' _n
	file write myfile `"# > ....... saved data ``function'_bim_short'_eigenvec.dta"' _n
	file write myfile `"# > ....... saved data ``function'_bim_short'_eigenval.dta"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > overview of `function'"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > 1. define ld pruned subset "' _n
	file write myfile `"# > 2. calculate pc"' _n
	file write myfile `"# > 3. process files to eigenvec and eigenval dta files"' _n
	file write myfile `"#########################################################################"' _n
	file close myfile	
	}
qui { // module 8 - clean folder
	files2dta, dir(`c(pwd)')
	gen keep = .
	replace keep = 1 if file == "``function'_bim_short'_eigenvec.dta"
	replace keep = 1 if file == "``function'_bim_short'_eigenval.dta"
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


	
