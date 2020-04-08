*! 1.0.4 Richard Anney 24apr2019
* clean up-front functions
*! 1.0.3 Richard Anney 21jan2019

program  bim2refid
local    function bim2refid
local    `function'_version v4
syntax , bim(string asis) ref(string asis) 

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
	ref_path_to_short, ref(`ref')
	local `function'_bim_short ${bim_short}
	local `function'_ref_short ${ref_short}
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
		noi checkfile, file(${init_root}/bim2dta/``function'_ref_short'/``function'_ref_short'_bim.dta)
		checkfile, file(${plink})
		checktabbed
		noi di as text"#########################################################################"
		noi di as text""	
		}
	}
qui { // module 3 - create _bim.dta files
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3: create _bim.dta files"
	noi di as text"#########################################################################"
	foreach i in bim ref {
		local already_created ${init_root}/bim2dta/``function'_`i'_short'/``function'_`i'_short'_bim.dta
		capture confirm file `already_created'
		if !_rc {
			noi checkfile, file(`already_created')
			}
		else {
			noi di as text"# > .................. create "as result "``function'_`i'_short'_bim.dta"
			foreach file in bim bed fam { 
				noi checkfile, file(`i'.`file')
				}
			noi bim2dta, bim(``i'')
			noi checkfile, file(`already_created')
			}		
		}
	noi di as text"# > .................. identify overlap"
	use snp loc_name using       ${init_root}/bim2dta/``function'_bim_short'/``function'_bim_short'_bim.dta, clear
	rename snp old_snp
	count
	local `function'_bimN `r(N)'
	noi di as text"# > ..... SNPs in bim "as result"``function'_bimN'"
	merge m:1 loc_name using ${init_root}/bim2dta/``function'_ref_short'/``function'_ref_short'_bim.dta
	rename snp new_snp
	count if _m != 1
	local `function'_refN `r(N)'
	noi di as text"# > ..... SNPs in ref "as result"``function'_refN'"
	keep if _m == 3
	keep old_snp new_snp loc_name
	count
	local `function'_mergeN `r(N)'
	noi di as text"# > ... SNPs in merge "as result"``function'_mergeN'"
	qui { // check and remove duplicates
		duplicates tag new_snp, gen(tag)
		egen keep = seq(),by(new_snp tag)
		keep if keep == 1
		}
	qui { // export extract list and update-name
		outsheet old_snp using         ``function'_bim_short'-`function'.extract, non noq replace
		outsheet old_snp new_snp using ``function'_bim_short'-`function'.update-name, non noq replace
		}
	qui { // process PLINK binaries
		!$plink --bfile `bim'                  --extract     ``function'_bim_short'-`function'.extract         --make-bed --out ``function'_bim_short'
		!$plink --bfile ``function'_bim_short' --update-name ``function'_bim_short'-`function'.update-name     --make-bed --out ``function'_bim_short'-refid
		}
	qui { // confirm creation
		noi checkfile, file(``function'_bim_short'-refid.bed)
		noi checkfile, file(``function'_bim_short'-refid.bim)
		noi checkfile, file(``function'_bim_short'-refid.fam)
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
	file write myfile `"# > ... input binaries ``function'_ref_short'"' _n
	file write myfile `"# > .... reference ids ``function'_ref_short'"' _n
	file write myfile `"# > .... SNPs in input ``function'_bimN'"' _n
	file write myfile `"# > ... SNPs in output ``function'_mergeN'"' _n
	file write myfile `"# > ... saved binaries ``function'_bim_short'-refid"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > overview of `function'"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > 1. import bim file"' _n
	file write myfile `"# > 2. merge with ref using loc_name"' _n
	file write myfile `"# > 3. remove duplicate"' _n
	file write myfile `"# > 4. create new plink binaries based on overlap"' _n
	file write myfile `"# > 5. update name of new binaries"' _n
	file write myfile `"#########################################################################"' _n
	file close myfile	
	}
qui { // module 5 - clean folder
	files2dta, dir(`c(pwd)')
	gen keep = .
	replace keep = 1 if file == "``function'_bim_short'-refid.bed"
	replace keep = 1 if file == "``function'_bim_short'-refid.bim"
	replace keep = 1 if file == "``function'_bim_short'-refid.fam"
	replace keep = 1 if file == "``function'_bim_short'-refid.log"
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
