*! 1.0.4 Richard Anney 24apr2019
* clean up-front functions
*! 1.0.3 Richard Anney 21jan2019

program  snp2refid
local    function snp2refid
local    `function'_version v2
syntax , ref(string asis) 

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

qui { // module 1 - map ref short name
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 1: map ref short name"
	noi di as text"#########################################################################"
	ref_path_to_short, ref(`ref')
	local `function'_ref_short ${ref_short}
	noi di as text""	
	}
qui { // module 2 - define files
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 2: define files"
	noi di as text"#########################################################################"
	qui { // check files
		noi checkfile, file(${init_root}/bim2dta/``function'_ref_short'/``function'_ref_short'_bim_noALLELE.dta)
		noi di as text"#########################################################################"
		noi di as text""	
		}
	}
qui { // module 3 - create _bim.dta files
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3: create _bim.dta files"
	noi di as text"#########################################################################"
	foreach i in ref {
		local already_created ${init_root}/bim2dta/``function'_`i'_short'/``function'_`i'_short'_bim_noALLELE.dta
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
			use ${init_root}/bim2dta/``function'_`i'_short'/``function'_`i'_short'_bim.dta, clear
			keep chr snp bp loc_name
			save `already_created', replace
			noi checkfile, file(`already_created')
			}		
		}
	noi di as text"# > .................. identify overlap"
	rename snp old_snp
	count
	local `function'_bimN `r(N)'
	checkloc_name
	duplicates drop loc_name, force
	merge 1:1 chr bp loc_name using `already_created'
	keep if _m == 3
	drop old_snp _m loc_name
	order snp
	}
noi di as text"#########################################################################"
noi di as text"# Completed: $S_DATE $S_TIME"
noi di as text"#########################################################################"
end;	
