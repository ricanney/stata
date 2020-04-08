*! 1.0.2 Richard Anney 17apr2019
* clean up-front functions
*! 1.0.1 Richard Anney 21dec2018

program  checkloc_name
local    function checkloc_name
local    `function'_version v2
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

qui { // module 1 - check whether the variable loc_name exists"
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 1: check whether the variable loc_name exists"
	noi di as text"#########################################################################"
	capture confirm string var loc_name
	if !_rc {
		noi di as text"# > .......... located "as result "loc_name"
		lab var loc_name 	"Marker loc_name (chr<#>:<pos>-<gt>)"
		}
	else {
		noi di as text"# > ...... not located "as text "loc_name"
		noi di as text"# > ......... generate "as result "loc_name"
		foreach varname in chr bp {
			capture confirm string var `varname'
			if !_rc {
				noi di as text"# > .......... located "as result "`varname'"as text" (string)"
				}
			else {				
				capture confirm numeric var `varname'
				if !_rc {
					noi di as text"# > .......... located "as result "`varname'"as text" (numeric)"
					noi di as text"# > .................. convert to string"
					tostring `varname', replace
					}
				else {
					noi di as text"# > ...... not located "as text "`varname'"
					noi di as text"# > .................. exit"
					exit
					}
				}
			}
		capture confirm string var gt
		if !_rc {
			noi di as text"# > .......... located "as result "gt"
			gen _gt = gt
			replace _gt = "R" if gt == "Y"
			replace _gt = "M" if gt == "K"
			gen loc_name = "chr" + chr + ":" + bp + "-" + _gt
			drop _gt
			}
		else {		
			noi di as text"# > ...... not located "as text "gt"
			noi di as text"# > ......... generate "as result "gt"		
			capture confirm string var a1
			if !_rc {
				noi di as text"# > .......... located "as result "a1"
				capture confirm string var a2
				if !_rc {
					noi di as text"# > .......... located "as result "a2"
					recodegenotype 
					gen gt = recodegenotype_gt
					replace recodegenotype_gt = "R" if gt == "Y"
					replace recodegenotype_gt = "M" if gt == "K"
					gen loc_name = "chr" + chr + ":" + bp + "-" + recodegenotype_gt
					drop recodegenotype_gt gt					
					}
				else {
					noi di as text"# > ...... not located "as text "a2"
					noi di as text"# > .................. exit"
					exit
					}
				}
			else {
				noi di as text"# > ...... not located "as text "a1"
				noi di as text"# > .................. exit"
				exit
				}
			}
		lab var loc_name 	"Marker loc_name (chr<#>:<pos>-<gt>)"
		}
	noi di as text"#########################################################################"
	noi di as text""
	}
noi di as text"#########################################################################"
noi di as text"# Completed: $S_DATE $S_TIME"
noi di as text"#########################################################################"
end;
