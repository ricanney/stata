*! 1.0.2 Richard Anney 24apr2019
* clean up-front functions
*! 1.0.1 Richard Anney 21dec2018

program  dir2dta
local    function dir2dta
local    `function'_version v2
syntax , dir(string asis)  

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

qui { // 
	noi di as text"# > ........ directory "as result "`dir'"
	clear
	set obs 1								
	gen folder = ""							
	save _`function'.dta,replace
	local myfiles: dir "`dir'" dirs "*" 	, respectcase				
	foreach folder of local myfiles {
		clear								
		set obs 1							
		gen folder = "`folder'" 					
		append using _`function'.dta						
		save _`function'.dta	,replace						
		}
	drop if folder == ""	
	save _`function'.dta	,replace	
	noi di as text"# > ......... saved to "as result "_`function'.dta"
	}

noi di as text"#########################################################################"
noi di as text"# Completed: $S_DATE $S_TIME"
noi di as text"#########################################################################"
end;	
