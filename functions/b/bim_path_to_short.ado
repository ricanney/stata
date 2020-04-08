*! 1.0.1 Richard Anney 15apr2019
*! 1.0.2 Richard Anney 30mar20202
* removed the creation of temp.do - *use global bim_short 	= bim1[1]* instead

program  bim_path_to_short
local    function bim_path_to_short
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
qui { // module 1 - define short name from path
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 1: prepare folders"
	noi di as text"#########################################################################"
	
	
	clear
	set obs 1
	gen bim = ""
	replace bim =  "`bim'" in 1
	replace bim = subinstr(bim,"\","%%",.)
	replace bim = subinstr(bim,"/","%%",.)
	split bim, p("%%")
	gen bim999 = ""
	for var bim1-bim999: replace bim1 = X if X !=""
	keep bim1
	replace bim1 = subinstr(bim1,"_bim.dta","",.)
	compress
	global bim_short 	= bim1[1]
	}
di as text"# > ............. file ${bim_short}"
restore
end; 

