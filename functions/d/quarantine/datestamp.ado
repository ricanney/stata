program datestamp
syntax 

local function datestamp
local `function'_version v1

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

preserve
local return `c(pwd)'
qui { //
	clear
	set obs 1
	gen a = "global DATE "
	gen b = "$S_DATE"
	replace b = subinstr(b, " ", "",.)
	replace b = strlower(b)	
	outsheet using datestamp.do, replace non noq
	do  datestamp.do
	erase  datestamp.do
	noi di as text"# > datestamp "as text" ...................... reporting \$DATE as" as result" $DATE"
	}
qui cd `return'
restore
noi di as text"#########################################################################"
noi di as text"# Completed: $S_DATE $S_TIME"
noi di as text"#########################################################################"
end;
