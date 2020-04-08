*! 1.0.2 Richard Anney 31mar2020

program  checktabbed
local    function checktabbed
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


preserve
local return `c(pwd)'
qui { // module 1 - check location of tabbed.pl - remove "perl "
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 1: check location of tabbed.pl"
	noi di as text"#########################################################################"
	cd ${init_root}
	clear
	set obs 1
	gen a = "${tabbed}"
	replace a = subinstr(a,"perl ","",.)
	local `function'_file = a[1]
	noi checkfile, file(``function'_file')
	noi di as text"#########################################################################"
	noi di as text""
	}
qui { // module 2 - test script
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 2: test script"
	noi di as text"#########################################################################"
	cd ${init_root}
	clear
	set obs 1
	gen a = "a b c d"
	outsheet a using `function'_test-perl-script.txt, noq replace
	!$tabbed  `function'_test-perl-script.txt
	capture confirm file "`function'_test-perl-script.txt.tabbed"
	if _rc==0 {
		noi di as text"# > ...... `function'" as result" tabbed.pl is set up correctly and working"
		}
	else {
		noi di as text"# > ...... `function'" as error" tabbed.pl is not set up correctly"
		noi di as error"# > active perl should be downloaded/installed on your computer (https://www.activestate.com/activeperl/downloads)"
		}
	noi di as text"#########################################################################"
	noi di as text""
	}
qui { // module 3 - clean up
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3: clean up temporary files"
	noi di as text"#########################################################################"
	erase `function'_test-perl-script.txt
	erase `function'_test-perl-script.txt.tabbed
	noi di as text""
	}
qui cd `return'
restore
noi di as text"#########################################################################"
noi di as text"# Completed: $S_DATE $S_TIME"
noi di as text"#########################################################################"
end;	
	
