*! 1.0.1 Richard Anney 1may2020

program  snp2ldpairs
local    function snp2ldpairs
local    `function'_version v1
syntax , ldref(string asis) snplist(string asis)

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
	ref_path_to_short, ref(`ldref')
	local `function'_ldref_short ${ref_short}
	cd ${init_root}	
	!mkdir `function'
	cd `function'
	ref_path_to_short, ref(`snplist')
	local `function'_snplist_short ${ref_short}
	noi di as text"# > ..... project name "as result"``function'_snplist_short'"
	!mkdir ``function'_snplist_short'
	cd ``function'_snplist_short'
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 2 - define files
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 1: check files"
	noi di as text"#########################################################################"
	foreach file in bim bed fam {
		noi checkfile, file(`ldref'.`file')
		}
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 3 - define ld-partners
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3: define ld-partners"
	noi di as text"#########################################################################"
	!${plink} --bfile `ldref' --r2 --ld-window 100000 --ld-window-kb 100000 --ld-snp-list `snplist'.snplist --ld-window-r2 0.8 --out ``function'_snplist_short'
	!${tabbed} ``function'_snplist_short'.ld
	import delim using ``function'_snplist_short'.ld.tabbed, clear varnames(1)
	rename (snp_a snp_b) (snp snp_proxy)
	for var chr_b bp_b: tostring X, replace
	gen position_proxy = "chr" + chr_b + ":" + bp_b
	keep snp snp_proxy position_proxy r2
	lab var snp 						"Marker Name"
	lab var snp_proxy 			"LD Partner Marker Name"
	lab var position_proxy 	"LD Partner Chromosome:Physical Location"
	lab var r2 							"Linkage Disequilibrium (R^2)"
	save ``function'_snplist_short'-`function'.dta, replace
	noi checkfile, file(``function'_snplist_short'-`function'.dta)
	noi di as text"#########################################################################"
	noi di as text""
	}
qui { // module 4 - clean folder
	files2dta, dir(`c(pwd)')
	gen keep = .
	replace keep = 1 if file == "``function'_snplist_short'-`function'.dta"
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
