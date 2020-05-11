*! 1.0.1 Richard Anney 30apr2020

program  summaryqc_update
local    function summaryqc
local    `function'_version v1
syntax ,  summaryqc(string asis) build_ref(string asis) id_ref(string asis) [png(string asis)]

noi di as text""
noi di as text"#########################################################################"
noi di as text"# > .................. `function'"
noi di as text"# > ........... author "as result"Richard Anney"
noi di as text"# > ... report bugs to "as result"anneyr@cardiff.ac.uk"
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
*
qui { // module 1 - prepare folders 
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 1: define short names"
	noi di as text"#########################################################################"
	ref_path_to_short, ref(`id_ref')
	local `function'_ref_short ${ref_short}
	ref_path_to_short, ref(`summaryqc')
	local `function'_summaryqc_short ${ref_short}
	noi di as text""	
	}
qui { // module 2 - define files
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 2: define files"
	noi di as text"#########################################################################"
	noi di as text"# > ......... locating files to process"
	noi checkfile, file(${init_root}/summaryqc/``function'_summaryqc_short'/``function'_summaryqc_short'-summaryqc.dta)
	noi checkfile, file(${init_root}/bim2dta/``function'_ref_short'/``function'_ref_short'_bim_noALLELE.dta)
	noi checkfile, file(`build_ref')
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 4 - update build
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 4: update name to ``function'_ref_short'"
	noi di as text"#########################################################################"	
	use ${init_root}/summaryqc/``function'_summaryqc_short'/``function'_summaryqc_short'-summaryqc.dta, clear
	snp2build, ref(`build_ref')
	gen build = "${snp2build_build}"
	if build == "hg19 +1" {
		noi di as text"# > ... study build is "as result "${snp2build_build}"
		snp2refid, ref(`id_ref')
		}
	else {
		noi di as text"# > ... study build is "as error "${snp2build_build}"
		exit
		}
	count
	local summaryqc_Nout `r(N)'
	noi di as text"# > ..... SNPs in file "as result "`summaryqc_Nout'"
	lab var chr 			"Chromosome"
	lab var snp 			"Marker Name"
	lab var bp 				"Physical Location"
	lab var a1 				"Allele 1 (ACGT(D/I))"
	lab var a2 				"Allele 2 (ACGT(D/I))"
	lab var beta      "Beta Coefficient"
	lab var se        "Standard Error"
	lab var z         "Z-statistic"
	lab var or        "Odds Ratio"
	lab var l95       "95% Confidence Interval (lower)"
	lab var u95       "95% Confidence Interval (upper)"
	lab var p         "P-Value"
	lab var n         "Sample Size"
	compress	
	order chr bp snp a1 a2 beta se z or l95 u95 p n
	keep  chr bp snp a1 a2 beta se z or l95 u95 p n
	for var chr bp: destring X, replace
	sort chr bp
	for var chr bp: tostring X, replace
	noi di as text"# > ..... save data as "as result"${init_root}/summaryqc/``function'_summaryqc_short'/``function'_summaryqc_short'-summaryqc.dta"
	save ${init_root}/summaryqc/``function'_summaryqc_short'/``function'_summaryqc_short'-summaryqc.dta, replace
	noi di as text"#########################################################################"	
	}
qui { // module 4 - plot updated manhattan
	noi di as text""
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 5: plot updated manhattan"
	noi di as text"#########################################################################"		
	graphmanhattan, summaryqc(${init_root}/summaryqc/``function'_summaryqc_short'/``function'_summaryqc_short') png(`png')
	checkfile, file(${init_root}/graphmanhattan/``function'_summaryqc_short'/``function'_summaryqc_short'_graphmanhattan.gph)
	noi di as text""
	}
qui { // module 8 - clean folder
	files2dta, dir(`c(pwd)')
	gen keep = .
	replace keep = 1 if file == "`output'-`function'.dta"
	replace keep = 1 if file == "`output'-`function'.log"
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

	
