*! 1.0.1 Richard Anney 30apr2019

program  intersect_genes
local    function intersect_genes
local    `function'_version v1
syntax  , generef(string asis) [range(string asis) range_list(string asis)] 

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
	local `function'_generef_short Homo_sapiens.GRCh37.87.gtf_exon
	cd ${init_root}	
	!mkdir `function'
	cd `function'
	noi di as text"# > ..... project name "as result"``function'_generef_short'"
	!mkdir ``function'_generef_short'
	cd ``function'_generef_short'
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 2 - define files
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 2: define files"
	noi di as text"#########################################################################"
	noi checkfile, file(`generef'/Homo_sapiens.GRCh37.87.gtf_exon.dta)
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 3 - define region to plot
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3: define region to plot"
	noi di as text"#########################################################################"
	clear
	set obs 1
	gen `function'_range = "`range'"
	gen `function'_range_list = "`range_list'"
	if `function'_range != "" {
		noi di as text"# > ....... processing "as result"single range defined"
		split `function'_range,p("chr"":""..")
		local chr 	= `function'_range2[1]
		local from 	= `function'_range3[1]
		local to 		= `function'_range4[1]
		noi di as text"# > define range as"
		noi di as text"# > .............. chr "as result"`chr'"
		noi di as text"# > .......... from bp "as result"`from'"
		noi di as text"# > ............ to bp "as result"`to'"
		noi di as text"# > define genes in region "
		use `generef'/Homo_sapiens.GRCh37.87.gtf_exon.dta, clear
		gen intersect = ""
		for var start end: destring X, replace
		qui { // transcripts within region
			replace intersect = "`range'" if chr == `chr' & inrange(start,`from',`to')
			replace intersect = "`range'" if chr == `chr' & inrange(end,`from',`to')
			}
		qui { // transcripts the encompass the region
			replace intersect = "`range'" if chr == `chr' & start <`from' & end >`to'
			}
		lab var intersect "overlapping range"
		keep if intersect == "`range'"
		count
		if `r(N)' == 0 {
			set obs 1
			replace biotype == "dummy"
			}
		save Homo_sapiens.GRCh37.87.gtf_exon-with-chr`chr'_`from'_`to'_`function'.dta, replace
		noi checkfile, file(Homo_sapiens.GRCh37.87.gtf_exon-with-chr`chr'_`from'_`to'_`function'.dta)
		}
	else if `function'_range_list != "" {
		noi di as text"# > ....... processing "as result"range file defined"
		qui { // define name of file
			ref_path_to_short, ref(`range_list')	
			clear
			set obs 1
			gen name = "Homo_sapiens.GRCh37.87.gtf_exon_${ref_short}"
			split name, p(".dta")
			replace name = name1 + "_`function'.dta"
			local `function'_name = name[1]
			}
		use `range_list', clear
		count
		foreach num of num 1 / `r(N)' {
			use `range_list', clear
			keep in `num'
			split `function'_range,p("chr"":""..")
			local chr 	= `function'_range2[1]
			local from 	= `function'_range3[1]
			local to 		= `function'_range4[1]
			noi di as text"# > define range as"
			noi di as text"# > .............. chr "as result"`chr'"
			noi di as text"# > .......... from bp "as result"`from'"
			noi di as text"# > ............ to bp "as result"`to'"
			noi di as text"# > define genes in region "
			use `generef'/Homo_sapiens.GRCh37.87.gtf_exon.dta, clear
			gen intersect = ""
			for var start end: destring X, replace
			qui { // transcripts within region
				replace intersect = "chr`chr':`from'..`to'" if chr == `chr' & inrange(start,`from',`to')
				replace intersect = "chr`chr':`from'..`to'" if chr == `chr' & inrange(end,`from',`to')
				}
			qui { // transcripts the encompass the region
				replace intersect = "chr`chr':`from'..`to'" if chr == `chr' & start <`from' & end >`to'
				}
			lab var intersect "overlapping range"
			keep if intersect == "chr`chr':`from'..`to'"
			count
			if `r(N)' == 0 {
				set obs 1
				replace biotype == "dummy"
				}
			compress
			save Homo_sapiens.GRCh37.87.gtf_exon-with-chr`chr'_`from'_`to'_`function'.dta, replace
			}		
		use `range_list', clear
		split `function'_range,p("chr"":""..")
		gen script = `"append using Homo_sapiens.GRCh37.87.gtf_exon-with-chr"' + `function'_range2 + "_" + `function'_range3 + "_" + `function'_range4 + "_`function'.dta"
		replace script = "use Homo_sapiens.GRCh37.87.gtf_exon-with-chr" + `function'_range2 + "_" + `function'_range3 + "_" + `function'_range4 + "_`function'.dta, clear" in 1
		outsheet script using script.do, non noq replace
		do script.do
		erase script.do
		compress
		save ``function'_name', replace
		noi checkfile, file(``function'_name')
		}
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui cd `return'
restore
noi di as text"#########################################################################"
noi di as text"# Completed: $S_DATE $S_TIME"
noi di as text"#########################################################################"
end;
