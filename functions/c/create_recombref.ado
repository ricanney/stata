*! 1.0.1 Richard Anney 2apr2020

program  create_recombref
local    function create_recombref
local    `function'_version v1
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
qui { // module 1 - prepare folders 
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 1: prepare folders"
	noi di as text"#########################################################################"
	cd ${init_root}	
	!mkdir `function'
	cd `function'
	noi di as text"# > ..... project name "as result"genetic_map_HapMapII_GRCh37"
	!mkdir  genetic_map_HapMapII_GRCh37
	cd      genetic_map_HapMapII_GRCh37
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 2 - download and convert master files
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 2: download and convert to dta"
	noi di as text"#########################################################################"
	local already_created ${init_root}/`function'/genetic_map_HapMapII_GRCh37/genetic_map_HapMapII_GRCh37-recombination-rate.dta
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		}
	else {
			clear
			set obs 1
			gen os = "`c(os)'"
			if os == "Unix" { 
				qui { // download files
					!wget   ftp://ftp.ncbi.nlm.nih.gov/hapmap/recombination/2011-01_phaseII_B37/genetic_map_HapMapII_GRCh37.tar.gz
					!tar -xvf genetic_map_HapMapII_GRCh37.tar.gz
					}
				}
			else if os == "Windows" { 	
				qui { // download files
					!bash -c "wget   ftp://ftp.ncbi.nlm.nih.gov/hapmap/recombination/2011-01_phaseII_B37/genetic_map_HapMapII_GRCh37.tar.gz"
					!bash -c "tar -xvf genetic_map_HapMapII_GRCh37.tar.gz"
					}
				}		
			import delim using genetic_map_GRCh37_chrX_par1.txt, clear
			rename (chr pos rate) (chr bp recomb_rate)
			keep chr bp recomb_rate
			replace chr = "25"
			destring chr, replace
			save genetic_map_GRCh37_chr25.dta, replace
			import delim using genetic_map_GRCh37_chrX_par2.txt, clear
			rename (chr pos rate) (chr bp recomb_rate)
			keep chr bp recomb_rate
			replace chr = "25"
			destring chr, replace
			append using genetic_map_GRCh37_chr25.dta
			save genetic_map_GRCh37_chr25.dta, replace
			import delim using genetic_map_GRCh37_chrX.txt, clear
			rename (chr pos rate) (chr bp recomb_rate)
			keep chr bp recomb_rate
			replace chr = "23"
			destring chr, replace
			save genetic_map_GRCh37_chr23.dta, replace	
			foreach chr of num 1/22 {
				import delim using genetic_map_GRCh37_chr`chr'.txt, clear
				rename (chr pos rate) (chr bp recomb_rate)
				keep chr bp recomb_rate
				replace chr = "`chr'"
				destring chr, replace
				save genetic_map_GRCh37_chr`chr'.dta, replace
				}
			use genetic_map_GRCh37_chr25.dta, clear
			foreach chr of num 1/23  {
				append using genetic_map_GRCh37_chr`chr'.dta
				}
			sort chr bp
			save genetic_map_HapMapII_GRCh37-recombination-rate.dta, replace
			
			noi checkfile, file(`already_created')
			}
	noi di as text"#########################################################################"
	noi di as text""			
	}
qui { // module 3 - clean folder
	files2dta, dir(`c(pwd)')
	gen keep = .
	replace keep = 1 if files == "genetic_map_HapMapII_GRCh37-recombination-rate.dta"
	drop if keep == 1
	gen script = "erase " + files
	outsheet script using temp.do, non noq replace 
	do temp.do
	erase temp.do
	clear 
	}
qui { //
	cd `return'
	}
restore
noi di as text"#########################################################################"
noi di as text"# Completed: $S_DATE $S_TIME"
noi di as text"#########################################################################"
end;
