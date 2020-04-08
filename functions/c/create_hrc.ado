*! 1.0.1 Richard Anney 3apr2020

program  create_hrc
local    function create_hrc
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
	noi di as text"# > ..... project name "as result"HRC.r1-1.GRCh37.wgs.mac5"
	!mkdir  HRC.r1-1.GRCh37.wgs.mac5
	cd      HRC.r1-1.GRCh37.wgs.mac5
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 2/6 - create all-hrc-1.1-chrall-mac5_bim.dta
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3/6: create all-hrc-1.1-chrall-mac5_bim.dta"
	noi di as text"#########################################################################"
	local main_output ${init_root}/bim2dta/all-hrc-1.1-chrall-mac5/all-hrc-1.1-chrall-mac5_bim.dta
	capture confirm file `main_output'
	if !_rc {
		noi checkfile, file(${init_root}/bim2dta/all-hrc-1.1-chrall-mac5/all-hrc-1.1-chrall-mac5_bim.dta)
		}
	else {	
		qui { // module 2 - download vcf
			noi di as text"#########################################################################"
			noi di as text"# SECTION - 2: download vcf and convert to plink binaries"
			noi di as text"#########################################################################"
			qui { // download HRC.r1-1.GRCh37.wgs.mac5.sites.sites.core
				noi di as text"# > ......... download "as result"HRC.r1-1.GRCh37.wgs.mac5"
				local already_created ${init_root}/`function'/HRC.r1-1.GRCh37.wgs.mac5/HRC.r1-1.GRCh37.wgs.mac5.sites.core
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
							!wget   ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.vcf.gz
							!zcat   HRC.r1-1.GRCh37.wgs.mac5.sites.vcf.gz | cut -f 1-5,8 > HRC.r1-1.GRCh37.wgs.mac5.sites.core
							}
						}
					else if os == "Windows" { 	
						qui { // download files
							!bash -c "wget   ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.vcf.gz"
							!bash -c "zcat    HRC.r1-1.GRCh37.wgs.mac5.sites.vcf.gz | cut -f 1-5,8 > HRC.r1-1.GRCh37.wgs.mac5.sites.core"	
							}
						}
					}
				}
			noi di as text"#########################################################################"
			noi di as text""	
			}
		qui { // module 3 - convert core to core by chr
			noi di as text"#########################################################################"
			noi di as text"# SECTION - 3: convert core to core by chr"
			noi di as text"#########################################################################"	
			foreach chr of num 1/22 {
				global chr `chr'
				local already_created ${init_root}/`function'/HRC.r1-1.GRCh37.wgs.mac5/HRC.r1-1.GRCh37.wgs.mac5.sites.core.chr${chr}
				capture confirm file `already_created'
				if !_rc {
					noi checkfile, file(`already_created')
					}
				else {	
					clear 
					set obs 1
					gen os = "`c(os)'"
					if os == "Unix" { 
						!cat cat HRC.r1-1.GRCh37.wgs.mac5.sites.core | grep -w '^${chr}' > HRC.r1-1.GRCh37.wgs.mac5.sites.core.chr${chr}
						}
					else if os == "Windows" { 	
						!bash -c "cat HRC.r1-1.GRCh37.wgs.mac5.sites.core | grep -w '^${chr}' > HRC.r1-1.GRCh37.wgs.mac5.sites.core.chr${chr}"
						}
					noi checkfile, file(`already_created')
					}
				}
			global chr 23
			local already_created ${init_root}/`function'/HRC.r1-1.GRCh37.wgs.mac5/HRC.r1-1.GRCh37.wgs.mac5.sites.core.chr${chr}
			capture confirm file `already_created'
			if !_rc {
					noi checkfile, file(`already_created')
					}
			else {	
				clear 
				set obs 1
				gen os = "`c(os)'"
				if os == "Unix" { 
					!cat cat HRC.r1-1.GRCh37.wgs.mac5.sites.core | grep -w '^X' > HRC.r1-1.GRCh37.wgs.mac5.sites.core.chr${chr}
					}
				else if os == "Windows" { 	
					!bash -c "cat HRC.r1-1.GRCh37.wgs.mac5.sites.core | grep -w '^X' > HRC.r1-1.GRCh37.wgs.mac5.sites.core.chr${chr}"
					}
				noi checkfile, file(`already_created')
				}		
			noi di as text"#########################################################################"
			noi di as text""
			}			
		qui { // module 4 - process core 
			noi di as text"#########################################################################"
			noi di as text"# SECTION - 4: process core to dta"
			noi di as text"#########################################################################"	
			foreach chr of num 1/23 {
				global chr `chr'
				local already_created ${init_root}/`function'/HRC.r1-1.GRCh37.wgs.mac5/HRC.r1-1.GRCh37.wgs.mac5.sites.core.chr`chr'.dta
				capture confirm file `already_created'
				if !_rc {
					noi checkfile, file(`already_created')
					}
				else {	
					import delim using HRC.r1-1.GRCh37.wgs.mac5.sites.core.chr`chr', clear
					rename (v1 - v5) (chr bp snp a2 a1)
					keep chr bp snp a2 a1 v6
					tostring chr, replace
					replace chr = "23" if chr == "X"
					split v6,p("AF=")
					split v62,p(";")
					rename v621 af
					keep chr bp snp a1 a2 af
					replace af = subinstr(af,"AF=","",.)
					destring af, replace
					checkloc_name
					replace snp = loc_name if snp == "."
					recodegenotype
					rename recode gt
					lab var chr 			"Chromosome"
					lab var snp 			"Marker Name"
					lab var bp 				"Physical Location"
					lab var a1 				"Allele 1 (ACGT(D/I))"
					lab var a2 				"Allele 2 (ACGT(D/I))"
					lab var gt 				"Genotype (RYMKSW(ID))"
					lab var af				"Allele Frequency (a1)" 
					lab var loc_name 	"Marker loc_name (chr<#>:<pos>-<gt>)"
					for var chr bp: tostring X, replace
					save HRC.r1-1.GRCh37.wgs.mac5.sites.core.chr`chr'.dta, replace
					noi checkfile, file(`already_created')
					}
				}
			noi di as text"#########################################################################"
			noi di as text""
			}	
		qui { // module 5 - combine core 
			noi di as text"#########################################################################"
			noi di as text"# SECTION - 5: combine core dta"
			noi di as text"#########################################################################"	
			local already_created ${init_root}/`function'/HRC.r1-1.GRCh37.wgs.mac5/HRC.r1-1.GRCh37.wgs.mac5.sites.core.chrall.dta
			capture confirm file `already_created'
			if !_rc {
				noi checkfile, file(`already_created')
				}
			else {	
				use HRC.r1-1.GRCh37.wgs.mac5.sites.core.chr23.dta, clear
				foreach chr of num 1/22 {
					append using HRC.r1-1.GRCh37.wgs.mac5.sites.core.chr`chr'.dta
					}
				for var chr bp: destring X, replace
				sort chr bp
				for var chr bp: tostring X, replace
				save HRC.r1-1.GRCh37.wgs.mac5.sites.core.chrall.dta, replace
				qui { // to limit conflicts - convert duplicate rsid to loc_name
					use HRC.r1-1.GRCh37.wgs.mac5.sites.core.chrall.dta, clear
					egen x = seq(),by(snp)
					keep if x != 1 
					keep snp 
					duplicates drop 
					merge 1:m snp  using HRC.r1-1.GRCh37.wgs.mac5.sites.core.chrall.dta
					replace snp = loc_name if _m == 3
					drop _m
					}
				for var chr bp: destring X, replace
				sort chr bp
				for var chr bp: tostring X, replace
				save HRC.r1-1.GRCh37.wgs.mac5.sites.core.chrall.dta, replace
				noi checkfile, file(`already_created')
				gen cm = "0"
				keep chr snp cm bp a1 a2
				outsheet chr snp cm bp a1 a2 using all-hrc-1.1-chrall-mac5.bim, non noq replace
				bim2dta, bim(${init_root}/`function'/HRC.r1-1.GRCh37.wgs.mac5/all-hrc-1.1-chrall-mac5)
				noi checkfile, file(${init_root}/bim2dta/all-hrc-1.1-chrall-mac5/all-hrc-1.1-chrall-mac5_bim.dta)
				noi di as text"#########################################################################"
				noi di as text""
				}
			}
		}
	}
qui { // module 7 - create bim_noALLELE.dta files
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 6: create bim_noALLELE.dta files"
	noi di as text"#########################################################################"	
	local already_created ${init_root}/bim2dta/all-hrc-1.1-chrall-mac5/all-hrc-1.1-chrall-mac5_bim_noALLELE.dta
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		}
	else {	
		use ${init_root}/bim2dta/all-hrc-1.1-chrall-mac5/all-hrc-1.1-chrall-mac5_bim.dta, clear
		keep  chr snp bp loc_name
		save ${init_root}/bim2dta/all-hrc-1.1-chrall-mac5/all-hrc-1.1-chrall-mac5_bim_noALLELE.dta, replace		
		noi checkfile, file(`already_created')
		}
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // clean-up
	files2dta, dir(`c(pwd)')
	gen keep = .
	replace keep = 1 if files == "HRC.r1-1.GRCh37.wgs.mac5.sites.core.chrall.dta"
	replace keep = 1 if files == "all-hrc-1.1-chrall-mac5.bim"
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
noi di as text"#########################################################################"
noi di as text"# Completed: $S_DATE $S_TIME"
noi di as text"#########################################################################"
end;
	
