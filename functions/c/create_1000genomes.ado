*! 1.0.1 Richard Anney 31mar2020

program  create_1000genomes
local    function create_1000genomes
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
	noi di as text"# > ..... project name "as result"ftp.1000genomes.ebi.ac.uk"
	!mkdir  ftp.1000genomes.ebi.ac.uk
	cd      ftp.1000genomes.ebi.ac.uk
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 2 - create all-1000g-phase3-chrall-mac5-v2.population
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 2: create all-1000g-phase3-chrall-mac5-v2.population file"
	noi di as text"#########################################################################"
	local already_created ${init_root}/`function'/ftp.1000genomes.ebi.ac.uk/all-1000g-phase3-chrall-mac5-v2.population
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		}
	else {
		local already_created ${init_root}/`function'/ftp.1000genomes.ebi.ac.uk/integrated_call_samples_v3.20130502.ALL.panel
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
					!wget   http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
					}
				}
			else if os == "Windows" { 	
				qui { // download files
					!bash -c "wget   http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
					}
				}
			import delim using integrated_call_samples_v3.20130502.ALL.panel, clear varnames(1)
			rename (sample gender) (iid sex)
			gen fid = iid
			replace sex = "1" if sex == "male"
			replace sex = "2" if sex == "female"
			destring sex, replace
			lab var sex				"Sex"
			lab define sex   1"Male" 2"Female" -9"Unknown"
			lab val sex sex
			lab var pop 			"Population Code"
			lab var super_pop	"Super-Population Code"
			lab var fid 			"Family Identifier"
			lab var iid				"Individual Identifier"	
			keep  fid iid pop super_pop sex
			order fid iid pop super_pop sex
			outsheet fid iid super using  ${init_root}/`function'/ftp.1000genomes.ebi.ac.uk/all-1000g-phase3-chrall-mac5-v2.super-population, non noq replace
			outsheet fid iid pop using  ${init_root}/`function'/ftp.1000genomes.ebi.ac.uk/all-1000g-phase3-chrall-mac5-v2.population, non noq replace
			outsheet fid iid sex using  ${init_root}/`function'/ftp.1000genomes.ebi.ac.uk/all-1000g-phase3-chrall-mac5-v2.update-sex, non noq nol replace
			save ${init_root}/`function'/ftp.1000genomes.ebi.ac.uk/all-1000g-phase3-chrall-mac5-v2-panel.dta, replace
			noi checkfile, file(`already_created')
			}
		}
	noi di as text"#########################################################################"
	noi di as text""
	}
qui { // module 3/6 - create all-1000g-phase3-chrall-mac5-v2 plink binaries
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3/6: create all-1000g-phase3-chrall-mac5-v2 plink binaries"
	noi di as text"#########################################################################"
	local main_output ${init_root}/`function'/ftp.1000genomes.ebi.ac.uk/all-1000g-phase3-chrall-mac5-v2.bed
	capture confirm file `main_output'
	if !_rc {
		noi checkfile, file(${init_root}/`function'/ftp.1000genomes.ebi.ac.uk/all-1000g-phase3-chrall-mac5-v2.bed)
		noi checkfile, file(${init_root}/`function'/ftp.1000genomes.ebi.ac.uk/all-1000g-phase3-chrall-mac5-v2.bim)
		noi checkfile, file(${init_root}/`function'/ftp.1000genomes.ebi.ac.uk/all-1000g-phase3-chrall-mac5-v2.fam)
		}
	else {
		noi di as text""	
		qui { // module 3 - download vcf
		noi di as text"#########################################################################"
		noi di as text"# SECTION - 3: download vcf and convert to plink binaries"
		noi di as text"#########################################################################"
		qui { // download ALL.chr`chr'.phase3_shapeit2_mvncall_integrated_v5a.20130502*
			foreach chr of num 22/1 {
				noi di as text"# > ......... download "as result"ALL.chr`chr'.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
				local already_created ${init_root}/`function'/ftp.1000genomes.ebi.ac.uk/ALL.chr`chr'.phase3_shapeit2_mvncall_integrated_v5a.20130502-mac5.bed
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
							!wget   http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr`chr'.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
							!${plink} --vcf ALL.chr`chr'.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --mac 5 --make-bed --out ALL.chr`chr'.phase3_shapeit2_mvncall_integrated_v5a.20130502
							}
						}
					else if os == "Windows" { 	
						qui { // download files
							!bash -c "wget   http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr`chr'.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
							!${plink} --vcf ALL.chr`chr'.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --mac 5 --make-bed --out ALL.chr`chr'.phase3_shapeit2_mvncall_integrated_v5a.20130502-mac5
							}
						}		
					}
				}
			foreach chr in X {
				noi di as text"# > ......... download "as result"http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz"
				local already_created ${init_root}/`function'/ftp.1000genomes.ebi.ac.uk/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502-mac5.bed
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
							!wget   http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz
							!${plink} --vcf ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz --mac 5 --make-bed --out ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502-mac5
							}
						}
					else if os == "Windows" { 	
						qui { // download files
							!bash -c "wget   http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz"
							!${plink} --vcf ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz --mac 5 --make-bed --out ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502-mac5
							}
						}		
					}
				}
			}
		noi di as text"#########################################################################"
		noi di as text""	
		}
		qui { // module 4 - create bim2dta files > update names using hrc (autosomes)
			noi di as text"#########################################################################"
			noi di as text"# SECTION - 4: create bim2dta files > update names using hrc (autosomes)"
			noi di as text"#########################################################################"	
			foreach chr of num 22/1 {
				local already_created ${init_root}/`function'/ftp.1000genomes.ebi.ac.uk/ALL.chr`chr'.phase3_shapeit2_mvncall_integrated_v5a.20130502-mac5-v2.bed
				capture confirm file `already_created'
				if !_rc {
					noi checkfile, file(`already_created')
					}
				else {	
					import delim using 	${init_root}/`function'/ftp.1000genomes.ebi.ac.uk/ALL.chr`chr'.phase3_shapeit2_mvncall_integrated_v5a.20130502-mac5.bim, clear
					rename (v1 v2 v3 v4 v5 v6) (chr snp cm bp a1 a2)
					drop if snp == "."
					qui { // drop duplicate snp names
						duplicates tag snp, gen(dup)
						keep if dup == 0
						drop dup
						}
					qui { // drop structural variants
						split a1,p("<")
						split a2,p("<")
						drop if a11 == ""
						drop if a21 == ""
						keep chr snp cm bp a1 a2
						}
					qui { // update indels
						gen recodegenotype_counta1 = length(a1)
						gen recodegenotype_counta2 = length(a2)
						replace a1 = "I" if recodegenotype_counta1 > recodegenotype_counta2
						replace a2 = "I" if recodegenotype_counta2 > recodegenotype_counta1
						replace a1 = "D" if a2 == "I"
						replace a2 = "D" if a1 == "I"
						compress
						keep chr snp cm bp a1 a2
						}
					qui { // update common length non-standard snps
						gen recodegenotype_count = length(a1)
						gen recodegenotype_b1 = ""
						gen recodegenotype_b2 = ""
						sum recodegenotype_count
						forval i = 1 / `r(max)' { 
							gen recodegenotype_a1_`i' = substr(a1, `i', 1) 
							gen recodegenotype_a2_`i' = substr(a2, `i', 1) 
							replace recodegenotype_b1 = recodegenotype_b1 + recodegenotype_a1_`i' if   recodegenotype_a1_`i' != recodegenotype_a2_`i'
							replace recodegenotype_b2 = recodegenotype_b2 + recodegenotype_a2_`i' if   recodegenotype_a1_`i' != recodegenotype_a2_`i'
							drop recodegenotype_a1_`i' recodegenotype_a2_`i'
							} 
						rename (a1 a2) (original_a1 original_a2)
						rename (recodegenotype_b1 recodegenotype_b2) (a1 a2)
						keep chr snp cm bp a1 a2 original_a1 original_a2
						}
					qui { // create loc_name 
						recodegenotype
						rename recode gt
						checkloc_name
						}
					qui { // remove duplicate loc_names
						rename snp original_snp
						egen x = seq(),by(loc_name)
						keep if x == 1
						drop x
						}
					qui { // limit to hrc
						merge 1:1 loc_name using ${init_root}/bim2dta/all-hrc-1.1-chrall-mac5/all-hrc-1.1-chrall-mac5_bim_noALLELE.dta
						keep if _m == 3
						drop _m
						}
					qui { // update plink binaries
						outsheet original_snp     													using ${init_root}/`function'/ftp.1000genomes.ebi.ac.uk/ALL.chr`chr'.phase3_shapeit2_mvncall_integrated_v5a.20130502-mac5.extract, non noq replace
						outsheet original_snp original_a1 original_a2 a1 a2 using ${init_root}/`function'/ftp.1000genomes.ebi.ac.uk/ALL.chr`chr'.phase3_shapeit2_mvncall_integrated_v5a.20130502-mac5.update-alleles, non noq replace
						outsheet original_snp snp 													using ${init_root}/`function'/ftp.1000genomes.ebi.ac.uk/ALL.chr`chr'.phase3_shapeit2_mvncall_integrated_v5a.20130502-mac5.update-name, non noq replace
						!${plink} --bfile 					${init_root}/`function'/ftp.1000genomes.ebi.ac.uk/ALL.chr`chr'.phase3_shapeit2_mvncall_integrated_v5a.20130502-mac5 ///
											--extract					${init_root}/`function'/ftp.1000genomes.ebi.ac.uk/ALL.chr`chr'.phase3_shapeit2_mvncall_integrated_v5a.20130502-mac5.extract ///
											--make-bed				///
											--out							${init_root}/`function'/ftp.1000genomes.ebi.ac.uk/ALL.chr`chr'.phase3_shapeit2_mvncall_integrated_v5a.20130502-mac5-extracted
						!${plink} --bfile 					${init_root}/`function'/ftp.1000genomes.ebi.ac.uk/ALL.chr`chr'.phase3_shapeit2_mvncall_integrated_v5a.20130502-mac5-extracted ///
											--update-alleles	${init_root}/`function'/ftp.1000genomes.ebi.ac.uk/ALL.chr`chr'.phase3_shapeit2_mvncall_integrated_v5a.20130502-mac5.update-alleles ///
											--make-bed				///
											--out							${init_root}/`function'/ftp.1000genomes.ebi.ac.uk/ALL.chr`chr'.phase3_shapeit2_mvncall_integrated_v5a.20130502-mac5-updated-alleles								
						!${plink} --bfile 					${init_root}/`function'/ftp.1000genomes.ebi.ac.uk/ALL.chr`chr'.phase3_shapeit2_mvncall_integrated_v5a.20130502-mac5-updated-alleles	 ///
											--update-name	${init_root}/`function'/ftp.1000genomes.ebi.ac.uk/ALL.chr`chr'.phase3_shapeit2_mvncall_integrated_v5a.20130502-mac5.update-name ///
											--make-bed				///
											--out							${init_root}/`function'/ftp.1000genomes.ebi.ac.uk/ALL.chr`chr'.phase3_shapeit2_mvncall_integrated_v5a.20130502-mac5-v2
						noi checkfile, file(`already_created')
						}
					}
				}
			noi di as text"#########################################################################"
			noi di as text""	
			}
		qui { // module 5 - create bim2dta files > update names using hrc (sex-chromosome)
			noi di as text"#########################################################################"
			noi di as text"# SECTION - 5: create bim2dta files > update names using hrc (sex-chromosome)"
			noi di as text"#########################################################################"	
			foreach chr of num 23 {
				local already_created ${init_root}/`function'/ftp.1000genomes.ebi.ac.uk/ALL.chr`chr'.phase3_shapeit2_mvncall_integrated_v1b.20130502-mac5-v2.bed
				capture confirm file `already_created'
				if !_rc {
					noi checkfile, file(`already_created')
					}
				else {	
					import delim using 	${init_root}/`function'/ftp.1000genomes.ebi.ac.uk/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502-mac5.bim, clear
					rename (v1 v2 v3 v4 v5 v6) (chr snp cm bp a1 a2)
					drop if snp == "."
					qui { // drop duplicate snp names
						duplicates tag snp, gen(dup)
						keep if dup == 0
						drop dup
						}
					qui { // drop structural variants
						split a1,p("<")
						split a2,p("<")
						drop if a11 == ""
						drop if a21 == ""
						keep chr snp cm bp a1 a2
						}
					qui { // update indels
						gen recodegenotype_counta1 = length(a1)
						gen recodegenotype_counta2 = length(a2)
						replace a1 = "I" if recodegenotype_counta1 > recodegenotype_counta2
						replace a2 = "I" if recodegenotype_counta2 > recodegenotype_counta1
						replace a1 = "D" if a2 == "I"
						replace a2 = "D" if a1 == "I"
						compress
						keep chr snp cm bp a1 a2
						}
					qui { // update common length non-standard snps
						gen recodegenotype_count = length(a1)
						gen recodegenotype_b1 = ""
						gen recodegenotype_b2 = ""
						sum recodegenotype_count
						forval i = 1 / `r(max)' { 
							gen recodegenotype_a1_`i' = substr(a1, `i', 1) 
							gen recodegenotype_a2_`i' = substr(a2, `i', 1) 
							replace recodegenotype_b1 = recodegenotype_b1 + recodegenotype_a1_`i' if   recodegenotype_a1_`i' != recodegenotype_a2_`i'
							replace recodegenotype_b2 = recodegenotype_b2 + recodegenotype_a2_`i' if   recodegenotype_a1_`i' != recodegenotype_a2_`i'
							drop recodegenotype_a1_`i' recodegenotype_a2_`i'
							} 
						rename (a1 a2) (original_a1 original_a2)
						rename (recodegenotype_b1 recodegenotype_b2) (a1 a2)
						keep chr snp cm bp a1 a2 original_a1 original_a2
						}
					qui { // create loc_name 
						recodegenotype
						rename recode gt
						checkloc_name
						}
					qui { // remove duplicate loc_names
						rename snp original_snp
						egen x = seq(),by(loc_name)
						keep if x == 1
						drop x
						}
					qui { // limit to hrc
						merge 1:1 loc_name using ${init_root}/bim2dta/all-hrc-1.1-chrall-mac5/all-hrc-1.1-chrall-mac5_bim_noALLELE.dta
						keep if _m == 3
						drop _m
						}
					qui { // update plink binaries
						outsheet original_snp     													using ${init_root}/`function'/ftp.1000genomes.ebi.ac.uk/ALL.chr`chr'.phase3_shapeit2_mvncall_integrated_v1b.20130502-mac5.extract, non noq replace
						outsheet original_snp original_a1 original_a2 a1 a2 using ${init_root}/`function'/ftp.1000genomes.ebi.ac.uk/ALL.chr`chr'.phase3_shapeit2_mvncall_integrated_v1b.20130502-mac5.update-alleles, non noq replace
						outsheet original_snp snp 													using ${init_root}/`function'/ftp.1000genomes.ebi.ac.uk/ALL.chr`chr'.phase3_shapeit2_mvncall_integrated_v1b.20130502-mac5.update-name, non noq replace
						!${plink} --bfile 					${init_root}/`function'/ftp.1000genomes.ebi.ac.uk/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502-mac5 ///
											--extract					${init_root}/`function'/ftp.1000genomes.ebi.ac.uk/ALL.chr`chr'.phase3_shapeit2_mvncall_integrated_v1b.20130502-mac5.extract ///
											--make-bed				///
											--out							${init_root}/`function'/ftp.1000genomes.ebi.ac.uk/ALL.chr`chr'.phase3_shapeit2_mvncall_integrated_v1b.20130502-mac5-extracted
						!${plink} --bfile 					${init_root}/`function'/ftp.1000genomes.ebi.ac.uk/ALL.chr`chr'.phase3_shapeit2_mvncall_integrated_v1b.20130502-mac5-extracted ///
											--update-alleles	${init_root}/`function'/ftp.1000genomes.ebi.ac.uk/ALL.chr`chr'.phase3_shapeit2_mvncall_integrated_v1b.20130502-mac5.update-alleles ///
											--make-bed				///
											--out							${init_root}/`function'/ftp.1000genomes.ebi.ac.uk/ALL.chr`chr'.phase3_shapeit2_mvncall_integrated_v1b.20130502-mac5-updated-alleles								
						!${plink} --bfile 					${init_root}/`function'/ftp.1000genomes.ebi.ac.uk/ALL.chr`chr'.phase3_shapeit2_mvncall_integrated_v1b.20130502-mac5-updated-alleles	 ///
											--update-name	${init_root}/`function'/ftp.1000genomes.ebi.ac.uk/ALL.chr`chr'.phase3_shapeit2_mvncall_integrated_v1b.20130502-mac5.update-name ///
											--make-bed				///
											--out							${init_root}/`function'/ftp.1000genomes.ebi.ac.uk/ALL.chr`chr'.phase3_shapeit2_mvncall_integrated_v1b.20130502-mac5-v2
						noi checkfile, file(`already_created')
						}
					}
				}
			noi di as text"#########################################################################"
			noi di as text""	
			}
		qui { // module 6 - check for duplicate 
		noi di as text"#########################################################################"
		noi di as text"# SECTION - 6: create all-1000g-phase3-chrall-mac5-v2"
		noi di as text"#########################################################################"	
			qui { // check for snpid duplication
			*	bim2dta, bim(${init_root}/`function'/ftp.1000genomes.ebi.ac.uk/ALL.chr`chr'.phase3_shapeit2_mvncall_integrated_v1b.20130502-mac5-v2)
			*	foreach chr of num 1/22 {
			*		bim2dta, bim(${init_root}/`function'/ftp.1000genomes.ebi.ac.uk/ALL.chr`chr'.phase3_shapeit2_mvncall_integrated_v5a.20130502-mac5-v2)
			*		}
			*	use ${init_root}/bim2dta/ALL.chr23.phase3_shapeit2_mvncall_integrated_v1b.20130502-mac5-v2/ALL.chr23.phase3_shapeit2_mvncall_integrated_v1b.20130502-mac5-v2_bim.dta, clear
			*	keep snp loc_name
			*	foreach chr of num 1/22 {
			*		append using ${init_root}/bim2dta/ALL.chr`chr'.phase3_shapeit2_mvncall_integrated_v5a.20130502-mac5-v2/ALL.chr`chr'.phase3_shapeit2_mvncall_integrated_v5a.20130502-mac5-v2_bim.dta
			*		keep snp loc_name
			*		}
				}
			clear
			set obs 22
			gen chr = _n
			tostring chr, replace
			gen file = "ALL.chr" + chr + ".phase3_shapeit2_mvncall_integrated_v5a.20130502-mac5-v2"
			foreach file in bed bim fam {
					gen `file' = file + ".`file'"
					}
			outsheet bed bim fam using ALL.chrall.phase3_shapeit2_mvncall_integrated_v5a.20130502-mac5-v2.merge-list, non noq replace
			!${plink} --bfile ALL.chr23.phase3_shapeit2_mvncall_integrated_v1b.20130502-mac5-v2 --update-sex all-1000g-phase3-chrall-mac5-v2.update-sex --merge-list ALL.chrall.phase3_shapeit2_mvncall_integrated_v5a.20130502-mac5-v2.merge-list --make-bed --out all-1000g-phase3-chrall-mac5-v2
			noi checkfile, file(${init_root}/`function'/ftp.1000genomes.ebi.ac.uk/all-1000g-phase3-chrall-mac5-v2.bed)
			noi checkfile, file(${init_root}/`function'/ftp.1000genomes.ebi.ac.uk/all-1000g-phase3-chrall-mac5-v2.bim)
			noi checkfile, file(${init_root}/`function'/ftp.1000genomes.ebi.ac.uk/all-1000g-phase3-chrall-mac5-v2.fam)
			noi di as text"#########################################################################"
			noi di as text""	
			}
		}
	noi di as text"#########################################################################"
	noi di as text""	
	}

qui { // module 7 - create all-1000g-phase3-chrall-mac5-v2.aims
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 7: create all-1000g-phase3-chrall-mac5-v2.aims file"
	noi di as text"#########################################################################"
	local already_created ${init_root}/`function'/ftp.1000genomes.ebi.ac.uk/all-1000g-phase3-chrall-mac5-v2.aims
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		}
	else {
		import delim using all-1000g-phase3-chrall-mac5-v2.population, clear varnames(nonames)
		rename (v1-v3) (fid iid pop)
		outsheet fid iid pop using all-1000g-phase3-chrall-mac5-v2.within, replace noq
		!${plink} --bfile all-1000g-phase3-chrall-mac5-v2 ///
		--fst ///
		--within all-1000g-phase3-chrall-mac5-v2.within ///
		--out all-1000g-phase3-chrall-mac5-v2
		import delim using all-1000g-phase3-chrall-mac5-v2.fst, clear
		keep snp fst
		drop if fst == "nan"
		destring fst, replace
		keep if fst > .05
		outsheet snp using all-1000g-phase3-chrall-mac5-v2.aims, non noq replace
		noi checkfile, file(`already_created')
		}
	noi di as text"#########################################################################"
	noi di as text""
	}
qui { // module 8 - create gbr-1000g-phase3-chrall-mac5-v2.bed
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 8: create gbr-1000g-phase3-chrall-mac5-v2.bed file"
	noi di as text"#########################################################################"
	local already_created ${init_root}/`function'/ftp.1000genomes.ebi.ac.uk/gbr-1000g-phase3-chrall-mac5-v2.bed
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		}
	else {
		import delim using all-1000g-phase3-chrall-mac5-v2.population, clear varnames(nonames)
		rename (v1-v3) (fid iid pop)
		keep if pop == "GBR"
		outsheet fid iid using gbr-1000g-phase3-chrall-mac5-v2.keep, replace noq non
		!${plink} --bfile all-1000g-phase3-chrall-mac5-v2 ///
		--keep gbr-1000g-phase3-chrall-mac5-v2.keep ///
		--make-bed ///
		--out gbr-1000g-phase3-chrall-mac5-v2
		noi checkfile, file(`already_created')
		}
	noi di as text"#########################################################################"
	noi di as text""
	}
qui { // module 9 - clean folder
	files2dta, dir(`c(pwd)')
	gen keep = .
	replace keep = 1 if files == "gbr-1000g-phase3-chrall-mac5-v2.bed"
	replace keep = 1 if files == "gbr-1000g-phase3-chrall-mac5-v2.bim"
	replace keep = 1 if files == "gbr-1000g-phase3-chrall-mac5-v2.fam"
	replace keep = 1 if files == "gbr-1000g-phase3-chrall-mac5-v2.log"
	replace keep = 1 if files == "all-1000g-phase3-chrall-mac5-v2.bed"
	replace keep = 1 if files == "all-1000g-phase3-chrall-mac5-v2.bim"
	replace keep = 1 if files == "all-1000g-phase3-chrall-mac5-v2.fam"
	replace keep = 1 if files == "all-1000g-phase3-chrall-mac5-v2.log"
	replace keep = 1 if files == "all-1000g-phase3-chrall-mac5-v2.aims"
	replace keep = 1 if files == "all-1000g-phase3-chrall-mac5-v2.population"
	replace keep = 1 if files == "all-1000g-phase3-chrall-mac5-v2.super-population"
	replace keep = 1 if files == "all-1000g-phase3-chrall-mac5-v2-panel.dta"
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
	