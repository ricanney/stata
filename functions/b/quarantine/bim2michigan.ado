*! 1.0.1 Richard Anney 19dec2018
* - all samples should be processed (genotypeqc) against common reference
* - need to adapt for unix 

program bim2michigan
global  bim2michigan_version v1
syntax , bim(string asis) ref(string asis)

noi di as text" "
noi di as text"#########################################################################"
noi di as text"# bim2michigan"
noi di as text"#########################################################################"
noi di as text"# Started: $S_DATE $S_TIME"
noi di as text"# Username:            `c(username)'"
noi di as text"# Version:             "as result"${bim2michigan_version}"
noi di as text"# Operating System:    `c(os)'"
noi di as text"#########################################################################"


qui { // module 1 - prepare folders 
	qui { // define short names
			clear
			set obs 2
			gen a = ""
			replace a =  "`bim'" in 1
			replace a =  "`ref'" in 2
			replace a = subinstr(a,"\","%%",.)
			replace a = subinstr(a,"/","%%",.)
			split a, p("%%")
			gen a999 = ""
			for var a1-a999: replace a1 = X if X !=""
			keep a1
			compress
			replace a1 = "global bim2michigan_s_bim " + a1 in 1
			replace a1 = "global bim2michigan_s_ref " + a1 in 2
			replace a1 = subinstr(a1,"_bim.dta","",.)
			outsheet a1 using tmpa.do, non noq replace
			do tmpa.do
			erase tmpa.do 
			}
	qui { // define working directory 
		noi di as text""	
		noi di as text"#########################################################################"
		noi di as text"# SECTION - 1: check folders and files"
		noi di as text"#########################################################################"
		noi di as text"# > bim2michigan ........................ project name is "as result"$bim2michigan_s_bim"
		noi di as text"# > bim2michigan ...................... current folder is "as result"`c(pwd)'"
		}
	qui { // create project directory
	  !mkdir $bim2michigan_s_bim
	  cd $bim2michigan_s_bim
		noi di as text"# > bim2michigan ....................... create folder is "as result"`c(pwd)'"
		}
	qui { // create temporary directory
		qui create_temp_dir	
		noi di as text"# > bim2michigan .................... temporary folder is "as result"`c(pwd)'"
		noi di as text"#########################################################################"
		}	
	}
qui { // module 2 - define files
	noi di as text""	
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 2: define files"
	noi di as text"#########################################################################"
	qui { // check files
		foreach file in bim bed fam { 
			noi checkfile, file(`bim'.`file')
			}
		noi checkfile, file(`ref')
		checkfile, file(${plink})
		noi di as text"#########################################################################"
		}
	}
qui { // module 3 - create vcf files
		noi di as text""	
		noi di as text"#########################################################################"
		noi di as text"# SECTION - 3: create vcf files"
		noi di as text"#########################################################################"
		qui { // identify _bim.dta file (loc_name)
			noi di as text"# > bim2michigan ........................................ identify _bim.dta files"
			noi checkfile,  file(`bim'_bim.dta)
			capture confirm file `bim'_bim.dta
			if !_rc {
				}
			else {
				bim2dta, bim(`bim')
				}
			}
		qui { // use data
			noi di as text"# > bim2michigan ........................................ use _bim.dta files"
			use `bim'_bim.dta, clear
			keep loc_name snp a1 a2 gt
			for var snp a1 a2 gt: rename X bim_X
			}
		qui { // merge with hrc 1.1 reference
			noi di as text"# > bim2michigan ........................... merging with "as result"`ref'"
			merge 1:1 loc_name using `ref'
			keep if _m == 3
			for var snp a1 a2 gt: rename X michigan_X
			}
		qui { // select snps
			noi di as text"# > bim2michigan ........................................ extracting overlapping snps"
			count
			noi di as text"# > bim2michigan ....... SNPs after merge / pre-clean (N) "as result `r(N)'
			outsheet bim_snp using bim2michigan.extract, non noq replace
			!$plink --bfile `bim' --extract bim2michigan.extract --make-bed --out bim2michigan-1
		  erase bim2michigan.extract
		  }
		qui { // flip strand
			noi di as text"# > bim2michigan ........................................ flipping strands"
			gen flip = .
			replace flip = 1 if bim_gt != michigan_gt
			outsheet bim_snp using bim2michigan.flip if flip == 1, non noq replace
			!$plink --bfile bim2michigan-1 --flip bim2michigan.flip --make-bed --out bim2michigan-2
			erase bim2michigan.flip
		}
		qui { // rename bim_snps to michigan_snps
			noi di as text"# > bim2michigan ........................................ renaming snps"
			outsheet bim_snp michigan_snp if bim_snp != michigan_snp using bim2michigan.update-name, non noq replace
			!$plink --bfile bim2michigan-2 --update-name bim2michigan.update-name --make-bed --out bim2michigan-3
			erase bim2michigan.update-name
			}
		qui { // set reference allele
			noi di as text"# > bim2michigan ........................................ setting reference allele"
			outsheet michigan_snp michigan_a2 using bim2michigan.a2-allele, non noq replace
			!$plink --bfile bim2michigan-3 --a2-allele bim2michigan.a2-allele --make-bed --out bim2michigan-4
			erase bim2michigan.a2-allele
			}
		qui { // convert to michigan vcf
			noi di as text"# > bim2michigan ........................................ converting to vcf"
			foreach num of num 1/23 { 
				!$plink --bfile bim2michigan-4 --keep-allele-order -chr `num' --recode vcf --out ${bim2michigan_s_bim}-chrom`num'_hg19_1
				}
			!bash -c "cat ${bim2michigan_s_bim}-chrom23_hg19_1.vcf | sed -e 's/^23/X/' > ${bim2michigan_s_bim}-chromX_hg19_1.vcf"
			erase         ${bim2michigan_s_bim}-chrom23_hg19_1.vcf
			}
		qui { // prepare for bgzip on rock
			noi di as text"# > bim2michigan ........................................ archiving vcf bundle"
			!bash -c "tar -zcvf ${bim2michigan_s_bim}.tar.gz *.vcf"
			noi di as text"# > bim2michigan ........................................ creating bash file to convert using bgzip"
			clear
			set obs 23
			gen N = _n
			tostring N, replace
			replace N = "X" in 23
			gen a = ""
			replace a = "/share/apps/bgzip ${bim2michigan_s_bim}-chrom" + N + "_hg19_1.vcf & "
			outsheet a using bim2michigan.sh, non noq replace
			!$dos2unix bim2michigan.sh
			!del *.vcf *.log
			noi di as text"#########################################################################"

			}
		}
qui { // module 3 - create bim2merge.log
	noi di as text""	
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3: create log file"
	noi di as text"#########################################################################"
	noi di as text"# > bim2michigan .................................... log "as result"bim2michigan.log"
	qui { // write meta log
		file open myfile using "bim2michigan.log", write replace
		file write myfile "#########################################################################" _n
		file write myfile "# bim2michigan"                                                         _n
		file write myfile "#########################################################################" _n
		file write myfile "# Started:             $S_DATE $S_TIME "                                   _n
		file write myfile "# Username:            `c(username)'"                                      _n
		file write myfile`"# Version:             ${bim2michigan_version}"'                        _n
		file write myfile`"# Operating System:    `c(os)'"'                                           _n
		file write myfile "#########################################################################" _n
		file write myfile "" _n
		file write myfile "#########################################################################" _n
  	file write myfile "# Dataset/s " _n
		file write myfile "# > bim2michigan ............................. dataset is ${bim2michigan_s_bim}" _n
		file write myfile "# > bim2michigan ........................... reference is ${bim2michigan_s_ref}" _n
   	file write myfile "#########################################################################" _n
  	file write myfile "# ${bim2michigan_s_bim}" _n
  	file write myfile "# pre-process " _n
		bim2count, bim(`bim')
		file write myfile "# > bim2michigan ............................. total SNPs ${bim2count_snp}" _n
		file write myfile "# > bim2michigan ............................ individuals ${bim2count_ind}" _n
   	file write myfile "#########################################################################" _n
		file write myfile "# post-process" _n
		bim2count, bim(bim2michigan-4)
		file write myfile "# > bim2michigan ............................. total SNPs ${bim2count_snp}" _n
		file write myfile "# > bim2michigan ............................ individuals ${bim2count_ind}" _n
		file write myfile "#########################################################################" _n
		file write myfile "# the ${bim2michigan_s_bim}.tar.gz bundle and bim2michigan.sh " _n
		file write myfile "# should be transferred to rocks, unpacked and repacked using the bgzip" _n
		file write myfile "# function from tabix. " _n
		file write myfile "#  " _n
		file write myfile "# the *.vcf.gz files are then ready to submit to the michigan imputation " _n
		file write myfile "# server at https://imputationserver.sph.umich.edu" _n
		file write myfile "# > michigan imputation server (options) ................ unphased" _n
		file write myfile "# > michigan imputation server (options) ................ hrc.r1.1.2016" _n
		file write myfile "# > michigan imputation server (options) ................ shapeit" _n
		file write myfile "#########################################################################" _n
		file close myfile
		}
	noi di as text"#########################################################################"
	}	
qui { // module 4 - move and clean
	noi di as text""	
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 4: move and clean up temporary files / folders"
	clear
	set obs 1
	gen os = "`c(os)'"
	if os == "Unix" { 
		!cp "bim2michigan.log" "../${bim2michigan_s_bim}-bim2michigan.log"
		!cp "${bim2michigan_s_bim}.tar.gz" "../${bim2michigan_s_bim}.tar.gz"
		!cp "bim2michigan.sh" "../bim2michigan.sh"
		}
	else if os == "Windows" {
		!copy /v /y "bim2michigan.log" "../${bim2michigan_s_bim}-bim2michigan.log"
		!copy /v /y "${bim2michigan_s_bim}.tar.gz" "../${bim2michigan_s_bim}.tar.gz"
		!copy /v /y "bim2michigan.sh" "../bim2michigan.sh"
		}
	delete_temp_dir
	}
noi di as text"#########################################################################"
noi di as text"# Completed: $S_DATE $S_TIME"
noi di as text"#########################################################################"
end;	
