*! 1.0.7 Richard Anney 05nov2018
*! 1.0.8 Richard Anney 06dec2018
* - update allele frequency pruning
* - update bim2ancestry method 
* - remove *.like file creation
* - check on unix server (rocks)
*! 1.0.9 Richard Anney 23jan2019
* - update code to account for the additional subroutine folders
* - update code to include new bim2graph* routines
*! 1.0.9 Richard Anney 5feb2019
* - add routine for admixed sample
*! 1.0.9 Richard Anney 30apr2019
* - update upfront functions

program  genotypeqc
local    function genotypeqc
local    `function'_version v9
syntax, bim(string asis) [known_array(string asis) admixed(string asis) bim2frq_compare_sd(real 4)]

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
*
qui { // <1> collect files
	qui { // module 1 - define thresholds
		noi di as text"#########################################################################"
		noi di as text"# SECTION - 1.1: define thresholds"
		noi di as text"#########################################################################"
		noi di as text"# > ............ --mac "as result"5"
		noi di as text"# > ........... --geno "as result"${geno1}"as text";"as result"${geno2}"
		noi di as text"# > ........... --mind "as result"${mind}"
		noi di as text"# > .......... --hardy "as result"1e-${hwep}"
		noi di as text"# > ... heterozygosity "as result"${hetsd}"as text" standard deviations from mean"
		noi di as text"# > .......... kinship "
		noi di as text"# > ....... duplicates "as result"${kin_d}"
		noi di as text"# > ....... 1st-degree "as result"${kin_f}"
		noi di as text"# > ....... 2nd-degree "as result"${kin_s}"
		noi di as text"# > ....... 3rd-degree "as result"${kin_t}"
		noi di as text"#########################################################################"
		noi di as text""	
		}
	qui { // module 2 - prepare folders 
		noi di as text"#########################################################################"
		noi di as text"# SECTION - 1.2: prepare folders"
		noi di as text"#########################################################################"
		bim_path_to_short, bim(`bim')
		local `function'_bim_short ${bim_short}
		ref_path_to_short, ref($bim2refid_ref)
		local `function'_bim2refid_short ${ref_short}
		ref_path_to_short, ref($bim2frq_compare_ref)
		local `function'_bim2frq_short ${ref_short}	
		ref_path_to_short, ref($bim2ancestry_ref)
		local `function'_bim2ancestry_short ${ref_short}			
		cd ${init_root}	
		!mkdir `function'
		cd `function'
		noi di as text"# > ..... project name "as result"``function'_bim_short'"
		!mkdir ``function'_bim_short'
		cd ``function'_bim_short'
		noi di as text"#########################################################################"
		noi di as text""	
		}
	qui { // module 3 - define files
		noi di as text"#########################################################################"
		noi di as text"# SECTION - 1.3: define files"
		noi di as text"#########################################################################"
		noi di as text"# > ......... locating files to process"
		foreach file in bim bed fam {
			noi checkfile, file(`bim'.`file')
			}
		noi di as text"# > ......... locating PLINK files"
		noi checkfile, file(${plink})
		noi checkfile, file(${plink2})
		checktabbed
		noi di as text"# > ......... locating bim2build files"
		noi checkfile, file(${bim2build_ref})
		noi di as text"# > ......... locating bim2refid reference files"
		bim_path_to_short, bim(${bim2refid_ref})
		local `function'_bim2refid_short ${bim_short}
		noi checkfile, file(${init_root}/bim2dta/``function'_bim2refid_short'/``function'_bim2refid_short'_bim.dta)
		noi di as text"# > ......... locating bim2frq_compare reference files"
		foreach file in bim bed fam {
			noi checkfile, file(${bim2frq_compare_ref}.`file')
			}
		noi di as text"# > ......... locating bim2ancestry reference files"
		foreach file in bim bed fam aims population { 
			noi checkfile, file(${bim2ancestry_ref}.`file')
			}		
		noi di as text"#########################################################################"
		noi di as text""	
		}
	qui { // module 4 - copy files to working folder
		noi di as text"#########################################################################"
		noi di as text"# SECTION - 1.4: copy genotypes to ${data_in}"
		noi di as text"#########################################################################"
		which_os
		if os == "Unix" { 
			foreach file in bim bed fam array { 
				capture confirm file ``function'_bim_short'.`file'
				if !_rc {
					noi checkfile, file(``function'_bim_short'.`file')
					}
				else {
					!cp "`bim'.`file'"   "``function'_bim_short'.`file'"
					noi checkfile, file(``function'_bim_short'.`file')
					}
				}
			}
		else if os == "Windows" {
			foreach file in bim bed fam array { 
				capture confirm file ``function'_bim_short'.`file'
				if !_rc {
					noi checkfile, file(``function'_bim_short'.`file')
					}
				else {
					!copy /v /y "`bim'.`file'.`file'" "``function'_bim_short'.`file'"
					noi checkfile, file(``function'_bim_short'.`file')
					}
				}
			}
		noi di as text"#########################################################################"
		noi di as text""
		}
	}
*
qui { // <2> apply bim2array
	qui { // module 1 - check for known_array
		noi di as text"#########################################################################"
		noi di as text"# SECTION - 2.1: check if array is known"
		noi di as text"#########################################################################"
		clear
		set obs 1
		gen known_array = "`known_array'"
		if known_array == "" {
			noi di as text"# > .................. "as error"array not specified"
			noi di as text"#########################################################################"
			noi di as text""
			noi di as text"#########################################################################"
			noi di as text"# SECTION - 2.2: array not known - check for *.array file"
			noi di as text"#########################################################################"
			local already_created ``function'_bim_short'.array
			capture confirm file `already_created'
			if !_rc {
					noi checkfile, file(`already_created')
					import delim  using ``function'_bim_short'.array, clear varnames(1) case(preserve)
					gsort -j
					local bim2array_array `=array[1]'
					local bim2array_jaccard `=jaccardindex[1]'
					noi di as text"# > . array defined as "as result"`bim2array_array'"
					noi di as text"# > .... jaccard index "as result"`bim2array_jaccard'"				
					}
			else {
					local already_created ${init_root}/bim2array/``function'_bim_short'/``function'_bim_short'.array
					capture confirm file `already_created'
					if !_rc {
						noi checkfile, file(`already_created')
						}
					else {
						noi checkfile, file(`already_created')
						noi di as text"# > .................. create "as result "``function'_bim_short'.array"
						noi bim2array, bim(`bim') dir($bim2array_dir)
						import delim  using `already_created', clear varnames(1) case(preserve)
						gsort -j
						local bim2array_array `=array[1]'
						local bim2array_jaccard `=jaccardindex[1]'
						noi di as text"# > . array defined as "as result"`bim2array_array'"
						noi di as text"# > .... jaccard index "as result"`bim2array_jaccard'"	
						}						
					}
			noi di as text"#########################################################################"
			noi di as text""
			}
		else { 
			noi di as text"#########################################################################"
			noi di as text"# SECTION - 2.2: define genotyping array from known_array"
			noi di as text"#########################################################################"
			local bim2array_array `known_array'
			local bim2array_jaccard NaN
			noi di as text"# > . array defined as "as result"`bim2array_array'"
			noi di as text"# > .... jaccard index "as result"`bim2array_jaccard'"
			file open myfile using "``function'_bim_short'.array", write replace
			file write myfile `"array	jaccardindex"' _n
			file write myfile `"`bim2array_array'	`bim2array_jaccard'"' _n
			file close myfile	
			noi di as text"#########################################################################"
			noi di as text""
			}
	

		}
	qui { // module 3 - update affyid arrays
		noi di as text"#########################################################################"
		noi di as text"# SECTION - 2.3: updating rsid if array is affyid"
		noi di as text"#########################################################################"	
		clear
		set obs 1
		gen array = "`bim2array_array'"
		gen update = .
		replace update = 1 if array == "genomewidesnp-6-na35-affyid"  
		replace update = 1 if array == "genomewidesnp-5-na35-affyid"
		if update == 1 { 
			noi di as text"# > .... array type is "as error"affyid"
			use ${bim2array_dir}/`bim2array_array'.dta, replace
			keep snp rsid 
			order snp rsid
			drop if snp == ""
			drop if rsid == ""
			duplicates tag snp, gen(tag)
			keep if tag == 0
			drop tag
			duplicates tag rsid, gen(tag)
			keep if tag == 0
			drop tag
			create_temp_dir
			outsheet rsid using affy-to-rsid.extract, non noq replace
			outsheet snp rsid using affy-to-rsid.update-name, non noq replace		
			!$plink --bfile `bim'        --extract     affy-to-rsid.extract     --make-bed --out tmp_affyid
			!$plink --bfile tmp_affyid   --update-name affy-to-rsid.update-name --make-bed --out ../``function'_bim_short'
			delete_temp_dir
			}
		else if update == . {
			noi di as text"# > .... array type is "as result"not affyid"
			}
		}		
	noi di as text"#########################################################################"
	noi di as text""
	}
*
qui { // <3> prepare files for QC
	qui { // module 1 - count individuals and markers
		noi di as text"#########################################################################"
		noi di as text"# SECTION - 3.1: count SNPs and individuals prior to QC "
		noi di as text"#########################################################################"
		bim2count, bim(${init_root}/genotypeqc/``function'_bim_short'/``function'_bim_short')
		noi di as text"# > ......... binaries "as result"``function'_bim_short'"
		noi di as text"# > ..... SNPs in file "as result "${bim2count_snp}"
		noi di as text"# > .... indiv in file "as result "${bim2count_ind}"
		noi di as text"#########################################################################"
		noi di as text""
		}
	qui { // module 2 - remove low count markers
		noi di as text"#########################################################################"
		noi di as text"# SECTION - 3.2: limit to --chr 1-23, --mac 5 and --geno 0.99 "
		noi di as text"#########################################################################"
		local already_created ${init_root}/genotypeqc/``function'_bim_short'/``function'_bim_short'-mac5.bed
		capture confirm file `already_created'
		if !_rc {
			noi checkfile, file(`already_created')
			}
		else {
			!$plink 	--bfile ``function'_bim_short' ///
					--chr 1-23 --mac 5 --geno 0.99  --make-bed ///
					--out ``function'_bim_short'-mac5
			noi checkfile, file(`already_created')
			}
		bim2count, bim(${init_root}/genotypeqc/``function'_bim_short'/``function'_bim_short'-mac5)
		noi di as text"# > ......... binaries "as result"``function'_bim_short'-mac5"
		noi di as text"# > ..... SNPs in file "as result "${bim2count_snp}"
		noi di as text"# > .... indiv in file "as result "${bim2count_ind}"
		noi di as text"#########################################################################"
		noi di as text""
		}
	}
*
qui { // <4> apply bim2build
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 4: define build from *.build file"
	noi di as text"#########################################################################"
	local already_created ${init_root}/genotypeqc/``function'_bim_short'/``function'_bim_short'-mac5-hg19.bed
	capture confirm file `already_created'
	if !_rc {
		noi di as text"# > ......... build is "as result"hg19 +1"
		bim2count, bim(${init_root}/genotypeqc/``function'_bim_short'/``function'_bim_short'-mac5-hg19)
		noi di as text"# > ......... binaries "as result"``function'_bim_short'-mac5-hg19"
		noi di as text"# > ..... SNPs in file "as result "${bim2count_snp}"
		noi di as text"# > .... indiv in file "as result "${bim2count_ind}"				
		}
	else {
		local already_created ${init_root}/bim2build/``function'_bim_short'-mac5/``function'_bim_short'-mac5.build
		capture confirm file `already_created'
		if !_rc {
			noi checkfile, file(`already_created')
			}
		else {
			bim2build, bim(${init_root}/genotypeqc/``function'_bim_short'/``function'_bim_short'-mac5) ref($bim2build_ref)
			noi checkfile, file(`already_created')
			}
		import delim using  `already_created', clear varnames(1) case(preserve)
		if build == "hg19 +1" {
			noi di as text"# > ......... build is "as result"hg19 +1"
			!$plink --bfile  ``function'_bim_short'-mac5 --make-bed --out ``function'_bim_short'-mac5-hg19	
			}
		else {
			global bim2build `=build[1]'
			noi di as text"# > ......... build is "as error"${bim2build}"
			noi di as text"# > ....... convert to "as result"hg19 +1"
			qui { // define markers in `out'
				local already_created ${init_root}/bim2dta/``function'_bim_short'-mac5/``function'_bim_short'-mac5_bim.dta
				capture confirm file `already_created'
				if !_rc {
					noi checkfile, file(`already_created')
					}
				else {
					noi di as text"# > .................. create "as result "``function'_bim_short'-mac5_bim.dta"
					bim2dta, bim(`out')
					noi checkfile, file(`already_created')
					}			
				}
			qui { // update chr map
				use `already_created', clear
				keep snp
				merge 1:1 snp using ${bim2array_dir}/`bim2array_array'.dta
				keep if _m == 3
				create_temp_dir
				outsheet snp using      ``function'_bim_short'-mac5-bim2build.extract, non noq replace
				outsheet snp chr using  ``function'_bim_short'-mac5-bim2build.update-chr, non noq replace
				outsheet snp bp  using  ``function'_bim_short'-mac5-bim2build.update-map, non noq replace
				!$plink --bfile  ${init_root}/genotypeqc/``function'_bim_short'/``function'_bim_short'-mac5          --extract    ``function'_bim_short'-mac5-bim2build.extract    --make-bed --out bim2build-v1.1
				!$plink --bfile  bim2build-v1.1 --update-chr ``function'_bim_short'-mac5-bim2build.update-chr --make-bed --out bim2build-v1.2 
				!$plink --bfile  bim2build-v1.2 --update-map ``function'_bim_short'-mac5-bim2build.update-map --make-bed --out bim2build-v1.3
				!$plink --bfile  bim2build-v1.3 --chr 1-23                                                    --make-bed --out ${init_root}/genotypeqc/``function'_bim_short'/``function'_bim_short'-mac5-hg19
				delete_temp_dir
				}
			qui { // re-check build
				bim2build, bim(${init_root}/genotypeqc/``function'_bim_short'/``function'_bim_short'-mac5-hg19) ref($bim2build_ref)
				noi di as text"# > ......... build is "as result"${bim2build}"
				bim2count, bim(${init_root}/genotypeqc/``function'_bim_short'/``function'_bim_short'-mac5-hg19)
				noi di as text"# > ......... binaries "as result"``function'_bim_short'-mac5-hg19"
				noi di as text"# > ..... SNPs in file "as result "${bim2count_snp}"
				noi di as text"# > .... indiv in file "as result "${bim2count_ind}"
				}
			}
		}
	noi di as text"#########################################################################"
	noi di as text""
	}
*
qui { // <5> apply bim2refid
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 5: update name (to reference) "
	noi di as text"#########################################################################"
	local already_created ${init_root}/bim2refid/``function'_bim_short'-mac5-hg19/``function'_bim_short'-mac5-hg19-refid.bed
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		}
	else {
		bim2refid, bim(${init_root}/genotypeqc/``function'_bim_short'/``function'_bim_short'-mac5-hg19) ref($bim2refid_ref)
		noi checkfile, file(`already_created')
		}
	noi di as text"# > ..... SNP names as "as result"``function'_bim2refid_short'"
	bim2count, bim(${init_root}/bim2refid/``function'_bim_short'-mac5-hg19/``function'_bim_short'-mac5-hg19-refid)
	noi di as text"# > ......... binaries "as result"``function'_bim_short'-mac5-hg19-refid"
	noi di as text"# > ..... SNPs in file "as result "${bim2count_snp}"
	noi di as text"# > .... indiv in file "as result "${bim2count_ind}"	
	noi di as text"#########################################################################"
	noi di as text""
	}
*
*error in module 6 related to admixed addition*
*
qui { // <6> apply bim2frq_compare
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 6: exclude divergent allele frequencies "
	noi di as text"#########################################################################"
	local already_created ${init_root}/genotypeqc/``function'_bim_short'/``function'_bim_short'-mac5-hg19-refid-bim2frq.bed
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		}
	else {
		clear
		set obs 1
		gen admixed = "`admixed'"
		if admixed == "" {
			noi di as text"# > .................. "as result "sample not admixed - select on allele frequency"
			qui { // define SNPs to exclude
				local already_created2 ${init_root}/bim2frq_compare/``function'_bim_short'-mac5-hg19-refid/``function'_bim_short'-mac5-hg19-refid-``function'_bim2frq_short'_bim2frq_compare.exclude
				capture confirm file `already_created2'
				if !_rc {
					noi checkfile, file(`already_created2')
					}
				else {
					bim2frq_compare, bim(${init_root}/bim2refid/``function'_bim_short'-mac5-hg19/``function'_bim_short'-mac5-hg19-refid) ref(${bim2frq_compare_ref}) sd(`bim2frq_compare_sd')
					noi checkfile, file(`already_created2')
					}
				}
			qui { // parse genotypes
				noi di as text"# > .................. "as result "update PLINK binaries"
				!$plink 	--bfile    ${init_root}/bim2refid/``function'_bim_short'-mac5-hg19/``function'_bim_short'-mac5-hg19-refid ///
						--exclude  ${init_root}/bim2frq_compare/``function'_bim_short'-mac5-hg19-refid/``function'_bim_short'-mac5-hg19-refid-``function'_bim2frq_short'_bim2frq_compare.exclude ///
						--make-bed ///
						--out ``function'_bim_short'-mac5-hg19-refid-bim2frq
				}
			}
		else {
			noi di as text"# > .................. "as result "sample is admixed - no selection on allele frequency"
			noi di as text"# > .................. "as result "update PLINK binaries"
			!$plink 	--bfile    ${init_root}/bim2refid/``function'_bim_short'-mac5-hg19/``function'_bim_short'-mac5-hg19-refid ///
					--make-bed ///
					--out ``function'_bim_short'-mac5-hg19-refid-bim2frq
			!mkdir ${init_root}/bim2frq_compare/``function'_bim_short'-mac5-hg19-refid
			tw scatteri 1 1, msymbol(i) ylab("") xlab("") ytitle("") xtitle("") yscale(off) xscale(off) plotregion(lpattern(blank)) nodraw saving(${init_root}/bim2frq_compare/``function'_bim_short'-mac5-hg19-refid/``function'_bim_short'-mac5-hg19-refid-``function'_bim2frq_short'_bim2frq_compare.gph, replace)
			}
		local already_created ${init_root}/genotypeqc/``function'_bim_short'/``function'_bim_short'-mac5-hg19-refid-bim2frq.bed
		noi checkfile, file(`already_created')
		}	
	bim2count, bim(${init_root}/genotypeqc/``function'_bim_short'/``function'_bim_short'-mac5-hg19-refid-bim2frq)
	noi di as text"# > ......... binaries "as result"``function'_bim_short'-mac5-hg19-refid-bim2frq"
	noi di as text"# > ........ reference "as result"``function'_bim2frq_short'"
	noi di as text"# > ..... SNPs in file "as result "${bim2count_snp}"
	noi di as text"# > .... indiv in file "as result "${bim2count_ind}"	
	noi di as text"#########################################################################"
	noi di as text""
	}
*		
qui { // <7> perform pre-QC metrics
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 7: calculate / plot pre-quality control metrics "
	noi di as text"#########################################################################"
	cd ${init_root}/genotypeqc/``function'_bim_short'
	local in  ${init_root}/genotypeqc/``function'_bim_short'/``function'_bim_short'-mac5-hg19-refid-bim2frq
	qui { // frq
		local graph bim2graphfrq
		noi di as text"# > ............ apply "as result "`graph'"
		local already_created ${init_root}/`graph'/``function'_bim_short'-mac5-hg19-refid-bim2frq/``function'_bim_short'-mac5-hg19-refid-bim2frq_`graph'.gph
		capture confirm file `already_created'
		if !_rc {
			noi checkfile, file(`already_created')
			}
		else {	
			`graph', bim(`in')
			noi checkfile, file(`already_created')
			}
		}
	qui { // het
		local graph bim2graphhet
		noi di as text"# > ............ apply "as result "`graph'"
		local already_created ${init_root}/`graph'/``function'_bim_short'-mac5-hg19-refid-bim2frq/``function'_bim_short'-mac5-hg19-refid-bim2frq_`graph'.gph
		capture confirm file `already_created'
		if !_rc {
			noi checkfile, file(`already_created')
			}
		else {	
			`graph', bim(`in')
			noi checkfile, file(`already_created')
			}
		}
	qui { // hwe
		local graph bim2graphhwe
		noi di as text"# > ............ apply "as result "`graph'"
		local already_created ${init_root}/`graph'/``function'_bim_short'-mac5-hg19-refid-bim2frq/``function'_bim_short'-mac5-hg19-refid-bim2frq_`graph'.gph
		capture confirm file `already_created'
		if !_rc {
			noi checkfile, file(`already_created')
			}
		else {	
			`graph', bim(`in')
			noi checkfile, file(`already_created')
			}
		}
	qui { // kin0	
		local graph  bim2graphkin0
		qui { // hist
			local graph2 hist
			noi di as text"# > ............ apply "as result "`graph'"
			local already_created ${init_root}/`graph'/``function'_bim_short'-mac5-hg19-refid-bim2frq/``function'_bim_short'-mac5-hg19-refid-bim2frq_`graph'-`graph2'.gph
			capture confirm file `already_created'
			if !_rc {
				noi checkfile, file(`already_created')
				}
			else {	
				`graph', bim(`in')
				noi checkfile, file(`already_created')
				}
			}
		qui { // scatter
			local graph2 scatter
			noi di as text"# > ............ apply "as result "`graph'"
			local already_created ${init_root}/`graph'/``function'_bim_short'-mac5-hg19-refid-bim2frq/``function'_bim_short'-mac5-hg19-refid-bim2frq_`graph'-`graph2'.gph
			capture confirm file `already_created'
			if !_rc {
				noi checkfile, file(`already_created')
				}
			else {	
				`graph', bim(`in')
				noi checkfile, file(`already_created')
				}
			}
		}
	qui { // miss	
		local graph  bim2graphmiss
		qui { // lmiss
			local graph2 lmiss
			noi di as text"# > ............ apply "as result "`graph'"
			local already_created ${init_root}/`graph'/``function'_bim_short'-mac5-hg19-refid-bim2frq/``function'_bim_short'-mac5-hg19-refid-bim2frq_`graph'-`graph2'.gph
			capture confirm file `already_created'
			if !_rc {
				noi checkfile, file(`already_created')
				}
			else {	
				`graph', bim(`in')
				noi checkfile, file(`already_created')
				}
			}
		qui { // imiss
			local graph2 imiss
			noi di as text"# > ............ apply "as result "`graph'"
			local already_created ${init_root}/`graph'/``function'_bim_short'-mac5-hg19-refid-bim2frq/``function'_bim_short'-mac5-hg19-refid-bim2frq_`graph'-`graph2'.gph
			capture confirm file `already_created'
			if !_rc {
				noi checkfile, file(`already_created')
				}
			else {	
				`graph', bim(`in')
				noi checkfile, file(`already_created')
				}
			}
		}
	noi di as text"#########################################################################"
	noi di as text""
	}
*
qui { // <8> apply quality control 
	local in  ${init_root}/genotypeqc/``function'_bim_short'/``function'_bim_short'-mac5-hg19-refid-bim2frq
	local qc  ${init_root}/genotypeqc/``function'_bim_short'/``function'_bim_short'-qc-``function'_version'
	local already_created `qc'.bed
	capture confirm file `already_created'
	if !_rc {
		noi di as text"#########################################################################"
		noi di as text"# SECTION - 8.1: apply quality control (completed)"
		noi di as text"#########################################################################"
		noi checkfile, file(`already_created')
		bim2count, bim(`qc')
		noi di as text"# > ..... binaries out "as result"`qc'"
		noi di as text"# > ..... SNPs in file "as result "${bim2count_snp}"
		noi di as text"# > .... indiv in file "as result "${bim2count_ind}"	
		noi di as text"#########################################################################"
		noi di as text""
		}
	else {
		qui { // module 1 - apply round 1
			noi di as text"#########################################################################"
			noi di as text"# SECTION - 8.1: apply quality control metrics (ROUND 1)"
			noi di as text"#########################################################################"
			local already_created `in'-round1.bed
			capture confirm file `already_created'
			if !_rc {
				noi checkfile, file(`already_created')
				local in `in'-round1
				}
			else {	
				create_temp_dir
				qui { // remove out of bound (heterozygosity)
					local out             ``function'_bim_short'-temp-r1.1
					bim_path_to_short, bim(`in')
					local in_file ${bim_short}
					noi di as text"# > ............ apply "as result "bim2graphhet.remove"
					noi di as text"# > ...... binaries in "as result"`in'"
					local het ${init_root}/bim2graphhet/`in_file'/`in_file'_bim2graphhet.remove
					import delim using `het',clear varn(noname)
					count
					if `r(N)' > 0 { 
						noi di as text"# > ........... remove "as result "`r(N)'" as text" individuals"
						!$plink --bfile `in' --remove `het' --set-hh-missing --make-bed --out `out'
						bim2count, bim(`c(pwd)'/`out')
						noi di as text"# > ..... binaries out "as result"`out'"
						noi di as text"# > ..... SNPs in file "as result "${bim2count_snp}"
						noi di as text"# > .... indiv in file "as result "${bim2count_ind}"	
						local in `out'
						}
					else { 
						noi di as text"# > ........... remove "as result "0" as text" individuals"
						}
					}
				qui { // remove out of bound (hardy-weinberg)
					local out ``function'_bim_short'-temp-r1.2
					noi di as text"# > ............ apply "as result "hwe"
					noi di as text"# > ...... binaries in "as result"`in'"
					!$plink --bfile `in' --hwe 1e-${hwep} --make-bed --out `out'
					bim2count, bim(`c(pwd)'/`out')
					noi di as text"# > ..... binaries out "as result"`out'"
					noi di as text"# > ..... SNPs in file "as result "${bim2count_snp}"
					noi di as text"# > .... indiv in file "as result "${bim2count_ind}"
					local in `out'
					}
				qui { // remove out of bound (missing)
					qui { // by-locus (1)
						local out ``function'_bim_short'-temp-r1.3
						noi di as text"# > ............ apply "as result "lmiss (1)"
						noi di as text"# > ...... binaries in "as result"`in'"
						!$plink --bfile `in' --geno ${geno1} --make-bed --out `out'
						bim2count, bim(`c(pwd)'/`out')
						noi di as text"# > ..... binaries out "as result"`out'"
						noi di as text"# > ..... SNPs in file "as result "${bim2count_snp}"
						noi di as text"# > .... indiv in file "as result "${bim2count_ind}"
						local in `out'
						}
					qui { // by-individual
						local out ``function'_bim_short'-temp-r1.4
						noi di as text"# > ............ apply "as result "imiss"
						noi di as text"# > ...... binaries in "as result"`in'"
						!$plink --bfile `in' --mind ${mind} --make-bed --out `out'
						bim2count, bim(`c(pwd)'/`out')
						noi di as text"# > ..... binaries out "as result"`out'"
						noi di as text"# > ..... SNPs in file "as result "${bim2count_snp}"
						noi di as text"# > .... indiv in file "as result "${bim2count_ind}"
						local in `out'
						}
					qui { // by-locus (2)
						local out ``function'_bim_short'-temp-r1.5
						noi di as text"# > ............ apply "as result "lmiss (2)"
						noi di as text"# > ...... binaries in "as result"`in'"
						!$plink --bfile `in' --geno ${geno2} --make-bed --out `out'
						bim2count, bim(`c(pwd)'/`out')
						noi di as text"# > ..... binaries out "as result"`out'"
						noi di as text"# > ..... SNPs in file "as result "${bim2count_snp}"
						noi di as text"# > .... indiv in file "as result "${bim2count_ind}"
						}
					}
				qui { // cryptic relatedness
					local already_created ${init_root}/bim2cryptic/``function'_bim_short'-temp-r1.5/``function'_bim_short'-temp-r1.5-bim2cryptic.remove
					noi di as text"# > ............ apply "as result "bim2cryptic.remove"
					capture confirm file `already_created'
					if !_rc {
						noi checkfile, file(`already_created')
						}
					else {	
						noi di as text"# > ...... binaries in "as result"`out'"
						bim2cryptic, bim(`c(pwd)'/`out')
						noi checkfile, file(`already_created')
						}
					import delim using `already_created',clear varn(noname)
					count
					if `r(N)' > 0 { 
						noi di as text"# > ........... remove "as result "`r(N)'" as text" individuals"
						!$plink --bfile `out' --remove `already_created' --set-hh-missing --make-bed --out ${init_root}/genotypeqc/``function'_bim_short'/``function'_bim_short'-mac5-hg19-refid-bim2frq-round1
						}
					else { 
						noi di as text"# > ........... remove "as result "0" as text" individuals"
						!$plink --bfile `out' --set-hh-missing --make-bed --out ${init_root}/genotypeqc/``function'_bim_short'/``function'_bim_short'-mac5-hg19-refid-bim2frq-round1
						}
					bim2count, bim(${init_root}/genotypeqc/``function'_bim_short'/``function'_bim_short'-mac5-hg19-refid-bim2frq-round1)
					noi di as text"# > ..... binaries out "as result"`in'-round1"
					noi di as text"# > ..... SNPs in file "as result "${bim2count_snp}"
					noi di as text"# > .... indiv in file "as result "${bim2count_ind}"	
					}
				local in ${init_root}/genotypeqc/``function'_bim_short'/``function'_bim_short'-mac5-hg19-refid-bim2frq-round1
				delete_temp_dir
				noi di as text"#########################################################################"
				noi di as text""
				}
			}
		qui { // module 2 - apply rounds 2 +
			noi di as text"#########################################################################"
			noi di as text"# SECTION - 8.2: apply quality control metrics (ROUND 2+)"
			noi di as text"#########################################################################"
			local final ${init_root}/genotypeqc/``function'_bim_short'/``function'_bim_short'-mac5-hg19-refid-bim2frq
			foreach round of num  2 / $rounds {
				clear
				set obs 1
				gen minus = `round' - 1
				gen plus  = `round' + 1
				foreach i in  plus minus {
					sum `i'
					local `i' `r(max)'
					}
				local in  `final'-round`minus'
				local out `final'-round`round'
				bim_path_to_short, bim(`in')
				local in_file ${bim_short}
				local already_created `out'.bed
				capture confirm file `already_created'
				create_temp_dir
				if !_rc {
					noi checkfile, file(`already_created')
					}
				else {	
					qui { // calculate heterozygosity
						noi di as text"# > ...... binaries in "as result"`in'"
						noi di as text"# > ............ apply "as result "bim2graphhet"	
						local already_created ${init_root}/bim2graphhet/`in_file'/`in_file'_bim2graphhet.gph
						capture confirm file `already_created'
						if !_rc {
							noi checkfile, file(`already_created')
							}
						else {	
							bim2graphhet, bim(`in')
							}
						}
					qui { // remove heterozygous out-of-bounds individuals
						local out   ``function'_bim_short'-temp-r`round'.1
						noi di as text"# > ............ apply "as result "bim2graphhet.remove"
						noi di as text"# > ...... binaries in "as result"`in'"
						local het ${init_root}/bim2graphhet/`in_file'/`in_file'_bim2graphhet.remove
						import delim using `het',clear varn(noname)
						count
						if `r(N)' > 0 { 
							noi di as text"# > ........... remove "as result "`r(N)'" as text" individuals"
							!$plink --bfile `in' --remove `het' --set-hh-missing --make-bed --out `out'
							bim2count, bim(`c(pwd)'/`out')
							noi di as text"# > ..... binaries out "as result"`out'"
							noi di as text"# > ..... SNPs in file "as result "${bim2count_snp}"
							noi di as text"# > .... indiv in file "as result "${bim2count_ind}"	
							local in `out'
							}
						else { 
							noi di as text"# > ........... remove "as result "0" as text" individuals"
							}
						}
					qui { // remove out of bound (hardy-weinberg)
						local out ``function'_bim_short'-temp-r`round'.2
						noi di as text"# > ............ apply "as result "hwe"
						noi di as text"# > ...... binaries in "as result"`in'"
						!$plink --bfile `in' --hwe 1e-${hwep} --make-bed --out `out'
						bim2count, bim(`c(pwd)'/`out')
						noi di as text"# > ..... binaries out "as result"`out'"
						noi di as text"# > ..... SNPs in file "as result "${bim2count_snp}"
						noi di as text"# > .... indiv in file "as result "${bim2count_ind}"
						local in `out'
						}
					qui { // remove out of bound (missing)
						qui { // by-locus (1)
							local out ``function'_bim_short'-temp-r`round'.3
							noi di as text"# > ............ apply "as result "lmiss (1)"
							noi di as text"# > ...... binaries in "as result"`in'"
							!$plink --bfile `in' --geno ${geno1} --make-bed --out `out'
							bim2count, bim(`c(pwd)'/`out')
							noi di as text"# > ..... binaries out "as result"`out'"
							noi di as text"# > ..... SNPs in file "as result "${bim2count_snp}"
							noi di as text"# > .... indiv in file "as result "${bim2count_ind}"
							local in `out'
							}
						qui { // by-individual
							local out ``function'_bim_short'-temp-r`round'.4
							noi di as text"# > ............ apply "as result "imiss"
							noi di as text"# > ...... binaries in "as result"`in'"
							!$plink --bfile `in' --mind ${mind} --make-bed --out `out'
							bim2count, bim(`c(pwd)'/`out')
							noi di as text"# > ..... binaries out "as result"`out'"
							noi di as text"# > ..... SNPs in file "as result "${bim2count_snp}"
							noi di as text"# > .... indiv in file "as result "${bim2count_ind}"
							local in `out'
							}
						qui { // by-locus (2)
							local out `final'-round`round'
							noi di as text"# > ............ apply "as result "lmiss (2)"
							noi di as text"# > ...... binaries in "as result"`in'"
							!$plink --bfile `in' --geno ${geno2} --make-bed --out `out'
							bim2count, bim(`out')
							noi di as text"# > ..... binaries out "as result"`out'"
							noi di as text"# > ..... SNPs in file "as result "${bim2count_snp}"
							noi di as text"# > .... indiv in file "as result "${bim2count_ind}"
							local in `out'
							}
						}	
					}		
				delete_temp_dir
				}
			noi di as text"#########################################################################"
			noi di as text""
			}
		qui { // module 4 - create quality control file
			noi di as text"#########################################################################"
			noi di as text"# SECTION - 8.3: create final quality control dataset"
			noi di as text"#########################################################################"
			!$plink --bfile `final'-round${rounds} --make-bed --out `qc'
			bim2count, bim(`qc')
			noi di as text"# > ..... binaries out "as result"`qc'"
			noi di as text"# > ..... SNPs in file "as result "${bim2count_snp}"
			noi di as text"# > .... indiv in file "as result "${bim2count_ind}"	
			noi di as text"#########################################################################"
			noi di as text""
			}
		}
	}
* 
qui { // <9> perform post-QC metrics
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 9: calculate / plot post-quality control metrics "
	noi di as text"#########################################################################"
	local qc ``function'_bim_short'-qc-``function'_version'
	local in  ${init_root}/genotypeqc/``function'_bim_short'/`qc'
	qui { // frq
		local graph bim2graphfrq
		noi di as text"# > ............ apply "as result "`graph'"
		local already_created ${init_root}/`graph'/`qc'/`qc'_`graph'.gph
		capture confirm file `already_created'
		if !_rc {
			noi checkfile, file(`already_created')
			}
		else {	
			`graph', bim(`in')
			noi checkfile, file(`already_created')
			}
		}
	qui { // het
		local graph bim2graphhet
		noi di as text"# > ............ apply "as result "`graph'"
		local already_created ${init_root}/`graph'/`qc'/`qc'_`graph'.gph
		capture confirm file `already_created'
		if !_rc {
			noi checkfile, file(`already_created')
			}
		else {	
			`graph', bim(`in')
			noi checkfile, file(`already_created')
			}
		}
	qui { // hwe
		local graph bim2graphhwe
		noi di as text"# > ............ apply "as result "`graph'"
		local already_created ${init_root}/`graph'/`qc'/`qc'_`graph'.gph
		capture confirm file `already_created'
		if !_rc {
			noi checkfile, file(`already_created')
			}
		else {	
			`graph', bim(`in')
			noi checkfile, file(`already_created')
			}
		}
	qui { // kin0	
		local graph  bim2graphkin0
		qui { // hist
			local graph2 hist
			noi di as text"# > ............ apply "as result "`graph'"
			local already_created ${init_root}/`graph'/`qc'/`qc'_`graph'-`graph2'.gph
			capture confirm file `already_created'
			if !_rc {
				noi checkfile, file(`already_created')
				}
			else {	
				`graph', bim(`in')
				noi checkfile, file(`already_created')
				}
			}
		qui { // scatter
			local graph2 scatter
			noi di as text"# > ............ apply "as result "`graph'"
			local already_created ${init_root}/`graph'/`qc'/`qc'_`graph'-`graph2'.gph
			capture confirm file `already_created'
			if !_rc {
				noi checkfile, file(`already_created')
				}
			else {	
				`graph', bim(`in')
				noi checkfile, file(`already_created')
				}
			}
		}
	qui { // miss	
		local graph  bim2graphmiss
		qui { // lmiss
			local graph2 lmiss
			noi di as text"# > ............ apply "as result "`graph'"
			local already_created ${init_root}/`graph'/`qc'/`qc'_`graph'-`graph2'.gph
			capture confirm file `already_created'
			if !_rc {
				noi checkfile, file(`already_created')
				}
			else {	
				`graph', bim(`in')
				noi checkfile, file(`already_created')
				}
			}
		qui { // imiss
			local graph2 imiss
			noi di as text"# > ............ apply "as result "`graph'"
			local already_created ${init_root}/`graph'/`qc'/`qc'_`graph'-`graph2'.gph
			capture confirm file `already_created'
			if !_rc {
				noi checkfile, file(`already_created')
				}
			else {	
				`graph', bim(`in')
				noi checkfile, file(`already_created')
				}
			}
		}
	noi di as text"#########################################################################"
	noi di as text""
	}
*
qui { // <10> determine ancestry metrics
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 10 : define eigenvec / ancestry plots "
	noi di as text"#########################################################################"
	local qc ``function'_bim_short'-qc-``function'_version'
	local in  ${init_root}/genotypeqc/``function'_bim_short'/`qc'
	noi di as text"# > ............ apply "as result "bim2ancestry"	
	local already_created ${init_root}/bim2ancestry/`qc'/`qc'-bim2ancestry-merged_population.dta
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		}
	else {	
		noi bim2ancestry, bim(`in') ref(${bim2ancestry_ref})
		noi checkfile, file(`already_created')
		}
	noi di as text"#########################################################################"
	noi di as text""
	}
*
qui { // <11> write log file
	bim2count, bim(`bim')
	local bim2count_snp_bim ${bim2count_snp}
	local bim2count_ind_bim ${bim2count_ind}
	bim2count, bim(${init_root}/genotypeqc/``function'_bim_short'/`qc')
	local bim2count_snp_qc  ${bim2count_snp}
	local bim2count_ind_qc  ${bim2count_ind}	
	file open myfile using "${init_root}/genotypeqc/``function'_bim_short'/``function'_bim_short'-`function'.log", write replace
	file write myfile`"#########################################################################"' _n
	file write myfile`"# > .................. `function'"' _n
	file write myfile`"# > ........... author Richard Anney"' _n
	file write myfile`"# > ... report bugs to anneyr@cardiff.ac.uk"' _n
	file write myfile`"# >  code available at https://github.com/ricanney"' _n
	file write myfile`"# > .......... license CC BY-NC (CC Attribution-NonCommercial)"' _n
	file write myfile`"# > .................. all functions considered beta "caveat emptor""' _n
	file write myfile`"#########################################################################"' _n
	file write myfile`"# Started:             $S_DATE $S_TIME "' _n
	file write myfile`"# Username:            `c(username)'"' _n
	file write myfile`"# Operating System:    `c(os)'"' _n
	file write myfile`"#########################################################################"' _n
	file write myfile `"# > ..... project name ``function'_bim_short'"' _n
	file write myfile `"# > .... plink version ${plink}"' _n
	file write myfile `"# > .... plink version ${plink2}"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > .... input dataset ``function'_bim_short'"' _n
	file write myfile `"# > ..... SNPs in file `bim2count_snp_bim'"' _n
	file write myfile `"# > ... indiv. in file `bim2count_ind_bim'"' _n
	file write myfile `"# > . array defined as `bim2array_array'"' _n
	file write myfile `"# > .... jaccard index `bim2array_jaccard'"' _n
	import delim using  ${init_root}/bim2build/``function'_bim_short'-mac5/``function'_bim_short'-mac5.build, clear varnames(1) case(preserve)
	gen script = "global bim2build " + build
	outsheet script using temp.do, non noq replace
	do temp.do
	erase temp.do
	file write myfile `"# > ......... build is ${bim2build}"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > ... output dataset ``function'_bim_short'-qc-``function'_version'"' _n
	file write myfile `"# > ..... SNPs in file `bim2count_snp_qc'"' _n
	file write myfile `"# > ... indiv. in file `bim2count_ind_qc'"' _n
	file write myfile `"# > ......... build is hg19 +1"' _n
	file write myfile `"# > ..... SNP names as ``function'_bim2refid_short'"' _n
	clear
	set obs 1
	gen admixed = "`admixed'"
	if admixed == "" {
		file write myfile `"# >  frq restricted to ``function'_bim2frq_short'"' _n
		}
	else {
		file write myfile `"# > frq not restricted *admixed sample*"' _n
		}
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > quality control paramenters"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > ..... rounds of qc $rounds"' _n
	file write myfile `"# > ............ --mac 5"' _n
	file write myfile `"# > ........... --geno ${geno1};${geno2}"' _n
	file write myfile `"# > ........... --mind ${mind}"' _n
	file write myfile `"# > .......... --hardy 1e-${hwep}"' _n
	file write myfile `"# > ... heterozygosity ${hetsd} standard deviations from mean"' _n
	file write myfile `"# > .......... kinship "' _n
	file write myfile `"# > ....... duplicates ${kin_d}"' _n
	file write myfile `"# > ....... 1st-degree ${kin_f}"' _n
	file write myfile `"# > ....... 2nd-degree ${kin_s}"' _n
	file write myfile `"# > ....... 3rd-degree ${kin_t}"' _n
	file write myfile `"# > ancestry mapped to ``function'_bim2ancestry_short'"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > overview of `function'"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > 1. import bim files"' _n
	file write myfile `"# > 2. define array using bim2array"' _n
	file write myfile `"# > 3. limit to mac 5"' _n
	file write myfile `"# > 4. define / update genome build using bim2build"' _n
	file write myfile `"# > 5. update marker identifiers bim2refid"' _n
	file write myfile `"# > 6. prune out-of-bounds allele frequencies using bim2frq_compare"' _n
	file write myfile `"# > 7. define pre-qc metrics "' _n
	file write myfile `"# > 8. apply quality control "' _n
	file write myfile `"# > 9. define post-qc metrics "' _n
	file write myfile `"# > 10. define ancestry metrics using bim2ancestry "' _n
	file write myfile `"#########################################################################"' _n
	file close myfile	
	}
*
qui { // <12> clean folder
	cd ${init_root}/genotypeqc/``function'_bim_short'/
	files2dta, dir(${init_root}/genotypeqc/``function'_bim_short'/)
	gen keep = .
	replace keep = 1 if file == "``function'_bim_short'-qc-``function'_version'.bed"
	replace keep = 1 if file == "``function'_bim_short'-qc-``function'_version'.bim"
	replace keep = 1 if file == "``function'_bim_short'-qc-``function'_version'.fam"
	replace keep = 1 if file == "``function'_bim_short'-qc-``function'_version'.log"
	replace keep = 1 if file == "``function'_bim_short'.bed"
	replace keep = 1 if file == "``function'_bim_short'.bim"
	replace keep = 1 if file == "``function'_bim_short'.fam"
	replace keep = 1 if file == "``function'_bim_short'.log"
	replace keep = 1 if file == "``function'_bim_short'-`function'.log"
	drop if keep == 1
	gen script = "erase " + file
	outsheet script using temp.do, non noq replace 
	do temp.do
	erase temp.do
	clear
	}
*
qui cd `return'
restore

noi di as text"#########################################################################"
noi di as text"# Completed: $S_DATE $S_TIME"
noi di as text"#########################################################################"
end;
