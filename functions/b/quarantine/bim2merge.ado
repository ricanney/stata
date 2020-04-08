*! 1.0.2 Richard Anney 19dec2018
* - all samples should be processed (genotypeqc) against common reference
* - limit to merging 2 datasets only

program  bim2merge
local    function bim2merge
local    `function'_version v2
syntax , bim1(string asis) bim2(string asis) 

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
	bim_path_to_short, bim(`bim1')
	local `function'_bim1_short ${bim_short}
	bim_path_to_short, bim(`bim2')
	local `function'_bim2_short ${bim_short}
	cd ${init_root}	
	!mkdir `function'
	cd `function'
	local `function'_project ``function'_bim1_short'-x-``function'_bim2_short'
	noi di as text"# > ..... project name "as result"``function'_project'"
	!mkdir ``function'_project'
	cd ``function'_project'
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 2 - define files
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 2: define files"
	noi di as text"#########################################################################"
	foreach file in bim bed fam {
		noi checkfile, file(`bim1'.`file')
		}
	foreach file in bim bed fam {
		noi checkfile, file(`bim2'.`file')
		}
	noi checkfile, file(${plink})
	checktabbed
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 3 - define common SNPs for merge
	noi di as text""	
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3: define common SNPs for merge"
	noi di as text"#########################################################################"
	foreach i in 1 2 {
		local already_created ${init_root}/bim2dta/``function'_bim`i'_short'/``function'_bim`i'_short'_bim.dta
		capture confirm file `already_created'
		if !_rc {
				noi checkfile, file(`already_created')
				}
		else {
				noi checkfile, file(`already_created')
				noi di as text"# > .................. create "as result "``function'_bim`i'_short'_bim.dta"
				noi bim2dta, bim(`bim`i'')
				}		
		}
	qui { // identify overlap
		noi di as text"# > .................. identify overlap"
		use loc_name using       ${init_root}/bim2dta/``function'_bim1_short'/``function'_bim1_short'_bim.dta, clear
		count
		local `function'_bim1N `r(N)'
		noi di as text"# > ..... SNPs in bim1 "as result"``function'_bim1N'"
		merge 1:1 loc_name using ${init_root}/bim2dta/``function'_bim2_short'/``function'_bim2_short'_bim.dta
		count if _m != 1
		local `function'_bim2N `r(N)'
		noi di as text"# > ..... SNPs in bim2 "as result"``function'_bim2N'"
		keep if _m == 3
		count
		local `function'_mergeN `r(N)'
		noi di as text"# > .... SNPs in merge "as result"``function'_mergeN'"
		}
	qui { // drop ambiguous markers
			noi di as text"# > ................. drop ambiguous SNPs"
			gen ambig = .
			replace ambig = 1 if gt == "S"
			replace ambig = 1 if gt == "W"
			count if ambig == 1
			local `function'_ambigN `r(N)'
			noi di as text"# > ............ drop "as result"``function'_ambigN'"as text" ambiguous SNPs"
			drop if ambig == 1
			drop ambig
			}
	qui { // define merge snps
		outsheet snp using ``function'_project'.extract, non noq replace
		}
	noi di as text"#########################################################################"
		}
qui { // module 4 - process and merge
		noi di as text""	
		noi di as text"#########################################################################"
		noi di as text"# SECTION - 4: process files"
		noi di as text"#########################################################################"
		qui { // process PLINK binaries
			qui { // limit to common SNPs
				noi di as text"# > .................. extracting common SNPs"
				!$plink --bfile `bim1' --make-bed --extract ``function'_project'.extract   --out ``function'_bim1_short'
				!$plink --bfile `bim2' --make-bed --extract ``function'_project'.extract   --out ``function'_bim2_short'
				}
			qui { // flip markers
				noi di as text"# > .................. identifying flip errors"
				!$plink --bfile ``function'_bim1_short' --bmerge ``function'_bim2_short' --make-bed --out ``function'_project'
				capture confirm file ``function'_project'-merge.missnp
				if !_rc {
					noi di as text"# > .................. fixing flip errors"
					!$plink --bfile ``function'_bim1_short' --flip ``function'_project'-merge.missnp --make-bed --out ``function'_bim1_short'-flip
					noi di as text"# > .................. merging datasets"
					!$plink --bfile ``function'_bim1_short'-flip --bmerge ``function'_bim2_short' --make-bed --out ``function'_project'
					}
				else {
					noi di as text"# > > .................. no flip errors / datasets merged"
					}
				noi di as text"#########################################################################"
				}
			}
		qui { // count snps and individuals in sets
			bim2count, bim(`bim1')
			local `function'_bim1snp ${bim2count_snp}
			local `function'_bim1ind ${bim2count_ind}
			bim2count, bim(`bim2')
			local `function'_bim2snp ${bim2count_snp}
			local `function'_bim2ind ${bim2count_ind}	
			bim2count, bim(${init_root}/bim2merge/``function'_project'/``function'_project')
			local `function'_bimpsnp ${bim2count_snp}
			local `function'_bimpind ${bim2count_ind}
			}
		}	
qui { // module 5 - write log file
	file open myfile using "``function'_project'-`function'.log", write replace
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
	file write myfile `"# > ..... project name ``function'_project'"' _n
	file write myfile `"# > .... plink version ${plink}"' _n
	file write myfile `"# > ............. bim1 ``function'_bim1_short'"' _n
	file write myfile `"# > ..... SNPs in bim1 ``function'_bim1snp'"' _n
	file write myfile `"# > ... indiv. in bim1 ``function'_bim1ind'"' _n
	file write myfile `"# > ............. bim2 ``function'_bim1_short'"' _n
	file write myfile `"# > ..... SNPs in bim2 ``function'_bim2snp'"' _n
	file write myfile `"# > ... indiv. in bim2 ``function'_bim2ind'"' _n
	file write myfile `"# > .. SNPs post merge ``function'_mergeN'"' _n
	file write myfile `"# > ............. drop ``function'_ambigN' ambiguous SNPs"' _n
	file write myfile `"# > ... SNPs in merged ``function'_bimpsnp'"' _n
	file write myfile `"# > . indiv. in merged ``function'_bimpind'"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > overview of `function'"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > 1. create/use and _bim.dta files from bim1 and bim2"' _n
	file write myfile `"# > 2. merge using loc_name to remove incompatibilities"' _n
	file write myfile `"# > 3. drop ambiguous SNPs"' _n
	file write myfile `"# > 5. merge and check for flip errors / re-merge"' _n
	file write myfile `"#########################################################################"' _n
	file close myfile	
	}
qui { // module 6 - clean folder
	files2dta, dir(`c(pwd)')
	gen keep = .
	replace keep = 1 if file == "``function'_project'.bim"
	replace keep = 1 if file == "``function'_project'.bed"
	replace keep = 1 if file == "``function'_project'.fam"
	replace keep = 1 if file == "``function'_project'.log"
	replace keep = 1 if file == "``function'_project'-`function'.log"
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
	
