*! 1.0.1 Richard Anney 22jan2020

program  summaryqc2gene
local    function summaryqc2gene
local    `function'_version v1

syntax , summaryqc(string asis) ldref(string asis) glist_hg19(string asis)

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
	ref_path_to_short, ref(`summaryqc')
	local `function'_summaryqc_short ${ref_short}
	cd ${init_root}	
	!mkdir `function'
	cd `function'
	noi di as text"# > ..... project name "as result"``function'_summaryqc_short'"
	!mkdir ``function'_summaryqc_short'
	cd ``function'_summaryqc_short'
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 2 - define files
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 2: define files"
	noi di as text"#########################################################################"
	noi checkfile, file(${plink})
	noi checkfile, file(${plink2})
	checktabbed
	noi di as text"# > ......... locating files to process"
	noi checkfile, file(`summaryqc'-summaryqc.dta)
	foreach i in _3prime_overlapping_ncrna IG_C_gene IG_C_pseudogene IG_D_gene IG_J_gene IG_J_pseudogene IG_V_gene IG_V_pseudogene Mt_rRNA Mt_tRNA TR_C_gene TR_D_gene TR_J_gene TR_J_pseudogene TR_V_gene antisense lincRNA miRNA misc_RNA polymorphic_pseudogene processed_transcript protein_coding pseudogene rRNA sense_intronic sense_overlapping snRNA snoRNA {
		noi checkfile, file(`glist_hg19'/Homo_sapiens.GRCh37.87.gtf-`i'.hg-list)
		}
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 3 - identify genes directly hit by SNPs
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3: map by SNP"
	noi di as text"#########################################################################"
	noi di as text"# > ..... loading from "as result"``function'_summaryqc_short'-summaryqc2top.dta"
	local already_created ${init_root}/summaryqc2top/``function'_summaryqc_short'/``function'_summaryqc_short'-summaryqc2top.dta
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		}
	else {	
		noi summaryqc2top, summaryqc(`summaryqc') ldref(`ldref') glist_hg19(`glist_hg19')
		noi checkfile, file(`already_created')
		}
	noi di as text"# > define LD pairs of top hits r^2"
	local already_created ${init_root}/`function'/``function'_summaryqc_short'/``function'_summaryqc_short'-ld-pairs.dta
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		}
	else {
		use 	${init_root}/summaryqc2top/``function'_summaryqc_short'/``function'_summaryqc_short'-summaryqc2top.dta, clear
		keep SNP
		outsheet using ``function'_summaryqc_short'.snp-in, non noq replace
		noi di as text"# > mapping LD-pairs using plink"
		!plink2 --bfile `ldref' --r2 --ld-window 100000 --ld-window-kb 100000 --ld-snp-list ``function'_summaryqc_short'.snp-in --ld-window-r2 0.8 --out ``function'_summaryqc_short'.snp-out
		!$tabbed ``function'_summaryqc_short'.snp-out.ld
		import delim using ``function'_summaryqc_short'.snp-out.ld.tabbed, clear varnames(1)
		rename (snp_a snp_b) (SNP PROXY_SNP)
		for var chr_b bp_b: tostring X, replace
		gen PROXY_POSITION = "chr" + chr_b + ":" + bp_b
		rename r2 R2
		keep SNP PROXY_SNP PROXY_POSITION R2
		save ${init_root}/`function'/``function'_summaryqc_short'/``function'_summaryqc_short'-ld-pairs.dta, replace
		noi checkfile, file(`already_created')
		}
	local already_created ${init_root}/`function'/``function'_summaryqc_short'/``function'_summaryqc_short'-`function'-snp-hit.dta
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		}
	else {
		use ${init_root}/`function'/``function'_summaryqc_short'/``function'_summaryqc_short'-ld-pairs.dta, clear
		merge m:1 SNP using ${init_root}/summaryqc2top/``function'_summaryqc_short'/``function'_summaryqc_short'-summaryqc2top.dta
		split PROXY_POSITION, p("chr"":")
		rename (PROXY_POSITION2 PROXY_POSITION3)(CHR BP)
		keep CHR SNP BP P
		order CHR SNP BP P
		outsheet using ``function'_summaryqc_short'-summaryqc2top.top, noq replace
		foreach i in _3prime_overlapping_ncrna IG_C_gene IG_C_pseudogene IG_D_gene IG_J_gene IG_J_pseudogene IG_V_gene IG_V_pseudogene Mt_rRNA Mt_tRNA TR_C_gene TR_D_gene TR_J_gene TR_J_pseudogene TR_V_gene antisense lincRNA miRNA misc_RNA polymorphic_pseudogene processed_transcript protein_coding pseudogene rRNA sense_intronic sense_overlapping snRNA snoRNA {
			!${plink} --annotate ``function'_summaryqc_short'-summaryqc2top.top ranges=`glist_hg19'/Homo_sapiens.GRCh37.87.gtf-`i'.hg-list minimal --out ``function'_summaryqc_short'-summaryqc2top.`i'
			}
		foreach i in _3prime_overlapping_ncrna IG_C_gene IG_C_pseudogene IG_D_gene IG_J_gene IG_J_pseudogene IG_V_gene IG_V_pseudogene Mt_rRNA Mt_tRNA TR_C_gene TR_D_gene TR_J_gene TR_J_pseudogene TR_V_gene antisense lincRNA miRNA misc_RNA polymorphic_pseudogene processed_transcript protein_coding pseudogene rRNA sense_intronic sense_overlapping snRNA snoRNA {
			import delim using ``function'_summaryqc_short'-summaryqc2top.`i'.annot, clear varnames(1) 
			keep snp pannot
			split pannot, p(" ""|")
			drop pannot1
			gen pannot99 = ""
			keep snp pannot2-pannot99
			duplicates drop
			drop pannot99
			egen x = seq(),by(snp)
			reshape long pannot, i(snp x) 
			keep snp pannot
			replace pannot = "." if pannot == ""
			duplicates drop
			sort pannot
			egen x = seq(), by(snp)
			reshape wide pannot, i(snp) j(x) 
			gen `i' = pannot1
			replace `i' = "" if `i' == "."
			capture confirm var pannot2 
			if !_rc {
				}
			else {
				gen pannot2 = ""
				}
			gen pannot99 = ""
			aorder
			for var pannot2 - pannot99: replace `i' = `i' + " " + X
			keep snp `i'
			replace `i' = "." if `i' == ""
			save ``function'_summaryqc_short'-summaryqc2top.`i'.annot.dta, replace
			}
		use ${init_root}/summaryqc2top/``function'_summaryqc_short'/``function'_summaryqc_short'-summaryqc2top.dta, clear
		rename SNP snp
		foreach i in _3prime_overlapping_ncrna IG_C_gene IG_C_pseudogene IG_D_gene IG_J_gene IG_J_pseudogene IG_V_gene IG_V_pseudogene Mt_rRNA Mt_tRNA TR_C_gene TR_D_gene TR_J_gene TR_J_pseudogene TR_V_gene antisense lincRNA miRNA misc_RNA polymorphic_pseudogene processed_transcript protein_coding pseudogene rRNA sense_intronic sense_overlapping snRNA snoRNA {
			merge 1:1 snp using ``function'_summaryqc_short'-summaryqc2top.`i'.annot.dta
			drop _m
			}
		sort P
		rename  snp SNP
		keep  SNP POSITION A1 A2 OR P _3prime_overlapping_ncrna IG_C_gene IG_C_pseudogene IG_D_gene IG_J_gene IG_J_pseudogene IG_V_gene IG_V_pseudogene Mt_rRNA Mt_tRNA TR_C_gene TR_D_gene TR_J_gene TR_J_pseudogene TR_V_gene antisense lincRNA miRNA misc_RNA polymorphic_pseudogene processed_transcript protein_coding pseudogene rRNA sense_intronic sense_overlapping snRNA snoRNA
		order SNP POSITION A1 A2 OR P _3prime_overlapping_ncrna IG_C_gene IG_C_pseudogene IG_D_gene IG_J_gene IG_J_pseudogene IG_V_gene IG_V_pseudogene Mt_rRNA Mt_tRNA TR_C_gene TR_D_gene TR_J_gene TR_J_pseudogene TR_V_gene antisense lincRNA miRNA misc_RNA polymorphic_pseudogene processed_transcript protein_coding pseudogene rRNA sense_intronic sense_overlapping snRNA snoRNA
		save ``function'_summaryqc_short'-`function'-snp-hit.dta,replace
		capture confirm file `already_created'
		}
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 4 - identify genes directly hit by ranges
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 4: map by regions"
	noi di as text"#########################################################################"
	noi di as text"# > ..... loading from "as result"``function'_summaryqc_short'-summaryqc2top.dta"
	local already_created ${init_root}/summaryqc2top/``function'_summaryqc_short'/``function'_summaryqc_short'-summaryqc2top.dta
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		}
	else {	
		noi summaryqc2top, summaryqc(`summaryqc') ldref(`ldref') glist_hg19(`glist_hg19')
		noi checkfile, file(`already_created')
		}
	local already_created ${init_root}/`function'/``function'_summaryqc_short'/``function'_summaryqc_short'-`function'-range-hit.dta
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		}
	else {
		use ${init_root}/summaryqc2top/``function'_summaryqc_short'/``function'_summaryqc_short'-summaryqc2top.dta, clear
		split LD_RANGE, p("chr"":""..")
		rename (LD_RANGE2 LD_RANGE3 LD_RANGE4) (CHR START END)
		for var CHR START END: destring X, replace
		sort CHR START END
		for var CHR START END: tostring X, replace
		replace CHR = "chr" + CHR
		keep SNP CHR START END
		order CHR START END SNP
		outsheet using ``function'_summaryqc_short'-summaryqc2top.bed, non noq replace
		foreach i in _3prime_overlapping_ncrna IG_C_gene IG_C_pseudogene IG_D_gene IG_J_gene IG_J_pseudogene IG_V_gene IG_V_pseudogene Mt_rRNA Mt_tRNA TR_C_gene TR_D_gene TR_J_gene TR_J_pseudogene TR_V_gene antisense lincRNA miRNA misc_RNA polymorphic_pseudogene processed_transcript protein_coding pseudogene rRNA sense_intronic sense_overlapping snRNA snoRNA {
			noi !bedtools intersect -a ``function'_summaryqc_short'-summaryqc2top.bed -b `glist_hg19'/Homo_sapiens.GRCh37.87.gtf-`i'.hg-list-35+5.bed -wa -wb > ``function'_summaryqc_short'-by-Homo_sapiens.GRCh37.87.gtf-`i'.hg-list-35+5.out
			import delim using ``function'_summaryqc_short'-by-Homo_sapiens.GRCh37.87.gtf-`i'.hg-list-35+5.out, clear varnames(nonames)
			count
			if `r(N)' == 0 {
				clear
				set obs 1
				gen SNP = ""
				gen `i' = ""
				}
			else {
				keep v4 v8
				rename (v4 v8) (SNP `i')
				duplicates drop
				sort SNP `i'
				egen x = seq(),by(SNP)
				reshape wide `i', i(SNP) j(x)	
				gen `i'999 = ""
				gen `i' = ""
				for var `i'1-`i'999: replace `i' = `i' + " " +  X
				foreach j of num 1/999 {
					replace `i' = subinstr(`i',"  "," ",.)
					}
				keep SNP `i'	
				}
			save ``function'_summaryqc_short'-summaryqc2top.`i'.annot.dta, replace
			}
		use ${init_root}/summaryqc2top/``function'_summaryqc_short'/``function'_summaryqc_short'-summaryqc2top.dta, clear
		foreach i in _3prime_overlapping_ncrna IG_C_gene IG_C_pseudogene IG_D_gene IG_J_gene IG_J_pseudogene IG_V_gene IG_V_pseudogene Mt_rRNA Mt_tRNA TR_C_gene TR_D_gene TR_J_gene TR_J_pseudogene TR_V_gene antisense lincRNA miRNA misc_RNA polymorphic_pseudogene processed_transcript protein_coding pseudogene rRNA sense_intronic sense_overlapping snRNA snoRNA {
			merge 1:1 SNP using ``function'_summaryqc_short'-summaryqc2top.`i'.annot.dta
			drop _m
			}
		sort P
		keep  SNP LD_RANGE A1 A2 OR P _3prime_overlapping_ncrna IG_C_gene IG_C_pseudogene IG_D_gene IG_J_gene IG_J_pseudogene IG_V_gene IG_V_pseudogene Mt_rRNA Mt_tRNA TR_C_gene TR_D_gene TR_J_gene TR_J_pseudogene TR_V_gene antisense lincRNA miRNA misc_RNA polymorphic_pseudogene processed_transcript protein_coding pseudogene rRNA sense_intronic sense_overlapping snRNA snoRNA
		order SNP LD_RANGE A1 A2 OR P _3prime_overlapping_ncrna IG_C_gene IG_C_pseudogene IG_D_gene IG_J_gene IG_J_pseudogene IG_V_gene IG_V_pseudogene Mt_rRNA Mt_tRNA TR_C_gene TR_D_gene TR_J_gene TR_J_pseudogene TR_V_gene antisense lincRNA miRNA misc_RNA polymorphic_pseudogene processed_transcript protein_coding pseudogene rRNA sense_intronic sense_overlapping snRNA snoRNA
		save ``function'_summaryqc_short'-`function'-range-hit.dta,replace
		capture confirm file `already_created'	
		}
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 4 - write log file
	file open myfile using "``function'_summaryqc_short'-`function'.log", write replace
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
	file write myfile `"# > ....... summaryqc ``function'_summaryqc_short'"' _n
	file write myfile `"# > ....... reference ``function'_ldref_short'"' _n
	file write myfile `"# > ...... saved data ``function'_summaryqc_short'-`function'.dta"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > overview of `function'"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > 1. identify top N SNPs"' _n
	file write myfile `"# > 2. extend LD to capture proxy_SNPs"' _n
	file write myfile `"# > 3. map genes"' _n
	file write myfile `"# > 4. map ranges"' _n
	file write myfile `"#########################################################################"' _n
	file close myfile	
	}
qui { // module 5 - clean folder
	files2dta, dir(`c(pwd)')
	gen keep = .
	replace keep = 1 if file == "``function'_summaryqc_short'-`function'.log"
	replace keep = 1 if file == "``function'_summaryqc_short'-`function'-snp-hit.dta"
	replace keep = 1 if file == "``function'_summaryqc_short'-`function'-range-hit.dta"
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

   
