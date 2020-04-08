*! 1.0.1 Richard Anney 29jan2020

program  summaryqc2annot
local    function summaryqc2annot
local    `function'_version v1

syntax , summaryqc(string asis) ldref(string asis) generef(string asis) eqtl(string asis) cadd(string asis)

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
	checktabbed
	noi di as text"# > ......... locating files to process"
	noi checkfile, file(`summaryqc'-summaryqc.dta)
	noi checkfile, file(${init_root}/summaryqc2top/``function'_summaryqc_short'/``function'_summaryqc_short'-summaryqc2top.dta)
	qui { // cheack ldref
		foreach file in bim bed fam {
			noi checkfile, file(`ldref'.`file')
			}
		}	
	qui { // check gene 
		foreach i in _3prime_overlapping_ncrna IG_C_gene IG_C_pseudogene IG_D_gene IG_J_gene IG_J_pseudogene IG_V_gene IG_V_pseudogene Mt_rRNA Mt_tRNA TR_C_gene TR_D_gene TR_J_gene TR_J_pseudogene TR_V_gene antisense lincRNA miRNA misc_RNA polymorphic_pseudogene processed_transcript protein_coding pseudogene rRNA sense_intronic sense_overlapping snRNA snoRNA {
			noi checkfile, file(`generef'/Homo_sapiens.GRCh37.87.gtf-`i'.hg-list)
			noi checkfile, file(`generef'/Homo_sapiens.GRCh37.87.gtf-`i'.hg-list-35+5.bed)
			}
		}
	qui { // check eqtl
		foreach i in adipose adrenal blood brain breast cardiovascular gut liver lung muscle nerve ovary pancreas pituitary prostate salivary skin spleen testis thyroid uterus vagina {
			noi checkfile, file(`eqtl'/GTEx_Analysis_v7_eQTL_`i'_1e5.dta)
			}
		}
	qui { // check cadd 
		foreach i in AnnoType CADDv1.3 ChromHMM Cons ENCODE fitCons mirSVR PolyPhen-SIFT protein segway TFBSmotif {
			noi checkfile, file(`cadd'/eur-1000g-phase3-chrall-mac5_`i'.dta)
			}
		}
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 3 - check summaryqc2top
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3: define LD partners for INDEX SNPs"
	noi di as text"#########################################################################"
	noi di as text"# > ......... checking "as result"``function'_summaryqc_short'-summaryqc2top.dta"
	local already_created ${init_root}/summaryqc2top/``function'_summaryqc_short'/``function'_summaryqc_short'-summaryqc2top.dta
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		}
	else {	
		noi summaryqc2top, summaryqc(`summaryqc') ldref(`ldref') 
		noi checkfile, file(`already_created')
		}
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 4 - define LD partners for INDEX SNPs
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 4: define LD partners for INDEX SNPs"
	noi di as text"#########################################################################"
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
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 5 - map genes (by SNP)
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 5: map genes (by SNP)"
	noi di as text"#########################################################################"
	local already_created ${init_root}/`function'/``function'_summaryqc_short'/``function'_summaryqc_short'-`function'-gene-by-snp.dta
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
			!${plink} --annotate ``function'_summaryqc_short'-summaryqc2top.top ranges=`generef'/Homo_sapiens.GRCh37.87.gtf-`i'.hg-list minimal --out ``function'_summaryqc_short'-summaryqc2top.`i'
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
		save ``function'_summaryqc_short'-`function'-gene-by-snp.dta,replace
		capture confirm file `already_created'
		}
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 6 - map genes (by LD Range)
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 6: map genes (by LD Range)"
	noi di as text"#########################################################################"
	local already_created ${init_root}/`function'/``function'_summaryqc_short'/``function'_summaryqc_short'-`function'-gene-by-range.dta
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
			noi !bedtools intersect -a ``function'_summaryqc_short'-summaryqc2top.bed -b `generef'/Homo_sapiens.GRCh37.87.gtf-`i'.hg-list-35+5.bed -wa -wb > ``function'_summaryqc_short'-by-Homo_sapiens.GRCh37.87.gtf-`i'.hg-list-35+5.out
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
		save ``function'_summaryqc_short'-`function'-gene-by-range.dta,replace
		capture confirm file `already_created'	
		}
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 7 - map eQTL (by SNP)
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 7: map eQTL (by SNP)"
	noi di as text"#########################################################################"
	local already_created ${init_root}/`function'/``function'_summaryqc_short'/``function'_summaryqc_short'-`function'-eqtl-by-snp-verbose.dta
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		}
	else {
		use `eqtl'/GTEx_Analysis_v7_eQTL_all_1e5.dta, clear
		renvars, upper
		rename SNP PROXY_SNP
		sort PVAL
		egen x = seq(),by(PROXY_SNP SYMBOL TISSUE SUPER_TISSUE)
		gen str5 X= string(BETA,"%06.5f")
		gen str4 Y= string(SE,"%05.4f")
		gen Z_= string(PVAL) 
		drop BETA
		gen BETA = X + " (" + Y + "); p = " + Z
		keep PROXY_SNP SYMBOL TISSUE SUPER_TISSUE BETA PVAL
		merge m:m PROXY_SNP using ``function'_summaryqc_short'-ld-pairs.dta
		keep if _m == 3
		sort PVAL
		egen x = seq(),by(SNP SYMBOL SUPER TISSUE)
		keep if x == 1
		egen y = seq(),by(SNP SYMBOL SUPER )
		gen str4 X= string(R2,"%05.4f")
		gen TOP = PROXY_SNP + " (" + BETA + "; R2 = " + X + ")"
		sort SNP SYMBOL SUPER TISSUE
		keep SNP SYMBOL SUPER TISSUE TOP y 
		order SNP SYMBOL SUPER TISSUE TOP 
		gen TISSUE_EQTL = SYMBOL + " [" + TISSUE + "; " + TOP + "]" 
		gen SUPER_EQTL = SYMBOL + " [" + SUPER + "; " + TOP + "]" 
		replace SUPER_EQTL = "" if y != 1
		keep SNP SYMBOL TISSUE_EQTL SUPER_EQTL
		save ``function'_summaryqc_short'-`function'-eqtl-by-snp-verbose.dta,replace
		noi checkfile, file(`already_created')
		}
	local already_created ${init_root}/`function'/``function'_summaryqc_short'/``function'_summaryqc_short'-`function'-eqtl-by-snp-compact.dta
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		}
	else {
		use ``function'_summaryqc_short'-`function'-eqtl-by-snp-verbose.dta, clear
		keep SNP SYMBOL SUPER_EQTL
		split SUPER, p(" ["";")
		keep SNP SYMBOL SUPER_EQTL2
		rename SUPER_EQTL2 TISSUE
		drop if TISSUE == ""
		sort TISSUE
		egen n = seq(),by(SNP SYMBOL)
		reshape wide TISSUE, i(SNP SYMBOL) j(n)
		gen TISSUE999 = ""
		gen TISSUE = ""
		for var TISSUE1 - TISSUE999: replace TISSUE = TISSUE + " " + X if X != ""
		keep SNP SYMBOL TISSUE
		gen EQTL = SYMBOL + "(" + TISSUE + ")"
		replace EQTL = subinstr(EQTL,"( ","(",.)
		keep SNP EQTL
		egen n = seq(),by(SNP )
		reshape wide EQTL, i(SNP ) j(n)
		gen EQTL999 = ""
		gen EQTL = ""
		for var EQTL1 - EQTL999: replace EQTL = EQTL + " " + X if X != ""
		keep SNP EQTL
		save ``function'_summaryqc_short'-`function'-eqtl-by-snp-compact.dta, replace
		noi checkfile, file(`already_created')
		}
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 8 - map CADD (by SNP)
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 8: map CADD (by SNP)"
	noi di as text"#########################################################################"
	qui { // map max consequence
		local already_created ${init_root}/`function'/``function'_summaryqc_short'/``function'_summaryqc_short'-`function'-AnnoType-by-snp.dta
		capture confirm file `already_created'		
		if !_rc {
			noi checkfile, file(`already_created')
			}
		else {
			use `cadd'/eur-1000g-phase3-chrall-mac5_AnnoType.dta, clear
			renvars, upper
			rename SNP PROXY_SNP
			merge m:m PROXY_SNP using ``function'_summaryqc_short'-ld-pairs.dta
			keep if _m == 3
			keep SNP CONSEQUENCE CONSSCORE GENENAME R2 PROXY_SNP
			rename GENENAME SYMBOL
			gen str4 X= string(R2,"%05.4f")
			gen Y = CONSS
			destring Y, replace
			drop if Y < 3
			replace CONSEQUENCE = CONSEQUENCE + "(" + CONSS + ")"
			gsort -Y
			egen x = seq(),by(SNP SYMBOL)
			keep if x == 1
			drop Y x
			gen Y = ""
			replace Y = CONSEQUENCE + "(" + PROXY + "; R2 = " + X + ")" if SYMBOL == "NA"
			replace Y = CONSEQUENCE if SNP == PROXY & SYMBOL == "NA"
			replace Y = CONSEQUENCE + ":" + SYMBOL + "(" + PROXY + "; R2 = " + X + ")" if SYMBOL != "NA"
			replace Y = CONSEQUENCE + ":" + SYMBOL if SNP == PROXY & SYMBOL != "NA"
			keep SNP Y
			egen n = seq(),by(SNP)
			reshape wide Y, i(SNP) j(n)
			gen Y999 = ""
			gen CONSEQUENCE = ""
			for var Y1 - Y999: replace CONSEQUENCE = CONSEQUENCE + " " + X if X != ""
			keep SNP CONSEQUENCE		
			save ``function'_summaryqc_short'-`function'-AnnoType-by-snp.dta,replace
			noi checkfile, file(`already_created')
			}
		}
	qui { // map max CADDv1.3 PHRED score
		local already_created ${init_root}/`function'/``function'_summaryqc_short'/``function'_summaryqc_short'-`function'-CADDv1.3-by-snp.dta
		capture confirm file `already_created'		
		if !_rc {
			noi checkfile, file(`already_created')
			}
		else {
			use `cadd'/eur-1000g-phase3-chrall-mac5_CADDv1.3.dta, clear
			renvars, upper
			rename SNP PROXY_SNP
			merge m:m PROXY_SNP using ``function'_summaryqc_short'-ld-pairs.dta
			keep if _m == 3
			keep SNP PHRED R2 PROXY_SNP
			destring PHRED, replace
			gen str4 CADD_PHRED= string(PHRED,"%05.4f")
			gen str4 X= string(R2,"%05.4f")
			gsort -PHRED
			egen x = seq(),by(SNP)
			keep if x == 1
			order SNP CADD
			gen Y = ""
			replace Y = CADD + "(" + PROXY + "; R2 = " + X + ")" 
			replace Y = CADD if SNP == PROXY 
			keep SNP Y 
			rename Y CADD_PHRED
			save ``function'_summaryqc_short'-`function'-CADDv1.3-by-snp.dta,replace
			noi checkfile, file(`already_created')
			}
		}
	qui { // map PolyPhen-SIFT
		local already_created ${init_root}/`function'/``function'_summaryqc_short'/``function'_summaryqc_short'-`function'-PolyPhen-SIFT-by-snp.dta
		capture confirm file `already_created'		
		if !_rc {
			noi checkfile, file(`already_created')
			}
		else {
			use `cadd'/eur-1000g-phase3-chrall-mac5_PolyPhen-SIFT.dta, clear
			renvars, upper
			rename SNP PROXY_SNP
			merge m:m PROXY_SNP using ``function'_summaryqc_short'-ld-pairs.dta
			keep if _m == 3
			keep SNP POLYPHENCAT SIFTCAT R2 PROXY_SNP GENENAME
			gen str4 X= string(R2,"%05.4f")
			for var POLYPHENCAT SIFTCAT: replace X = "" if X == "NA"
			gen Y = ""
			replace Y = POLYPHENCAT + ":" + GENENAME + "(" + PROXY + "; R2 = " + X + ")" if POLYPHENCAT != ""
			replace Y = POLYPHENCAT + ":" + GENENAME  if POLYPHENCAT != "" & SNP == PROXY 
			gen Z = ""
			replace Z = SIFTCAT + ":" + GENENAME + "(" + PROXY + "; R2 = " + X + ")" if SIFTCAT != ""
			replace Z = SIFTCAT + ":" + GENENAME  if SIFTCAT != "" & SNP == PROXY 
			keep SNP Y Z
			rename (Y Z) (POLYPHENCAT SIFTCAT)
			for var POLYPHENCAT SIFTCAT : replace X = "NA" if X == ""
			save tmp.dta, replace
		
			use tmp.dta,replace
			keep SNP POLYPHENCAT
			duplicates drop
			sort POLYPHENCAT
			egen n = seq(),by(SNP)
			reshape wide POLYPHENCAT, i(SNP) j(n)
			gen POLYPHENCAT999 = ""
			gen POLYPHENCAT = ""
			for var POLYPHENCAT1 - POLYPHENCAT999: replace POLYPHENCAT = POLYPHENCAT + " " + X if X != "" & X != "NA"
			keep SNP POLYPHENCAT
			replace PO = "NA" if PO == ""
			save tmp2.dta,replace
			
			use tmp.dta,replace
			keep SNP SIFTCAT
			duplicates drop
			sort SIFTCAT
			egen n = seq(),by(SNP)
			reshape wide SIFTCAT, i(SNP) j(n)
			gen SIFTCAT999 = ""
			gen SIFTCAT = ""
			for var SIFTCAT1 - SIFTCAT999: replace SIFTCAT = SIFTCAT + " " + X if X != "" & X != "NA"
			keep SNP SIFTCAT
			replace SIFTCAT = "NA" if SIFTCAT == ""
			merge 1:1 SNP using tmp2.dta
			keep if _m == 3
			drop if SIFTCAT == "NA" & PO == "NA"
			keep SNP PO SI
			save ``function'_summaryqc_short'-`function'-PolyPhen-SIFT-by-snp.dta, replace
			noi checkfile, file(`already_created')
			}
		}
	qui { // map max ChromHMM
		local already_created ${init_root}/`function'/``function'_summaryqc_short'/``function'_summaryqc_short'-`function'-ChromHMM-by-snp.dta
		capture confirm file `already_created'		
		if !_rc {
			noi checkfile, file(`already_created')
			}
		else {
			use `cadd'/eur-1000g-phase3-chrall-mac5_ChromHMM.dta, clear
			renvars, upper
			rename SNP PROXY_SNP
			merge m:m PROXY_SNP using ``function'_summaryqc_short'-ld-pairs.dta
			keep if _m == 3
			keep SNP CHMMTSSA - CHMMQUIES R2 PROXY_SNP
			for var CHMMTSSA -CHMMQUIES: destring X, replace
			gen str4 X= string(R2,"%05.4f")
			save tmp.dta,replace
			foreach i in CHMMTSSA CHMMTSSAFLNK CHMMTXFLNK CHMMTX CHMMTXWK CHMMENHG CHMMENH CHMMZNFRPTS CHMMHET CHMMTSSBIV CHMMBIVFLNK CHMMENHBIV CHMMREPRPC CHMMREPRPCWK CHMMQUIES {
				use tmp.dta, clear
				collapse (max) `i', by(SNP)
				save tmp_`i'.dta, replace
				merge 1:m SNP `i' using tmp.dta
				keep if _m == 3
				keep SNP PROXY_SNP X `i' R2
				gsort -R2
				egen x = seq(),by(SNP `i')
				keep if x == 1
				gen str4 Y= string(`i',"%05.4f")
				gen max`i' = ""
				replace max`i' = Y + "(" + PROXY + "; R2 = " + X + ")" 
				replace max`i' = Y if SNP == PROXY 
				replace max`i' = Y if Y == "0.00"
				keep SNP max`i'
				save tmp_`i'.dta, replace
				}
			use tmp_CHMMTSSA.dta, clear
			foreach i in  CHMMTSSAFLNK CHMMTXFLNK CHMMTX CHMMTXWK CHMMENHG CHMMENH CHMMZNFRPTS CHMMHET CHMMTSSBIV CHMMBIVFLNK CHMMENHBIV CHMMREPRPC CHMMREPRPCWK CHMMQUIES {
				merge 1:1 SNP using tmp_`i'.dta
				keep if _m == 3
				drop _m
				}
			save ``function'_summaryqc_short'-`function'-ChromHMM-by-snp.dta, replace
			noi checkfile, file(`already_created')
			}
		}
	qui { // map max Cons
		/* Conservation (GERP) - not included
		   The Genomic Evolutionary Rate Profiling (GERP) score was obtained 
		   from the GERP website in September of 2011. It ranges from -12.3 to
		   6.17, with 6.17 being the most conserved. 
		   http://mendel.stanford.edu/SidowLab/downloads/gerp/index.html
		   http://genome.cshlp.org/content/15/7/901.full
		 */
		local already_created ${init_root}/`function'/``function'_summaryqc_short'/``function'_summaryqc_short'-`function'-Cons-by-snp.dta
		capture confirm file `already_created'		
		if !_rc {
			noi checkfile, file(`already_created')
			}
		else {
			use `cadd'/eur-1000g-phase3-chrall-mac5_Cons.dta, clear
			renvars, upper
			rename SNP PROXY_SNP
			merge m:m PROXY_SNP using ``function'_summaryqc_short'-ld-pairs.dta
			keep if _m == 3
			keep SNP PRIPHCONS MAMPHCONS VERPHCONS PRIPHYLOP MAMPHYLOP VERPHYLOP GERPN GERPS GERPRS GERPRSPVAL  FITCONS R2 PROXY_SNP
			for var  PRIPHCONS MAMPHCONS VERPHCONS PRIPHYLOP MAMPHYLOP VERPHYLOP GERPN GERPS GERPRS GERPRSPVAL  FITCONS R2 : destring X, replace force
			replace GERPRSPVAL = 1e-100 if GERPRSPVAL <=1e-100
			gen str4 X= string(R2,"%05.4f")
			save tmp.dta,replace
			foreach i in PRIPHCONS MAMPHCONS VERPHCONS PRIPHYLOP MAMPHYLOP VERPHYLOP FITCONS {
				use tmp.dta, clear
				keep SNP `i'
				sort `i'
				drop if `i' == .
				save tmp_`i'.dta, replace
				collapse (min) `i', by(SNP)
				save tmp_`i'_min.dta, replace
				use tmp_`i'.dta, clear
				collapse (max) `i', by(SNP)
				save tmp_`i'_max.dta, replace
				use tmp_`i'_min.dta, clear
				rename `i' min
				merge 1:1 SNP using tmp_`i'_max.dta
				rename `i' max
				gen str4 X= string(min,"%05.4f")
				gen str4 Y= string(max,"%05.4f")
				gen `i' = X + "-" + Y
				keep SNP `i'
				save tmp_`i'.dta, replace
				}
			use tmp_PRIPHCONS.dta, clear
			foreach i in  MAMPHCONS VERPHCONS PRIPHYLOP MAMPHYLOP VERPHYLOP FITCONS {
				merge 1:1 SNP using tmp_`i'.dta
				keep if _m == 3
				drop _m
				}
			save ``function'_summaryqc_short'-`function'-Cons-by-snp.dta, replace
			noi checkfile, file(`already_created')
			}
		}
	qui { // map ENCODE
		local already_created ${init_root}/`function'/``function'_summaryqc_short'/``function'_summaryqc_short'-`function'-ENCODE-by-snp.dta
		capture confirm file `already_created'		
		if !_rc {
			noi checkfile, file(`already_created')
			}
		else {
			use `cadd'/eur-1000g-phase3-chrall-mac5_ENCODE.dta, clear
			renvars, upper
			rename SNP PROXY_SNP
			merge m:m PROXY_SNP using ``function'_summaryqc_short'-ld-pairs.dta
			keep if _m == 3
			keep SNP ENCOCCOMBPVAL ENCOCDNASEPVAL ENCOCFAIREPVAL ENCOCPOLIIPVAL ENCOCCTCFPVAL ENCOCMYCPVAL R2 PROXY_SNP
			for var ENCOCCOMBPVAL ENCOCDNASEPVAL ENCOCFAIREPVAL ENCOCPOLIIPVAL ENCOCCTCFPVAL ENCOCMYCPVAL: destring X, replace force
			gen str4 X= string(R2,"%05.4f")
			save tmp.dta,replace 
			foreach i in ENCOCCOMBPVAL ENCOCDNASEPVAL ENCOCFAIREPVAL ENCOCPOLIIPVAL ENCOCCTCFPVAL ENCOCMYCPVAL {
				use tmp.dta, clear
				keep PROXY `i' SNP R2 X
				gen str4 `i'_PHRED= string(`i',"%05.4f")
				drop if `i' == .
				gsort -`i'
				egen x = seq(),by(SNP)
				keep if x == 1
				order SNP `i'
				gen Y = ""
				replace Y = `i'_PHRED + "(" + PROXY + "; R2 = " + X + ")" 
				replace Y = `i'_PHRED if SNP == PROXY 
				keep SNP Y 
				rename Y `i'_PHRED
				save tmp_`i'.dta, replace
				}
			use tmp_ENCOCCOMBPVAL.dta, clear
			foreach i in  ENCOCDNASEPVAL ENCOCFAIREPVAL ENCOCPOLIIPVAL ENCOCCTCFPVAL ENCOCMYCPVAL {
				merge 1:1 SNP using tmp_`i'.dta
				keep if _m == 3
				drop _m
				}
			save ``function'_summaryqc_short'-`function'-ENCODE-by-snp.dta, replace
			noi checkfile, file(`already_created')
			}
		}		
	qui { // map mirSVR
		local already_created ${init_root}/`function'/``function'_summaryqc_short'/``function'_summaryqc_short'-`function'-mirSVR-by-snp.dta
		capture confirm file `already_created'		
		if !_rc {
			noi checkfile, file(`already_created')
			}
		else {
			use `cadd'/eur-1000g-phase3-chrall-mac5_mirSVR.dta, clear
			renvars, upper
			rename SNP PROXY_SNP
			merge m:m PROXY_SNP using ``function'_summaryqc_short'-ld-pairs.dta
			keep if _m == 3
			keep SNP MIRSVRSCORE R2 PROXY_SNP
			for var MIRSVRSCORE : destring X, replace force
			gen str4 X= string(R2,"%05.4f")
			keep PROXY MIRSVRSCORE SNP R2 X
			gen str6 Y= string(MIRSVRSCORE,"%06.5f")
			drop if MIRSVRSCORE == .
			gsort -MIRSVRSCORE
			egen x = seq(),by(SNP)
			keep if x == 1
			order SNP MIRSVRSCORE
			drop MIRSVRSCORE
			gen MIRSVRSCORE = ""
			replace MIRSVRSCORE = Y + "(" + PROXY + "; R2 = " + X + ")" 
			replace MIRSVRSCORE = Y if SNP == PROXY 
			keep SNP MIRSVRSCORE 
			save ``function'_summaryqc_short'-`function'-mirSVR-by-snp.dta, replace
			noi checkfile, file(`already_created')
			}
		}
	qui { // map TFBSmotif
		local already_created ${init_root}/`function'/``function'_summaryqc_short'/``function'_summaryqc_short'-`function'-TFBSmotif-by-snp.dta
		capture confirm file `already_created'		
		if !_rc {
			noi checkfile, file(`already_created')
			}
		else {
			use `cadd'/eur-1000g-phase3-chrall-mac5_TFBSmotif.dta, clear
			renvars, upper
			rename SNP PROXY_SNP
			merge m:m PROXY_SNP using ``function'_summaryqc_short'-ld-pairs.dta
			keep if _m == 3
			keep SNP TFBS R2 PROXY_SNP
			for var TFBS: destring X, replace force
			gen str4 X= string(R2,"%05.4f")
			keep PROXY TFBS SNP R2 X
			gsort -TFBS
			tostring TFBS, replace
			drop if TFBS == "."
			egen x = seq(),by(SNP)
			keep if x == 1
			gen Y = ""
			replace Y = TFBS + "(" + PROXY + "; R2 = " + X + ")" 
			replace Y = TFBS if SNP == PROXY 
			keep SNP Y 
			rename Y TFBS
			save ``function'_summaryqc_short'-`function'-TFBSmotif-by-snp.dta, replace
			noi checkfile, file(`already_created')
			}
		}
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 9 - create summary_excel_sheet
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 9: create Supplementary Tables"
	noi di as text"#########################################################################"	
	local already_created ${init_root}/`function'/``function'_summaryqc_short'/``function'_summaryqc_short'-`function'-supplementary-tables.xlsx
	capture confirm file `already_created'		
	if !_rc {
		noi checkfile, file(`already_created')
		}
	else {
		clear
		set obs 100
		gen n = _n
		qui { // define sheets
			gen sheet     = ""
			replace sheet = "gene_protein_coding" 			if n == 1
			replace sheet = "gene_other(1)" 					if n == 2
			replace sheet = "gene_other(2)" 					if n == 3
			replace sheet = "eqtl_compact" 					if n == 4
			replace sheet = "eqtl_tissue(full)"				if n == 5
			replace sheet = "eqtl_super(full)" 				if n == 6
			replace sheet = "consequence"			 			if n == 7
			replace sheet = "epigenetic_marks"				if n == 8
			drop if sheet == ""
			}
		qui { // add description
			gen description = ""
			gen link = ""
			qui { // gene
				replace description = "Annotation of transcripts (protein coding genes) with physical overlap with the associatied region (LD_RANGE). Data Source: (ftp://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz)" if sheet == "gene_protein_coding"
				replace description = "Annotation of transcripts (pseudogene lincRNA antisense miRNA misc_RNA snRNA snoRNA) with physical overlap with the associatied region (LD_RANGE). Data Source: (ftp://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz)" if sheet == "gene_other(1)"
				replace description = "Annotation of transcripts (sense_intronic rRNA processed_transcript sense_overlapping IG_V_pseudogene IG_V_gene TR_V_gene TR_J_gene polymorphic_pseudogene IG_D_gene  Mt_tRNA _3prime_overlapping_ncrna IG_J_gene IG_C_gene IG_C_pseudogene TR_C_gene TR_J_pseudogene TR_D_gene IG_J_pseudogene Mt_rRNA) with physical overlap with the associatied region (LD_RANGE). Data Source: (ftp://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz)" if sheet == "gene_other(2)"
				}
			qui { // eqtl
				replace description = "Annotation of cis-eQTL P<1e-5 to SNP and LD partners (R^2 > 0.8) - summary report gene + super tissue only. Data Source: GTEx_v7" if sheet == "eqtl_compact"
				replace description = "Annotation of cis-eQTL P<1e-5 to SNP and LD partners (R^2 > 0.8) - full report for individual tissue (max eQTL reported per index SNP). Data Source: GTEx_v7" if sheet == "eqtl_tissue(full)"
				replace description = "Annotation of cis-eQTL P<1e-5 to SNP and LD partners (R^2 > 0.8) - full report for super tissue (max eQTL reported per index SNP). Data Source: GTEx_v7" if sheet == "eqtl_super(full)"
				}
			qui { // consequence
				replace description = "Annotation of variant consequence (Consequence, CADD Phred Score, PolyPhen Category and SIFT Category) to SNP and LD partners (R^2 > 0.8) - (Consequence report limited to CONSSCORE >2; 5PRIME_UTR(3) REGULATORY(4) NONCODING_CHANGE,SPLIC_SITE,SYNONYMOUS,UNKNOWN(5) CANONICAL_SPLICE(6) FRAME_SHIFT,INFRAME,NON_SYNONYMOUS,STOP_LOST(7) STOP_GAINED(8)); report limited to max. CONSSCORE per gene per index SNP. Data Source: cadd.gs.washington.edu" if sheet == "consequence"
				}
			qui { // epigenetic_marks	
				replace description = "Annotation of ChromHMM and ENCODE epigenetic marks to SNP and LD partners (R^2 > 0.8). Data Source: cadd.gs.washington.edu" if sheet == "epigenetic_marks"
				}
		
			replace link  = `"=HYPERLINK(""' + sheet + `"!A1",""' + sheet + `"")"' 
			export excel sheet description using ``function'_summaryqc_short'-`function'-supplementary-tables.xlsx, firstrow(variables) sheet("index") replace
			}
		qui { // gene
			qui { // notes 
				/*
				* class = 1
				protein_coding	20327
				* class = 2
				pseudogene	13920
				lincRNA	7109
				antisense	5273
				miRNA	3049
				misc_RNA	2033
				snRNA	1916
				snoRNA	1457
				* class = 3
				sense_intronic	741
				rRNA	526
				processed_transcript	514
				sense_overlapping	202
				IG_V_pseudogene	187
				IG_V_gene	138
				TR_V_gene	97
				TR_J_gene	74
				polymorphic_pseudogene	45
				IG_D_gene	37
				TR_V_pseudogene	27
				Mt_tRNA	22
				3prime_overlapping_ncrna	21
				IG_J_gene	18
				IG_C_gene	14
				IG_C_pseudogene	9
				TR_C_gene	5
				TR_J_pseudogene	4
				TR_D_gene	3
				IG_J_pseudogene	3
				Mt_rRNA	2
				*/
				}
			qui { // gene_protein_coding
				use ${init_root}/summaryqc2top/``function'_summaryqc_short'/``function'_summaryqc_short'-summaryqc2top.dta, clear
				merge 1:1 SNP using ``function'_summaryqc_short'-`function'-gene-by-range.dta
				drop if _me == 2
				sort P
				keep  SNP POSITION A1 A2 OR P LD_RANGE protein_coding
				export excel using ``function'_summaryqc_short'-`function'-supplementary-tables.xlsx, firstrow(variables) sheet("gene_protein_coding") sheetreplace
			}
			qui { // gene_other(1)
				use ${init_root}/summaryqc2top/``function'_summaryqc_short'/``function'_summaryqc_short'-summaryqc2top.dta, clear
				merge 1:1 SNP using ``function'_summaryqc_short'-`function'-gene-by-range.dta
				drop if _me == 2
				sort P
				keep  SNP POSITION A1 A2 OR P LD_RANGE pseudogene	lincRNA antisense miRNA misc_RNA snRNA snoRNA
				export excel using ``function'_summaryqc_short'-`function'-supplementary-tables.xlsx, firstrow(variables) sheet("gene_other(1)") sheetreplace
				}
			qui { // gene_other(2)
				use ${init_root}/summaryqc2top/``function'_summaryqc_short'/``function'_summaryqc_short'-summaryqc2top.dta, clear
				merge 1:1 SNP using ``function'_summaryqc_short'-`function'-gene-by-range.dta
				drop if _me == 2
				sort P
				keep  SNP POSITION A1 A2 OR P LD_RANGE sense_intronic rRNA processed_transcript sense_overlapping IG_V_pseudogene IG_V_gene TR_V_gene TR_J_gene polymorphic_pseudogene IG_D_gene  Mt_tRNA _3prime_overlapping_ncrna IG_J_gene IG_C_gene IG_C_pseudogene TR_C_gene TR_J_pseudogene TR_D_gene IG_J_pseudogene Mt_rRNA
				export excel using ``function'_summaryqc_short'-`function'-supplementary-tables.xlsx, firstrow(variables) sheet("gene_other(2)") sheetreplace
			}		
			}
		qui { // eqtl
			qui { // eqtl_compact
				use ${init_root}/summaryqc2top/``function'_summaryqc_short'/``function'_summaryqc_short'-summaryqc2top.dta, clear
				merge 1:m SNP using ``function'_summaryqc_short'-`function'-eqtl-by-snp-compact.dta
				drop if _me == 2
				sort P		
				export excel SNP POS A1 A2 OR P EQTL using ``function'_summaryqc_short'-`function'-supplementary-tables.xlsx, firstrow(variables) sheet("eqtl_compact") sheetreplace
				}
			qui { // eqtl_tissue(full)
				use ${init_root}/summaryqc2top/``function'_summaryqc_short'/``function'_summaryqc_short'-summaryqc2top.dta, clear
				merge 1:m SNP using ``function'_summaryqc_short'-`function'-eqtl-by-snp-verbose.dta
				drop if _me == 2
				sort P		
				drop if TISSUE_EQTL == ""
				export excel SNP POS A1 A2 OR P TISSUE_EQTL using ``function'_summaryqc_short'-`function'-supplementary-tables.xlsx, firstrow(variables) sheet("eqtl_tissue(full)") sheetreplace
				}
			qui { // eqtl_super(full)
				use ${init_root}/summaryqc2top/``function'_summaryqc_short'/``function'_summaryqc_short'-summaryqc2top.dta, clear
				merge 1:m SNP using ``function'_summaryqc_short'-`function'-eqtl-by-snp-verbose.dta
				drop if _me == 2
				sort P		
				drop if SUPER == ""
				export excel SNP POS A1 A2 OR P SUPER using ``function'_summaryqc_short'-`function'-supplementary-tables.xlsx, firstrow(variables) sheet("eqtl_super(full)") sheetreplace
				}
			}
		qui { // consequence
			use ${init_root}/summaryqc2top/``function'_summaryqc_short'/``function'_summaryqc_short'-summaryqc2top.dta, clear
			merge 1:1 SNP using ``function'_summaryqc_short'-`function'-AnnoType-by-snp.dta
			drop if _me == 2
			keep SNP POS A1 A2 OR P CONSEQ
			sort P
			merge 1:1 SNP using ``function'_summaryqc_short'-`function'-CADDv1.3-by-snp.dta
			drop if _me == 2
			keep SNP POS A1 A2 OR P CONSEQ CADD
			sort P		
			merge 1:1 SNP using ``function'_summaryqc_short'-`function'-PolyPhen-SIFT-by-snp.dta
			drop if _me == 2
			sort P				
			export excel SNP POS A1 A2 OR P CONSEQ CADD POLYPHEN SIFT using ``function'_summaryqc_short'-`function'-supplementary-tables.xlsx, firstrow(variables) sheet("consequence") sheetreplace
			}
		qui { // epigenetic_marks
			use ${init_root}/summaryqc2top/``function'_summaryqc_short'/``function'_summaryqc_short'-summaryqc2top.dta, clear
			merge 1:1 SNP using ``function'_summaryqc_short'-`function'-ChromHMM-by-snp.dta
			drop if _me == 2
			sort P		
			drop _m
			merge 1:1 SNP using ``function'_summaryqc_short'-`function'-ENCODE-by-snp.dta
			drop if _me == 2
			sort P		
			drop _m	
			export excel SNP POS A1 A2 OR P maxCHMMTSSA - ENCOCMYCPVAL_PHRED  using ``function'_summaryqc_short'-`function'-supplementary-tables.xlsx, firstrow(variables) sheet("epigenetic_marks") sheetreplace
			}
		noi checkfile, file(`already_created')
		}	
	noi di as text"#########################################################################"
	noi di as text""	
	}		
		
qui { // module 9 - write log file
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
	file write myfile `"# > ...... saved data ``function'_summaryqc_short'-`function'-gene-by-snp.dta"' _n
	file write myfile `"# > ...... saved data ``function'_summaryqc_short'-`function'-gene-by-range.dta"' _n
	file write myfile `"# > ...... saved data ``function'_summaryqc_short'-`function'-eqtl-by-snp-verbose.dta"' _n
	file write myfile `"# > ...... saved data ``function'_summaryqc_short'-`function'-eqtl-by-snp-compact.dta"' _n
	file write myfile `"# > ...... saved data ``function'_summaryqc_short'-`function'-AnnoType-by-snp.dta"' _n
	file write myfile `"# > ...... saved data ``function'_summaryqc_short'-`function'-CADDv1.3-by-snp.dta"' _n
	file write myfile `"# > ...... saved data ``function'_summaryqc_short'-`function'-ChromHMM-by-snp.dta"' _n
	file write myfile `"# > ...... saved data ``function'_summaryqc_short'-`function'-Cons-by-snp.dta"' _n
	file write myfile `"# > ...... saved data ``function'_summaryqc_short'-`function'-ENCODE-by-snp.dta"' _n
	file write myfile `"# > ...... saved data ``function'_summaryqc_short'-`function'-mirSVR-by-snp.dta"' _n
	file write myfile `"# > ...... saved data ``function'_summaryqc_short'-`function'-PolyPhen-SIFT-by-snp.dta"' _n
	file write myfile `"# > ...... saved data ``function'_summaryqc_short'-`function'-TFBSmotif-by-snp.dta"' _n
	file write myfile `"# > ...... saved data ``function'_summaryqc_short'-`function'-supplementary-tables.xlsx"' _n
	
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > overview of `function'"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > 1. identify top N SNPs"' _n
	file write myfile `"# > 2. extend LD to capture proxy_SNPs"' _n
	file write myfile `"# > 3. map annotations"' _n
	file write myfile `"#########################################################################"' _n
	file close myfile	
	}
qui { // module 5 - clean folder
	files2dta, dir(`c(pwd)')
	gen keep = .
	replace keep = 1 if file == "``function'_summaryqc_short'-`function'.log"
	replace keep = 1 if file == "``function'_summaryqc_short'-`function'-gene-by-snp.dta"
	replace keep = 1 if file == "``function'_summaryqc_short'-`function'-gene-by-range.dta"
	replace keep = 1 if file == "``function'_summaryqc_short'-`function'-eqtl-by-snp-verbose.dta"
	replace keep = 1 if file == "``function'_summaryqc_short'-`function'-eqtl-by-snp-compact.dta"
	replace keep = 1 if file == "``function'_summaryqc_short'-`function'-AnnoType-by-snp.dta"
	replace keep = 1 if file == "``function'_summaryqc_short'-`function'-CADDv1.3-by-snp.dta"
	replace keep = 1 if file == "``function'_summaryqc_short'-`function'-ChromHMM-by-snp.dta"
	replace keep = 1 if file == "``function'_summaryqc_short'-`function'-Cons-by-snp.dta"
	replace keep = 1 if file == "``function'_summaryqc_short'-`function'-ENCODE-by-snp.dta"
	replace keep = 1 if file == "``function'_summaryqc_short'-`function'-mirSVR-by-snp.dta"
	replace keep = 1 if file == "``function'_summaryqc_short'-`function'-PolyPhen-SIFT-by-snp.dta"
	replace keep = 1 if file == "``function'_summaryqc_short'-`function'-TFBSmotif-by-snp.dta"
	replace keep = 1 if file == "``function'_summaryqc_short'-`function'-supplementary-tables.xlsx"
	
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

   
