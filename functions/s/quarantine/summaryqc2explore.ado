*! 1.0.1 Richard Anney 29jan2020

program  summaryqc2explore
local    function summaryqc2explore
local    `function'_version v1

syntax , summaryqc(string asis) ldref(string asis) generef(string asis) [snp(string asis) snp_list(string asis) gene(string asis) gene_list(string asis) range(string asis) range_name(string asis) range_list(string asis) recombref(string asis)]

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
	ref_path_to_short, ref(`summaryqc')
	local `function'_sqc_short ${ref_short}
	ref_path_to_short, ref(`ldref')
	local `function'_ldref_short ${ref_short}
	cd ${init_root}	
	!mkdir `function'
	cd `function'
	noi di as text"# > ..... project name "as result"``function'_sqc_short'"
	!mkdir ``function'_sqc_short'
	cd ``function'_sqc_short'
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
	qui { // cheack ldref
		foreach file in bim bed fam {
			noi checkfile, file(`ldref'.`file')
			}
		}	
	qui { // check gene 
		foreach i in _3prime_overlapping_ncrna IG_C_gene IG_C_pseudogene IG_D_gene IG_J_gene IG_J_pseudogene IG_V_gene IG_V_pseudogene Mt_rRNA Mt_tRNA TR_C_gene TR_D_gene TR_J_gene TR_J_pseudogene TR_V_gene antisense lincRNA miRNA misc_RNA polymorphic_pseudogene processed_transcript protein_coding pseudogene rRNA sense_intronic sense_overlapping snRNA snoRNA {
			noi checkfile, file(`generef'/Homo_sapiens.GRCh37.87.gtf_exon-`i'.hg-list)
			noi checkfile, file(`generef'/Homo_sapiens.GRCh37.87.gtf-`i'.hg-list)
			noi checkfile, file(`generef'/Homo_sapiens.GRCh37.87.gtf-`i'.hg-list-35+5.bed)
			}
		}
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 3 - explore associations
	qui { // by snp
		qui { // define index and LD partners
			clear 
			set obs 1
			gen snp = "`snp'"
			if snp == "" {
				}
			else  {
				noi di as text"#########################################################################"
				noi di as text"# SECTION - 3: defining association by "as result "`snp'" as text " and LD partners"
				noi di as text"#########################################################################"
				qui { // define ld-pairs
					noi di as text"# > define LD pairs of "as result"`snp'"
					local already_created ${init_root}/snp2ldpair/``function'_ldref_short'/``function'_ldref_short'-`snp'-snp2ldpair.dta
					capture confirm file `already_created'
					if !_rc {
						noi checkfile, file(`already_created')
						}
					else {
						snp2ldpair, snp(`snp') ldref(`ldref')
						noi checkfile, file(`already_created')
						}
					}
				qui { // process test data-sets
					noi di as text"# > ........processing "as result"``function'_sqc_short'"
					local already_created ``function'_sqc_short'-index.dta
					capture confirm file `already_created'
					if !_rc {
						noi checkfile, file(`already_created')
						}
					else {
						use snp a1 a2 beta se p using `summaryqc'-summaryqc.dta, clear
						gen GWAS = "``function'_sqc_short'"
						rename (snp a1 a2) (PROXY_SNP INDEX_A1 INDEX_A2)
						renvars, upper
						save ``function'_sqc_short'-index.dta, replace
						noi checkfile, file(`already_created')
						}
					local already_created ``function'_sqc_short'-ref.dta
					capture confirm file `already_created'
					if !_rc {
						noi checkfile, file(`already_created')
						}
					else {
						use snp a1 a2 using `summaryqc'-summaryqc.dta, clear
						recodegenotype 
						rename (a1 recode) (REF INDEX_GT)
						renvars, upper
						keep SNP REF INDEX_GT
						save ``function'_sqc_short'-ref.dta, replace
						noi checkfile, file(`already_created')
						}
					}
				qui { // merge ld and test data
					noi di as text"# > merge LD data with "as result"``function'_sqc_short'"
					local already_created ``function'_sqc_short'-`snp'-`function'.dta
					capture confirm file `already_created'
					if !_rc {
						noi checkfile, file(`already_created')
						}
					else {
						use ``function'_sqc_short'-index.dta, clear
						merge 1:1 PROXY_SNP using ${init_root}/snp2ldpair/``function'_ldref_short'/``function'_ldref_short'-`snp'-snp2ldpair.dta
						keep if _m == 3
						gen KEEP = 1 if PROXY_SNP == SNP
						sort KEEP SNP R2
						egen x = seq(),by(SNP)
						keep if x == 1
						drop _m KEEP x
						merge m:1 SNP using ``function'_sqc_short'-ref.dta
						keep if _m == 3
						sort SNP GWAS
						gen TAG = ""
						replace TAG = "SNP"   if SNP == PROXY_SNP
						replace TAG = "PROXY" if SNP != PROXY_SNP
						replace TAG = "INDEX" if GWAS == "``function'_sqc_short'"
						gen str4 X= string(R2,"%05.4f")
						gen SNP_ = ""
						replace SNP_ = SNP + "-" + REF if TAG == "INDEX"
						replace SNP_ = SNP + "-" + SNP_A1 if TAG == "SNP"
						replace SNP_ = PROXY_SNP  + "-" + PROXY_A1 + "(" + X + ")" if TAG == "PROXY" & REF == SNP_A1 
						replace SNP_ = PROXY_SNP  + "-" + PROXY_A2 + "(" + X + ")" if TAG == "PROXY" & REF == SNP_A2 
						gen BETA_ = .
						replace BETA_ =  BETA if TAG == "INDEX"
						replace BETA_ =  BETA if TAG == "SNP" & REF == SNP_A1
						replace BETA_ = -BETA if TAG == "SNP" & REF == SNP_A2
						replace BETA_ =  BETA if TAG == "PROXY" & REF == SNP_A1
						replace BETA_ = -BETA if TAG == "PROXY" & REF == SNP_A2
						keep SNP GWAS SNP_ BETA_ SE P TAG
						sort SNP GWAS
						renvars, upper
						gen OR   = exp(BETA)
						gen LB  = OR - (1.96 * SE)
						gen UB  = OR + (1.96 * SE)
						gen str4 X= string(OR,"%05.4f")
						gen str4 Y= string(LB,"%05.4f")
						gen str4 Z_= string(UB,"%05.4f") 
						drop OR
						gen OR = X + " (" + Y + "-" + Z_ + ")"
						keep  SNP SNP_ TAG GWAS OR P
						order SNP SNP_ TAG GWAS OR P
						lab var SNP   "Index SNP"
						lab var SNP_  "Reported SNP"
						lab var TAG   "Relationship to INDEX"
						lab var GWAS  "GWAS"
						lab var OR    "Odds Ratio (95%CI)"
						lab var P     "P-value"
						save ``function'_sqc_short'-`snp'-`function'.dta, replace
						noi checkfile, file(`already_created')
						}
					}
				qui { // report to excel
					noi di as text"# > ........ report to "as result"``function'_sqc_short'-`function'.xlsx"
					use `already_created', clear
					export excel ``function'_sqc_short'-`function'.xlsx, firstrow(varlabels) sheet("`snp'") sheetreplace
					}
				noi di as text"#########################################################################"
				noi di as text""
				}
			}
		}
	qui { // by snp-list
		qui { // define index and LD partners
			clear 
			set obs 1
			gen snp_list = "`snp_list'"
			if snp_list == "" {
				}
			else  {
				noi di as text"#########################################################################"
				noi di as text"# SECTION - 3: defining association by "as result "`snp-list'" as text " (snp-list) and LD partners"
				noi di as text"#########################################################################"
				ref_path_to_short, ref(`snp_list')
				local `function'_sl_short ${ref_short}
				qui { // define ld-pairs
					noi di as text"# > define LD pairs of "as result"``function'_sl_short'"
					local already_created ${init_root}/snp2ldpair/``function'_ldref_short'/``function'_ldref_short'-``function'_sl_short'-snp2ldpair.dta
					capture confirm file `already_created'
					if !_rc {
						noi checkfile, file(`already_created')
						}
					else {
						noi snp2ldpair, snp_list(`snp_list') ldref(`ldref')
						noi checkfile, file(`already_created')
						}
					}
				qui { // process test data-sets
					noi di as text"# > ........processing "as result"``function'_sqc_short'"
					local already_created ``function'_sqc_short'-index.dta
					capture confirm file `already_created'
					if !_rc {
						noi checkfile, file(`already_created')
						}
					else {
						use snp a1 a2 beta se p using `summaryqc'-summaryqc.dta, clear
						gen GWAS = "``function'_sqc_short'"
						rename (snp a1 a2) (PROXY_SNP INDEX_A1 INDEX_A2)
						renvars, upper
						save ``function'_sqc_short'-index.dta, replace
						noi checkfile, file(`already_created')
						}
					local already_created ``function'_sqc_short'-ref.dta
					capture confirm file `already_created'
					if !_rc {
						noi checkfile, file(`already_created')
						}
					else {
						use snp a1 a2 using `summaryqc'-summaryqc.dta, clear
						recodegenotype 
						rename (a1 recode) (REF INDEX_GT)
						renvars, upper
						keep SNP REF INDEX_GT
						save ``function'_sqc_short'-ref.dta, replace
						noi checkfile, file(`already_created')
						}
					}
				qui { // merge ld and test data
					noi di as text"# > merge LD data with "as result"``function'_sqc_short'"
					local already_created ``function'_sqc_short'-``function'_sl_short'-`function'.dta
					capture confirm file `already_created'
					if !_rc {
						noi checkfile, file(`already_created')
						}
					else {
						use ``function'_sqc_short'-index.dta, clear
						merge 1:1 PROXY_SNP using ${init_root}/snp2ldpair/``function'_ldref_short'/``function'_ldref_short'-``function'_sl_short'-snp2ldpair.dta
						keep if _m == 3
						gen KEEP = 1 if PROXY_SNP == SNP
						sort KEEP SNP R2
						egen x = seq(),by(SNP)
						keep if x == 1
						drop _m KEEP x
						merge m:1 SNP using ``function'_sqc_short'-ref.dta
						keep if _m == 3
						sort SNP GWAS
						gen TAG = ""
						replace TAG = "SNP"   if SNP == PROXY_SNP
						replace TAG = "PROXY" if SNP != PROXY_SNP
						replace TAG = "INDEX" if GWAS == "``function'_sqc_short'"
						gen str4 X= string(R2,"%05.4f")
						gen SNP_ = ""
						replace SNP_ = SNP + "-" + REF if TAG == "INDEX"
						replace SNP_ = SNP + "-" + SNP_A1 if TAG == "SNP"
						replace SNP_ = PROXY_SNP  + "-" + PROXY_A1 + "(" + X + ")" if TAG == "PROXY" & REF == SNP_A1 
						replace SNP_ = PROXY_SNP  + "-" + PROXY_A2 + "(" + X + ")" if TAG == "PROXY" & REF == SNP_A2 
						gen BETA_ = .
						replace BETA_ =  BETA if TAG == "INDEX"
						replace BETA_ =  BETA if TAG == "SNP" & REF == SNP_A1
						replace BETA_ = -BETA if TAG == "SNP" & REF == SNP_A2
						replace BETA_ =  BETA if TAG == "PROXY" & REF == SNP_A1
						replace BETA_ = -BETA if TAG == "PROXY" & REF == SNP_A2
						keep SNP GWAS SNP_ BETA_ SE P TAG
						sort SNP GWAS
						renvars, upper
						gen OR   = exp(BETA)
						gen LB  = OR - (1.96 * SE)
						gen UB  = OR + (1.96 * SE)
						gen str4 X= string(OR,"%05.4f")
						gen str4 Y= string(LB,"%05.4f")
						gen str4 Z_= string(UB,"%05.4f") 
						drop OR
						gen OR = X + " (" + Y + "-" + Z_ + ")"
						keep  SNP SNP_ TAG GWAS OR P
						order SNP SNP_ TAG GWAS OR P
						lab var SNP   "Index SNP"
						lab var SNP_  "Reported SNP"
						lab var TAG   "Relationship to INDEX"
						lab var GWAS  "GWAS"
						lab var OR    "Odds Ratio (95%CI)"
						lab var P     "P-value"
						save ``function'_sqc_short'-``function'_sl_short'-`function'.dta, replace
						noi checkfile, file(`already_created')
						}
					}
				qui { // report to excel
					noi di as text"# > ........ report to "as result"``function'_sqc_short'-``function'_sl_short'-`function'.xlsx"
					use `already_created', clear
					export excel ``function'_sqc_short'-`function'.xlsx, firstrow(varlabels) sheet("``function'_sl_short'") sheetreplace
					}
				noi di as text"#########################################################################"
				noi di as text""
				}
			}
		}
	qui { // by gene 
		clear 
		set obs 1
		gen gene = "`gene'"
		if gene == "" {
				}
		else  {
			qui { // plot association in gene region
				noi di as text"#########################################################################"
				noi di as text"# SECTION - 3a: plotting association around "as result "`gene'"
				noi di as text"#########################################################################"
				local already_created ``function'_sqc_short'-`gene'_graphlocuszoom.png
				capture confirm file `already_created'
				if !_rc {
					noi checkfile, file(`already_created')
					}
				else {
					qui { // define range
						noi di as text"# > define location of "as result"`gene'"
						foreach i in _3prime_overlapping_ncrna IG_C_gene IG_C_pseudogene IG_D_gene IG_J_gene IG_J_pseudogene IG_V_gene IG_V_pseudogene Mt_rRNA Mt_tRNA TR_C_gene TR_D_gene TR_J_gene TR_J_pseudogene TR_V_gene antisense lincRNA miRNA misc_RNA polymorphic_pseudogene processed_transcript protein_coding pseudogene rRNA sense_intronic sense_overlapping snRNA snoRNA {
							import delim using `generef'/Homo_sapiens.GRCh37.87.gtf-`i'.hg-list-35+5.bed, clear
							keep if v4 == "`gene'"
							count
							if `r(N)' != 0 {
								local chr 	= v1[1]
								local from 	= v2[1]
								local to 		= v3[1]
								}
							}
						noi di as text"# > .............. chr "as result"`chr'"
						noi di as text"# > .......... from bp "as result"`from'"
						noi di as text"# > ............ to bp "as result"`to'"
						clear
						set obs 1
						gen start = `from'
						gen end = `to'
						gen size = end - start
						gen _sta = .
						gen _end = .
						replace _sta = round(start - 50000,1000)
						replace _end = round(end + 50000,1000)
						local from 	= _sta[1]
						local to 		= _end[1]
						}
					qui { //plot association around gene
						noi di as text"# > .......... running "as result"`gene'"as text " through graphlocuszoom"
						graphlocuszoom, generef(`generef') ldref(`ldref') recombref(`recombref') summaryqc(`summaryqc') range(`chr':`from'..`to') label(on) biotype(all) 
						graph use ${init_root}/graphlocuszoom/``function'_sqc_short'/``function'_sqc_short'-`chr'_`from'_`to'_graphlocuszoom.gph
						graph save ``function'_sqc_short'-`gene'_graphlocuszoom.gph, replace
						graph export ``function'_sqc_short'-`gene'_graphlocuszoom.png, as(png) height(2000) width(6000) replace
						noi checkfile, file(`already_created')
						}
					}
				}
			qui { // calculate gene-based association (magma)
				noi di as text"#########################################################################"
				noi di as text"# SECTION - 3b: defining gene-based association for "as result "`gene'"
				noi di as text"#########################################################################"
				local already_created ``function'_sqc_short'-`gene'-`function'.dta
				capture confirm file `already_created'
				if !_rc {
					noi checkfile, file(`already_created')
					}
				else {
					noi di as text"#########################################################################"
					noi di as text"# SECTION - 3c: run magma by "as result "`gene'"
					noi di as text"#########################################################################"
					qui { // define .snploc .pvaloc n
						noi di as text"# > ......... locating files to process"
						noi di as text"# >  parse snploc from "as result"``function'_sqc_short'"
						noi di as text"# >  parse pvaloc from "as result"``function'_sqc_short'"
						use snp chr bp p n using `summaryqc'-summaryqc.dta, clear
						outsheet snp chr bp using ``function'_sqc_short'.snploc, non noq replace
						outsheet snp p n using    ``function'_sqc_short'.pvaloc, noq replace
						}
					qui { // define geneloc
						noi di as text"# >  parse genloc from "as result"Homo_sapiens.GRCh37.87.gtf"
						foreach i in _3prime_overlapping_ncrna IG_C_gene IG_C_pseudogene IG_D_gene IG_J_gene IG_J_pseudogene IG_V_gene IG_V_pseudogene Mt_rRNA Mt_tRNA TR_C_gene TR_D_gene TR_J_gene TR_J_pseudogene TR_V_gene antisense lincRNA miRNA misc_RNA polymorphic_pseudogene processed_transcript protein_coding pseudogene rRNA sense_intronic sense_overlapping snRNA snoRNA {
							import delim using `generef'/Homo_sapiens.GRCh37.87.gtf_exon-`i'.hg-list, clear
							keep if v4 == "`gene'"
							count
							if `r(N)' > 1 {
								rename (v1 v4 v5 v6) (chr geneid start stop)
								split v7,p(" ")
								gen strand = ""
								replace strand = "+" if v72 == ">"
								replace strand = "-" if v71 == "<"
								order geneid chr start stop strand
								gen genename = geneid
								keep geneid chr start stop strand genename
								duplicates drop
								count
								if `r(N)' != 0 {
									egen x = seq(), by(geneid)
									tostring x,replace
									replace geneid = geneid + "." + x
									drop x
									}
								outsheet geneid chr start stop strand genename using ``function'_sqc_short'-`gene'.geneloc, non noq replace
								}
							}
						}
					qui { // define gene-based p
						noi di as text"# > map snps to genes (+/- 35,10)"
						!magma --annotate window=35,10 --snp-loc ``function'_sqc_short'.snploc --gene-loc ``function'_sqc_short'-`gene'.geneloc --out ``function'_sqc_short'-`gene'
						noi di as text"# > define gene-based pval using magma"
						noi !magma 	--bfile 		`ldref' ///
									--pval 		``function'_sqc_short'.pvaloc use=snp,p ncol=n ///
									--gene-model multi ///
									--gene-annot 	``function'_sqc_short'-`gene'.genes.annot ///
									--genes-only ///
									--gene-settings adap-permp=1000000,10 ///
									--seed 987654321 ///
									--out 		``function'_sqc_short'-`gene'
						!${tabbed} ``function'_sqc_short'-`gene'.genes.out
						import delim using ``function'_sqc_short'-`gene'.genes.out.tabbed, clear varnames(1)
						lab var gene 			"Gene Name"
						lab var chr  			"Chromosome"
						lab var start 			"Transcript Start (hg19) -35kb"
						lab var stop 			"Transcript End (hg19) +10kb"
						lab var nsnps  		"SNPs in model"
						lab var p_snpwise_mean 	"MAGMA P-val (mean)"
						lab var p_snpwise_top1 	"MAGMA P-val (top)"
						lab var zstat            "Z-Score (gene)"
						lab var permp			"MAGMA P-val (perm)"
						lab var nperm            "N Permutations"
						keep  gene chr start stop nsnps p_snpwise_mean p_snpwise_top1 zstat perm nperm
						order gene chr start stop nsnps p_snpwise_mean p_snpwise_top1 zstat perm nperm
						save ``function'_sqc_short'-`gene'-`function'.dta, replace
						export excel ``function'_sqc_short'-`function'.xlsx, firstrow(varlabels) sheet("`gene'") sheetreplace	
						noi checkfile, file(`already_created')
						}
					}
				noi di as text"#########################################################################"
				noi di as text""
				}
			}
		}
	qui { // by gene-list
		clear 
		set obs 1
		gen gene_list = "`gene_list'"
		if gene_list == "" {
				}
		else  {
			qui { // plot association in gene region
				import delim using `gene_list', clear varname(nonames)
				sxpose, clear
				gen _var9999 = ""
				gen gene_list = ""
				for var _var1 - _var9999: replace gene_list = gene_list + " " + X
				local gene_list2 = gene_list[1]
				foreach gene in `gene_list2' {			
					noi di as text"#########################################################################"
					noi di as text"# SECTION - 3a: plotting association around "as result "`gene'"
					noi di as text"#########################################################################"
					local already_created ``function'_sqc_short'-`gene'_graphlocuszoom.png
					capture confirm file `already_created'
					if !_rc {
						noi checkfile, file(`already_created')
						}
					else {
						qui { // define range
							noi di as text"# > define location of "as result"`gene'"
							foreach i in _3prime_overlapping_ncrna IG_C_gene IG_C_pseudogene IG_D_gene IG_J_gene IG_J_pseudogene IG_V_gene IG_V_pseudogene Mt_rRNA Mt_tRNA TR_C_gene TR_D_gene TR_J_gene TR_J_pseudogene TR_V_gene antisense lincRNA miRNA misc_RNA polymorphic_pseudogene processed_transcript protein_coding pseudogene rRNA sense_intronic sense_overlapping snRNA snoRNA {
								import delim using `generef'/Homo_sapiens.GRCh37.87.gtf-`i'.hg-list-35+5.bed, clear
								keep if v4 == "`gene'"
								count
								if `r(N)' > 1 {
									local chr 	= v1[1]
									local from 	= v2[1]
									local to 		= v3[1]
									}
								}
							clear
							set obs 1
							gen from = "`from'"
							if from == "" {
								noi di as text"# > .... could not map "as error"`gene'"
								}
							else {
								noi di as text"# > .............. chr "as result"`chr'"
								noi di as text"# > .......... from bp "as result"`from'"
								noi di as text"# > ............ to bp "as result"`to'"
								clear
								set obs 1
								gen start = `from'
								gen end = `to'
								gen size = end - start
								gen _sta = .
								gen _end = .
								replace _sta = round(start - 50000,1000)
								replace _end = round(end + 50000,1000)
								local from 	= _sta[1]
								local to 		= _end[1]
								gen chr = "`chr'"
								foreach num of num 1/22 {
									if chr == "chr`num'" {
										noi di as text"# > .......... running "as result"`gene'"as text " through graphlocuszoom"
										noi di as text"# > .............. chr "as result"`chr'"
										noi di as text"# > .......... from bp "as result"`from'"
										noi di as text"# > ............ to bp "as result"`to'"	
										graphlocuszoom, generef(`generef') ldref(`ldref') recombref(`recombref') summaryqc(`summaryqc') range(`chr':`from'..`to') label(on) biotype(all)
										graph use ${init_root}/graphlocuszoom/``function'_sqc_short'/``function'_sqc_short'-`chr'_`from'_`to'_graphlocuszoom.gph
										graph save ``function'_sqc_short'-`gene'_graphlocuszoom.gph, replace
										graph export ``function'_sqc_short'-`gene'_graphlocuszoom.png, as(png) height(2000) width(6000) replace
										noi checkfile, file(`already_created')
										}
									}
								if chr == "chr23" {
									noi di as text"# > ...... not running "as error"`gene'"as text " through graphlocuszoom"
									noi di as text"# > .............. chr "as error"`chr'"
									noi di as text"# > .......... from bp "as result"`from'"
									noi di as text"# > ............ to bp "as result"`to'"								
									}
								}
							}
						}
					}
				noi di as text"#########################################################################"
				noi di as text""
				}
			qui { // calculate gene-based association (magma)
				noi di as text"#########################################################################"
				noi di as text"# SECTION - 3b: defining gene-based association for "as result "`gene_list'"
				noi di as text"#########################################################################"
				ref_path_to_short, ref(`gene_list')
				local `function'_gl_short ${ref_short}		
				local already_created ``function'_sqc_short'-``function'_gl_short'-`function'.dta
				capture confirm file `already_created'
				if !_rc {
					noi checkfile, file(`already_created')
					}
				else {
					noi di as text"#########################################################################"
					noi di as text"# SECTION - 3c: run magma by "as result "``function'_gl_short'"
					noi di as text"#########################################################################"
					qui { // define .snploc .pvaloc n
						noi di as text"# > ......... locating files to process"
						noi di as text"# >  parse snploc from "as result"``function'_sqc_short'"
						noi di as text"# >  parse pvaloc from "as result"``function'_sqc_short'"
						use snp chr bp p n using `summaryqc'-summaryqc.dta, clear
						outsheet snp chr bp using ``function'_sqc_short'.snploc, non noq replace
						outsheet snp p n using    ``function'_sqc_short'.pvaloc, noq replace
						}
					qui { // define geneloc
						import delim using `gene_list', clear varname(nonames)
						rename v1 gene
						duplicates drop
						save ``function'_sqc_short'-``function'_gl_short'.dta, replace				
						noi di as text"# >  parse genloc from "as result"Homo_sapiens.GRCh37.87.gtf"
						foreach i in _3prime_overlapping_ncrna IG_C_gene IG_C_pseudogene IG_D_gene IG_J_gene IG_J_pseudogene IG_V_gene IG_V_pseudogene Mt_rRNA Mt_tRNA TR_C_gene TR_D_gene TR_J_gene TR_J_pseudogene TR_V_gene antisense lincRNA miRNA misc_RNA polymorphic_pseudogene processed_transcript protein_coding pseudogene rRNA sense_intronic sense_overlapping snRNA snoRNA {
							import delim using `generef'/Homo_sapiens.GRCh37.87.gtf_exon-`i'.hg-list, clear
							gen gene = v4
							merge m:1 gene using ``function'_sqc_short'-``function'_gl_short'.dta
							keep if _m == 3
							count
							if `r(N)' > 1 {
								rename (v1 v4 v5 v6) (chr geneid start stop)
								split v7,p(" ")
								gen strand = ""
								replace strand = "+" if v72 == ">"
								replace strand = "-" if v71 == "<"
								order geneid chr start stop strand
								gen genename = geneid
								keep geneid chr start stop strand genename
								duplicates drop
								count
								if `r(N)' > 1 {
									egen x = seq(), by(geneid)
									tostring x,replace
									replace geneid = geneid + "." + x
									drop x
									}
								save ``function'_sqc_short'-`i'-``function'_gl_short'.geneloc.dta, replace
								}
							}
						files2dta, dir(`c(pwd)')
						split files,p("gene-list.geneloc.")
						keep if files2 == "dta"
						gen script = ""
						replace script = "append using " + files 
						replace script = "use " + files + ", clear " in 1
						outsheet script using script.do, replace non noq
						do script.do
						drop if chr ==23 
						outsheet using ``function'_sqc_short'-``function'_gl_short'.geneloc, non noq replace	 
						}
					qui { // define gene-based p
						noi di as text"# > map snps to genes (+/- 35,10)"
						!magma --annotate window=35,10 --snp-loc ``function'_sqc_short'.snploc --gene-loc ``function'_sqc_short'-``function'_gl_short'.geneloc --out ``function'_sqc_short'-``function'_gl_short'
						noi di as text"# > define gene-based pval using magma"
						noi !magma 	--bfile 		`ldref' ///
									--pval 		``function'_sqc_short'.pvaloc use=snp,p ncol=n ///
									--gene-model multi ///
									--gene-annot 	``function'_sqc_short'-``function'_gl_short'.genes.annot ///
									--genes-only ///
									--gene-settings adap-permp=1000000,10 ///
									--seed 987654321 ///
									--out 		``function'_sqc_short'-``function'_gl_short'
						!${tabbed} ``function'_sqc_short'-``function'_gl_short'.genes.out
						import delim using ``function'_sqc_short'-``function'_gl_short'.genes.out.tabbed, clear varnames(1)
						lab var gene 			"Gene Name"
						lab var chr  			"Chromosome"
						lab var start 			"Transcript Start (hg19) -35kb"
						lab var stop 			"Transcript End (hg19) +10kb"
						lab var nsnps  		"SNPs in model"
						lab var p_snpwise_mean 	"MAGMA P-val (mean)"
						lab var p_snpwise_top1 	"MAGMA P-val (top)"
						lab var zstat            "Z-Score (gene)"
						lab var permp			"MAGMA P-val (perm)"
						lab var nperm            "N Permutations"
						keep  gene chr start stop nsnps p_snpwise_mean p_snpwise_top1 zstat perm nperm
						order gene chr start stop nsnps p_snpwise_mean p_snpwise_top1 zstat perm nperm
						save ``function'_sqc_short'-``function'_gl_short'-`function'.dta, replace
						export excel ``function'_sqc_short'-`function'.xlsx, firstrow(varlabels) sheet("``function'_gl_short'") sheetreplace	
						noi checkfile, file(`already_created')
						}
					}
				noi di as text"#########################################################################"
				noi di as text""
				}
			}
		}
	qui { // by range
		clear 
		set obs 1
		gen range = "`range'"
		if range == "" {
				}
		else  {
			qui { // plot association in gene region
				noi di as text"#########################################################################"
				noi di as text"# SECTION - 3a: plotting association around "as result "`range'"
				noi di as text"#########################################################################"
				split range, p("chr"":""..")
				local chr 	= range2[1]
				local from 	= range3[1]
				local to 		= range4[1]
				clear
				set obs 1
				gen start = `from'
				gen end = `to'
				gen size = end - start
				gen _sta = .
				gen _end = .
				replace _sta = round(start - 50000,1000)
				replace _end = round(end + 50000,1000)
				local from 	= _sta[1]
				local to 		= _end[1]		
				noi di as text"# > .............. chr "as result"chr`chr'"
				noi di as text"# > .......... from bp "as result"`from'"
				noi di as text"# > ............ to bp "as result"`to'"
				gen range_name = ""
				replace range_name = "`range_name'" if range_name == ""
				replace range_name = "chr`chr'_`from'_`to'" if range_name == ""
				local range_name2 = range_name[1]
				local already_created ``function'_sqc_short'-`range_name2'_graphlocuszoom.png
				capture confirm file `already_created'
				if !_rc {
					noi checkfile, file(`already_created')
					}
				else {
					qui { // plot association around gene
						noi di as text"# > .......... running "as result"chr`chr':`from'..`to'"as text " through graphlocuszoom"
						graphlocuszoom, generef(`generef') ldref(`ldref') recombref(`recombref') summaryqc(`summaryqc') range(chr`chr':`from'..`to') label(on) biotype(all) 
						graph use ${init_root}/graphlocuszoom/``function'_sqc_short'/``function'_sqc_short'-chr`chr'_`from'_`to'_graphlocuszoom.gph
						graph save ``function'_sqc_short'-`range_name2'_graphlocuszoom.gph, replace
						graph export ``function'_sqc_short'-`range_name2'_graphlocuszoom.png, as(png) height(2000) width(6000) replace
						noi checkfile, file(`already_created')
						}
					}
				noi di as text"#########################################################################"
				noi di as text""
				}
			qui { // calculate gene-based association (magma)
				noi di as text"#########################################################################"
				noi di as text"# SECTION - 3b: defining gene-based association for "as result "`range'"
				noi di as text"#########################################################################"
				local already_created ``function'_sqc_short'-`range_name2'-`function'.dta
				capture confirm file `already_created'
				if !_rc {
					noi checkfile, file(`already_created')
					}
				else {
					noi di as text"#########################################################################"
					noi di as text"# SECTION - 3c: run magma by "as result "`range'"
					noi di as text"#########################################################################"
					qui { // define .snploc .pvaloc n
						noi di as text"# > ......... locating files to process"
						noi di as text"# >  parse snploc from "as result"``function'_sqc_short'"
						noi di as text"# >  parse pvaloc from "as result"``function'_sqc_short'"
						use snp chr bp p n using `summaryqc'-summaryqc.dta, clear
						outsheet snp chr bp using ``function'_sqc_short'.snploc, non noq replace
						outsheet snp p n using    ``function'_sqc_short'.pvaloc, noq replace
						}
					qui { // define geneloc
						noi di as text"# >  parse genloc from "as result"`range'"
						clear
						set obs 1
						gen range = "`range'"
						split range, p("chr"":""..")
						gen chr	= range2[1]
						gen start	= range3[1]
						gen stop	= range4[1]
						gen strand = "+"
						gen range_name = "`range_name'"
						replace range_name = "chr`chr'_`from'_`to'" if range_name == ""
						local range_name2 = range_name[1]
						local range_name 3
						duplicates drop 
						outsheet range_name chr start stop strand range using ``function'_sqc_short'-`range_name2'.geneloc, non noq replace
						}
					qui { // define gene-based p
						noi di as text"# > map snps to range (+/- 0,0)"
						!magma --annotate  --snp-loc ``function'_sqc_short'.snploc --gene-loc ``function'_sqc_short'-`range_name2'.geneloc --out ``function'_sqc_short'-`range_name2'
						noi di as text"# > define gene-based pval using magma"
						noi !magma 	--bfile 		`ldref' ///
									--pval 		``function'_sqc_short'.pvaloc use=snp,p ncol=n ///
									--gene-model multi ///
									--gene-annot 	``function'_sqc_short'-`range_name2'.genes.annot ///
									--genes-only ///
									--gene-settings adap-permp=1000000,10 ///
									--seed 987654321 ///
									--out 		``function'_sqc_short'-`range_name2'
						!${tabbed} ``function'_sqc_short'-`range_name2'.genes.out
						import delim using ``function'_sqc_short'-`range_name2'.genes.out.tabbed, clear varnames(1)
						lab var gene 			"Locus"
						lab var chr  			"Chromosome"
						lab var start 			"Start (hg19)"
						lab var stop 			"End (hg19)"
						lab var nsnps  		"SNPs in model"
						lab var p_snpwise_mean 	"MAGMA P-val (mean)"
						lab var p_snpwise_top1 	"MAGMA P-val (top)"
						lab var zstat            "Z-Score (gene)"
						lab var permp			"MAGMA P-val (perm)"
						lab var nperm            "N Permutations"
						rename gene locus
						keep  locus chr start stop nsnps p_snpwise_mean p_snpwise_top1 zstat perm nperm
						order locus chr start stop nsnps p_snpwise_mean p_snpwise_top1 zstat perm nperm
						save ``function'_sqc_short'-`range_name2'-`function'.dta, replace
						export excel ``function'_sqc_short'-`function'.xlsx, firstrow(varlabels) sheet("`range_name2'") sheetreplace	
						noi checkfile, file(`already_created')
						}
					noi di as text"#########################################################################"
					noi di as text""
					}
				}
			noi di as text"#########################################################################"
			noi di as text""
			}
		}	
	qui { // by range-list
		clear 
		set obs 1
		gen range = "`range_list'"
		if range == "" {
				}
		else  {
			qui { // plot association in region
				noi di as text"#########################################################################"
				noi di as text"# SECTION - 3a: plotting association around "as result "`range_list'"
				noi di as text"#########################################################################"
				import delim using `range_list', clear varname(nonames)
				split v1, p("chr"":""..")
				gen chr 	= v12
				gen start = v13
				gen end 	= v14
				for var start end: destring X, replace
				gen size = end - start
				gen _sta = .
				gen _end = .
				replace _sta = round(start - 50000,1000)
				replace _end = round(end + 50000,1000)
				rename v2 range_name
				for var _sta _end: tostring X, replace
				egen x = seq()
				sum x
				global `function'_range_n `r(max)'
				tostring x, replace
				gen local_chr  = "global `function'_r_chr"  + x + " " + chr
				gen local_sta  = "global `function'_r_sta"  + x + " " + _sta
				gen local_end  = "global `function'_r_end"  + x + " " + _end
				gen local_name = "global `function'_r_name" + x + " " + range_name
				foreach i in chr sta end name {
					outsheet local_`i' using script.do, non noq replace
					do script.do
					}
				count
				foreach num of num 1 / `r(N)' {
					capture confirm file ``function'_sqc_short'-${`function'_r_name`num'}_graphlocuszoom.png
					if !_rc {
						noi checkfile, file(``function'_sqc_short'-${`function'_r_name`num'}_graphlocuszoom.png)
						}
					else {
						noi di as text"# > .......... running "as result"chr${`function'_r_chr`num'}:${`function'_r_sta`num'}..${`function'_r_end`num'}"as text " through graphlocuszoom"
						graphlocuszoom, generef(`generef') ldref(`ldref') recombref(`recombref') summaryqc(`summaryqc') range(chr${`function'_r_chr`num'}:${`function'_r_sta`num'}..${`function'_r_end`num'}) label(on) biotype(all)
						graph use ${init_root}/graphlocuszoom/``function'_sqc_short'/``function'_sqc_short'-chr${`function'_r_chr`num'}_${`function'_r_sta`num'}_${`function'_r_end`num'}_graphlocuszoom.gph
						graph save ``function'_sqc_short'-${`function'_r_name`num'}_graphlocuszoom.gph, replace
						graph export ``function'_sqc_short'-${`function'_r_name`num'}_graphlocuszoom.png, as(png) height(2000) width(6000) replace
						noi checkfile, file(``function'_sqc_short'-${`function'_r_name`num'}_graphlocuszoom.png)
						}
					noi di as text"#########################################################################"
					}
				noi di as text""
				}
			qui { // calculate region-based association (magma)
				noi di as text"#########################################################################"
				noi di as text"# SECTION - 3b: defining region-based association for "as result "`range_list'"
				noi di as text"#########################################################################"		
				ref_path_to_short, ref(`range_list')
				local `function'_rl_short ${ref_short}
				local already_created ``function'_sqc_short'-``function'_rl_short'-`function'.dta
				capture confirm file `already_created'
				if !_rc {
					noi checkfile, file(`already_created')
					}
				else {
					noi di as text"#########################################################################"
					noi di as text"# SECTION - 3c: run magma by "as result "``function'_rl_short'"
					noi di as text"#########################################################################"
					qui { // define .snploc .pvaloc n
						noi di as text"# > ......... locating files to process"
						noi di as text"# >  parse snploc from "as result"``function'_sqc_short'"
						noi di as text"# >  parse pvaloc from "as result"``function'_sqc_short'"
						use snp chr bp p n using `summaryqc'-summaryqc.dta, clear
						outsheet snp chr bp using ``function'_sqc_short'.snploc, non noq replace
						outsheet snp p n using    ``function'_sqc_short'.pvaloc, noq replace
						}
					qui { // define geneloc
						import delim using `range_list', clear varname(nonames)
						split v1,p("chr"":""..")
						gen chr	= v12
						gen start	= v13
						gen stop	= v14
						gen strand = "+"
						gen range_name = v2
						outsheet range_name chr start stop strand range using ``function'_sqc_short'-``function'_rl_short'.geneloc, non noq replace			 
						}
					qui { // define gene-based p
						noi di as text"# > map snps to genes (+/- 0,0)"
						!magma --annotate window=0,0 --snp-loc ``function'_sqc_short'.snploc --gene-loc ``function'_sqc_short'-``function'_rl_short'.geneloc --out ``function'_sqc_short'-``function'_rl_short'
						noi di as text"# > define gene-based pval using magma"
						noi !magma 	--bfile 		`ldref' ///
									--pval 		``function'_sqc_short'.pvaloc use=snp,p ncol=n ///
									--gene-model multi ///
									--gene-annot 	``function'_sqc_short'-``function'_rl_short'.genes.annot ///
									--genes-only ///
									--gene-settings adap-permp=1000000,10 ///
									--seed 987654321 ///
									--out 		``function'_sqc_short'-``function'_rl_short'
						!${tabbed} ``function'_sqc_short'-``function'_rl_short'.genes.out
						import delim using ``function'_sqc_short'-``function'_rl_short'.genes.out.tabbed, clear varnames(1)
						lab var gene 			"Locus"
						lab var chr  			"Chromosome"
						lab var start 			"Start (hg19)"
						lab var stop 			"End (hg19)"
						lab var nsnps  		"SNPs in model"
						lab var p_snpwise_mean 	"MAGMA P-val (mean)"
						lab var p_snpwise_top1 	"MAGMA P-val (top)"
						lab var zstat            "Z-Score (gene)"
						lab var permp			"MAGMA P-val (perm)"
						lab var nperm            "N Permutations"
						rename gene locus
						keep  locus chr start stop nsnps p_snpwise_mean p_snpwise_top1 zstat perm nperm
						order locus chr start stop nsnps p_snpwise_mean p_snpwise_top1 zstat perm nperm
						save ``function'_sqc_short'-``function'_rl_short'-`function'.dta, replace
						export excel ``function'_sqc_short'-`function'.xlsx, firstrow(varlabels) sheet("``function'_rl_short'") sheetreplace	
						noi checkfile, file(`already_created')
						}
					}
				noi di as text"#########################################################################"
				noi di as text""
				}
			}
		}
	}
qui { // module 4 - write log file
	file open myfile using "``function'_sqc_short'-`function'.log", write replace
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
	file write myfile `"# > ...... saved data ``function'_sqc_short'-`function'.xlsx"' _n
	file write myfile `"#########################################################################"' _n
	file close myfile	
	}
qui { // module 5 - clean folder
	files2dta, dir(`c(pwd)')
	gen keep = .
	replace keep = 1 if files == "``function'_sqc_short'-`function'.xlsx"
	replace keep = 1 if files == "``function'_sqc_short'-`snp'-`function'.dta"
	replace keep = 1 if files == "``function'_sqc_short'-``function'_sl_short'-`function'.dta"
	replace keep = 1 if files == "``function'_sqc_short'-`gene'-`function'.dta"
	replace keep = 1 if files == "``function'_sqc_short'-``function'_gl_short'-`function'.dta"
	replace keep = 1 if files == "``function'_sqc_short'-range'-`function'.dta"
	replace keep = 1 if files == "``function'_sqc_short'-``function'_rl_short'-`function'.dta"
	split files,p(".log"".png")
	replace keep = 1 if files != files1
	keep files keep
	split files,p("`function'")
	replace keep = 1 if files2 == ".dta"
	drop if keep == 1
	keep files
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

   
