*! 1.0.1 Richard Anney 29jan2020
*! 1.0.2 Richard Anney 5may2020

program  summaryqc2annot
local    function summaryqc2annot
local    `function'_version v2

syntax , summaryqc(string asis) ldref(string asis) generef(string asis) eqtlref(string asis) caddref(string asis) [minp(real 1e-5)]

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
	qui { // define short names
		ref_path_to_short, ref(`summaryqc')
		local `function'_summaryqc_short ${ref_short}
		ref_path_to_short, ref(`ldref')
		local `function'_ldref_short ${ref_short}
		local `function'_glist_short Homo_sapiens.GRCh37.87.gtf
		ref_path_to_short, ref(`eqtlref')
		local `function'_elist_short ${ref_short}
		ref_path_to_short, ref(`caddref')
		local `function'_clist_short ${ref_short}	
		}
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
	foreach file in bim bed fam {
		noi checkfile, file(`ldref'.`file')
		}
	noi checkfile, file(`generef'/Homo_sapiens.GRCh37.87.gtf_exon.dta)
	noi checkfile, file(`eqtlref')
	noi checkfile, file(`caddref')
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 3 - apply annotations
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3: apply annotations"
	noi di as text"#########################################################################"
	noi di as text"# > ......... checking "as result"``function'_summaryqc_short'-summaryqc2top.dta"
	local already_created ${init_root}/summaryqc2top/``function'_summaryqc_short'/``function'_summaryqc_short'-summaryqc2top.dta
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		}
	else {	
		noi summaryqc2top, summaryqc(`summaryqc') ldref(`ldref') minp(`minp')
		noi checkfile, file(`already_created')
		}
	noi di as text"# > ......... checking "as result"``function'_summaryqc_short'-summaryqc2gene.dta"
	local already_created ${init_root}/summaryqc2gene/``function'_summaryqc_short'/``function'_summaryqc_short'-summaryqc2gene.dta
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		}
	else {	
		noi summaryqc2gene, summaryqc(`summaryqc') ldref(`ldref') generef(`generef') minp(`minp')
		noi checkfile, file(`already_created')
		}
	noi di as text"# > ......... checking "as result"``function'_summaryqc_short'-summaryqc2eqtl.dta"
	local already_created ${init_root}/summaryqc2eqtl/``function'_summaryqc_short'/``function'_summaryqc_short'-summaryqc2eqtl.dta
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		}
	else {	
		noi summaryqc2eqtl, summaryqc(`summaryqc') ldref(`ldref') eqtlref(`eqtlref') minp(`minp')
		noi checkfile, file(`already_created')
		}
	noi di as text"# > ......... checking "as result"``function'_summaryqc_short'-summaryqc2cadd.dta"
	local already_created ${init_root}/summaryqc2cadd/``function'_summaryqc_short'/``function'_summaryqc_short'-summaryqc2cadd.dta
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		}
	else {	
		noi summaryqc2cadd, summaryqc(`summaryqc') ldref(`ldref') caddref(`caddref') minp(`minp')
		noi checkfile, file(`already_created')
		}

	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 4 - create summary_excel_sheet
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 4: create Supplementary Tables"
	noi di as text"#########################################################################"	
	local already_created ${init_root}/`function'/``function'_summaryqc_short'/``function'_summaryqc_short'-`function'.xlsx
	qui { // create index sheet
		clear
		set obs 3
		gen n = _n
		qui { // define sheets
			gen sheet     = ""
			replace sheet = "gene_annotation" 			if n == 1
			replace sheet = "eqtl_annotation" 			if n == 2
			replace sheet = "cadd_annotation" 			if n == 3
			drop if sheet == ""
			}
		qui { // add description
			gen description = ""
			replace description = "Annotation of transcripts with physical overlap with the associatied region. Data Source: (ftp://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz)" if sheet == "gene_annotation"
			replace description = "Annotation of cis-eQTL P<1e-5 to SNP and LD partners (R^2 > 0.8). Data Source: GTEx_v7" if sheet == "eqtl_annotation"
			replace description = "Annotation of variant consequence (Consequence, CADD Phred Score, PolyPhen Category and SIFT Category) to SNP and LD partners (R^2 > 0.8) - (max cadd_conscore reported; 5PRIME_UTR(3) REGULATORY(4) NONCODING_CHANGE,SPLIC_SITE,SYNONYMOUS,UNKNOWN(5) CANONICAL_SPLICE(6) FRAME_SHIFT,INFRAME,NON_SYNONYMOUS,STOP_LOST(7) STOP_GAINED(8)); report limited to max. CONSSCORE per gene per index SNP. Data Source: cadd.gs.washington.edu" if sheet == "cadd_annotation"
			}
		qui { // add link
			gen link = ""
			replace link  = `"=HYPERLINK(""' + sheet + `"!A1",""' + sheet + `"")"' 
			}
			export excel sheet description  using ``function'_summaryqc_short'-`function'.xlsx, firstrow(variables) sheet("index") replace
		}
	qui { // create gene_annotation
		use ${init_root}/summaryqc2gene/``function'_summaryqc_short'/``function'_summaryqc_short'-summaryqc2gene.dta, clear
		sort p
		export excel using ``function'_summaryqc_short'-`function'.xlsx, firstrow(variables) sheet("gene_annotation") sheetreplace
		}
	qui { // create eqtl_annotation
		use ${init_root}/summaryqc2eqtl/``function'_summaryqc_short'/``function'_summaryqc_short'-summaryqc2eqtl.dta, clear
		sort p
		export excel using ``function'_summaryqc_short'-`function'.xlsx, firstrow(variables) sheet("eqtl_annotation") sheetreplace
		}
	qui { // create cadd_annotation
		use ${init_root}/summaryqc2cadd/``function'_summaryqc_short'/``function'_summaryqc_short'-summaryqc2cadd.dta, clear
		sort p
		export excel using ``function'_summaryqc_short'-`function'.xlsx, firstrow(variables) sheet("cadd_annotation") sheetreplace
		}
	noi checkfile, file(`already_created')
	noi di as text"#########################################################################"
	noi di as text""	
	}		
		
qui { // module 5 - write log file
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
	file write myfile `"# > .... LD reference ``function'_ldref_short'"' _n
	file write myfile `"# > .. eQTL reference ``function'_elist_short'"' _n
	file write myfile `"# > .. gene reference ``function'_glist_short'"' _n
	file write myfile `"# > .. CADD reference ``function'_clist_short'"' _n
	file write myfile `"# > .. min-P reported `: display %10.4e `minp''"' _n
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
qui { // module 6 - clean folder
	files2dta, dir(`c(pwd)')
	gen keep = .
	replace keep = 1 if file == "``function'_summaryqc_short'-`function'.log"
	replace keep = 1 if file == "``function'_summaryqc_short'-`function'.xlsx"
	
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

   
