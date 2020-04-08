*! 1.0.1 Richard Anney 31mar2020

program  create_generef
local    function create_generef
local    `function'_version v1
syntax  ,  [version(real 87)]

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
	noi di as text"# > ..... project name "as result"Homo_sapiens.GRCh37.`version'.gtf"
	!mkdir  Homo_sapiens.GRCh37.`version'.gtf
	cd      Homo_sapiens.GRCh37.`version'.gtf
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 2 - download and convert master files
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 2: download and convert master files"
	noi di as text"#########################################################################"
	local already_created ${init_root}/`function'/Homo_sapiens.GRCh37.`version'.gtf/Homo_sapiens.GRCh37.`version'.gtf.dta
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
					!wget   ftp://ftp.ensembl.org/pub/grch37/release-`version'/gtf/homo_sapiens/Homo_sapiens.GRCh37.`version'.gtf.gz
					}
				qui { // unpack files
					!gunzip Homo_sapiens.GRCh37.`version'.gtf.gz
					}
				}
			else if os == "Windows" { 	
				qui { // download files
					!bash -c "wget   ftp://ftp.ensembl.org/pub/grch37/release-`version'/gtf/homo_sapiens/Homo_sapiens.GRCh37.`version'.gtf.gz"
					}
				qui { // unpack files
					!bash -c "gunzip Homo_sapiens.GRCh37.`version'.gtf.gz"
					}
				}		
			import delim using Homo_sapiens.GRCh37.`version'.gtf, clear varnames(noname) rowrange(6:)
			rename (v1 - v9) (chr source feature start end score strand frame attribute)
			replace chr = "23" if chr == "X"
			replace chr = "24" if chr == "Y"
			replace chr = "25" if chr == "XY"
			replace chr = "26" if chr == "MT"
			destring chr, replace force
			keep chr source feature start end score strand frame attribute
			save Homo_sapiens.GRCh37.`version'.gtf.dta, replace
			noi checkfile, file(`already_created')
			}
	noi di as text"#########################################################################"
	noi di as text""			
	}
qui { // module 3 - process gene files
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3: process gene files"
	noi di as text"#########################################################################"
	local already_created ${init_root}/`function'/Homo_sapiens.GRCh37.`version'.gtf/Homo_sapiens.GRCh37.`version'.gtf_gene.dta
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		}
	else {
		foreach chr of num 1/23 {
			use Homo_sapiens.GRCh37.`version'.gtf.dta if chr == `chr' & feature == "gene", clear
			split attribute,p("; ")
			drop attribute
			gen ensembl_geneID = subinstr(attribute1, "gene_id ", "",.) 
			gen symbol         = subinstr(attribute3, "gene_name ", "",.)
			gen biotype        = subinstr(attribute5, "gene_biotype ", "",.)
			foreach i in ensembl symbol biotype { 
				replace `i' = subinstr(`i', `"""',"",.)
				replace `i' = subinstr(`i', ";","",.)
				}
			keep  ensembl_geneID chr start end strand symbol biotype
			order ensembl_geneID chr start end strand symbol biotype
			sort  ensembl_geneID
			compress
			gen label = ""
			replace label =        symbol + " >" if strand == "+"
			replace label = "< " + symbol        if strand == "-"
			compress
			save Homo_sapiens.GRCh37.`version'.chr`chr'-gene.gtf.dta, replace
			}
		foreach chr of num 1/22 {
			append using Homo_sapiens.GRCh37.`version'.chr`chr'-gene.gtf.dta
			}
		duplicates drop
		sort chr start
		lab var ensembl_geneID 	"gene_id (ensembl): The stable identifier for the gene"
		lab var chr 						"chromosome"
		lab var start 					"start position of the feature, with sequence numbering starting at 1"
		lab var end 						"end position of the feature, with sequence numbering starting at 1."
		lab var strand 					"defined as + (forward) or - (reverse)"
		lab var symbol 					"official symbol of this gene"
		lab var biotype 				"biotype of this gene"
		lab var label 					"symbol + strand label (+ = >; - = <)"  
		save Homo_sapiens.GRCh37.`version'.gtf_gene.dta, replace
		noi checkfile, file(`already_created')
		foreach chr of num 1/23 {
			erase Homo_sapiens.GRCh37.`version'.chr`chr'-gene.gtf.dta
			}
		}
	noi di as text"#########################################################################"
	noi di as text""		
	}
qui { // module 4 - process exon files
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 4: process exon files"
	noi di as text"#########################################################################"
	local already_created ${init_root}/`function'/Homo_sapiens.GRCh37.`version'.gtf/Homo_sapiens.GRCh37.`version'.gtf_exon.dta
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		}
	else {
		foreach chr of num 1/23 {
			use Homo_sapiens.GRCh37.`version'.gtf.dta if chr == `chr' & feature == "exon", clear
			split attribute,p("; ")
			drop attribute
			gen ensembl_geneID  = subinstr(attribute1, "gene_id ", "",.) 
			gen ensembl_txID    = subinstr(attribute3, "transcript_id ", "",.) 
			gen exonCount       = subinstr(attribute5, "exon_number ", "",.)
			gen symbol          = subinstr(attribute6, "gene_name ", "",.)
			gen ensembl_geneID2 = subinstr(attribute9, "transcript_name ", "",.)
			foreach i in ensembl_geneID ensembl_txID ensembl_geneID2 exonCount symbol { 
				replace `i' = subinstr(`i', `"""',"",.)
				replace `i' = subinstr(`i', ";","",.)
				}
			for var start end: rename X txX
			keep  ensembl_geneID ensembl_geneID2 ensembl_txID exonCount txstart txend symbol
			order ensembl_geneID ensembl_geneID2 ensembl_txID exonCount txstart txend symbol
			sort symbol
			merge m:m symbol using ${init_root}/`function'/Homo_sapiens.GRCh37.`version'.gtf/Homo_sapiens.GRCh37.`version'.gtf_gene.dta
			keep if _m == 3
			keep if chr == `chr'
			keep  symbol ensembl_geneID ensembl_txID exonCount chr start end txstart txend symbol ensembl_geneID2 biotype label
			order ensembl_geneID chr start end symbol biotype
			sort  chr start end 
			compress
			save Homo_sapiens.GRCh37.`version'.chr`chr'-exon.gtf.dta, replace
			}
		foreach chr of num 1/22 {
			append using Homo_sapiens.GRCh37.`version'.chr`chr'-exon.gtf.dta
			}
		duplicates drop
		sort chr start
		lab var ensembl_geneID 	"gene_id (ensembl): The stable identifier for the gene"
		lab var chr 						"chromosome"
		lab var start 					"start position of the gene, with sequence numbering starting at 1"
		lab var end 						"end position of the gene, with sequence numbering starting at 1."
		lab var symbol 					"official symbol of this gene"
		lab var biotype 				"biotype of this gene"
		lab var label 					"symbol + strand label (+ = >; - = <)"  
		lab var ensembl_geneID2 "gene_id (ensembl): The stable identifier for the gene + version"
    lab var ensembl_txID 		"The stable identifier for this transcript"               
		lab var exonCount    		"exon count"                   
		lab var txstart      		"start position of the exon, with sequence numbering starting at 1"                
		lab var txend        		"end position of the exon, with sequence numbering starting at 1"               
	
		
		save Homo_sapiens.GRCh37.`version'.gtf_exon.dta, replace
		noi checkfile, file(`already_created')
		foreach chr of num 1/23 {
			erase Homo_sapiens.GRCh37.`version'.chr`chr'-exon.gtf.dta
			}
		noi checkfile, file(`already_created')
		}			
	noi di as text"#########################################################################"
	noi di as text""		
	}
qui { // module 5 - create hg-list files
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 5: create hg-list files"
	noi di as text"#########################################################################"
	use  Homo_sapiens.GRCh37.`version'.gtf_gene.dta, clear
			for var start end: tostring X, replace
			foreach biotype in _3prime_overlapping_ncrna IG_C_gene IG_C_pseudogene IG_D_gene IG_J_gene IG_J_pseudogene IG_V_gene IG_V_pseudogene Mt_rRNA Mt_tRNA TR_C_gene TR_D_gene TR_J_gene TR_J_pseudogene TR_V_gene antisense lincRNA miRNA misc_RNA polymorphic_pseudogene processed_transcript protein_coding pseudogene rRNA sense_intronic sense_overlapping snRNA snoRNA {
				outsheet chr start end symbol if biotype == "`biotype'" using Homo_sapiens.GRCh37.`version'.gtf-`biotype'.hg-list, non noq replace
				}
			use  Homo_sapiens.GRCh37.`version'.gtf_exon.dta, clear
			for var txstart txend start end: tostring X, replace
			foreach biotype in _3prime_overlapping_ncrna IG_C_gene IG_C_pseudogene IG_D_gene IG_J_gene IG_J_pseudogene IG_V_gene IG_V_pseudogene Mt_rRNA Mt_tRNA TR_C_gene TR_D_gene TR_J_gene TR_J_pseudogene TR_V_gene antisense lincRNA miRNA misc_RNA polymorphic_pseudogene processed_transcript protein_coding pseudogene rRNA sense_intronic sense_overlapping snRNA snoRNA {
				outsheet chr txstart txend symbol start end label if biotype == "`biotype'" using Homo_sapiens.GRCh37.`version'.gtf_exon-`biotype'.hg-list, non noq replace
				}
	noi di as text"#########################################################################"
	noi di as text""							
	}
qui { //
	cd `return'
	}
restore
noi di as text"#########################################################################"
noi di as text"# Completed: $S_DATE $S_TIME"
noi di as text"#########################################################################"
end;
