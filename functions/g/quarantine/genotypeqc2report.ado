*! 1.0.9 Richard Anney 10dec2018
* - update bim2path to set folder - 1.0.9 - 10dec2018 - Richard Anney

program  genotypeqc2report
local    function genotypeqc2report
local    `function'_version v9
syntax,  bim(string asis)

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
	bim_path_to_short, bim(`bim')
	local `function'_bim_short ${bim_short}
	ref_path_to_short, ref($bim2frq_compare_ref)
	local `function'_bim2frq_short ${ref_short}
	cd ${init_root}	
	!mkdir `function'
	cd `function'
	noi di as text"# > ..... project name "as result"``function'_bim_short'"
	!mkdir ``function'_bim_short'
	cd ``function'_bim_short'
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 2 - define parameters from log file
		noi di as text"#########################################################################"
		noi di as text"# SECTION - 2: define parameters from log file"
		noi di as text"#########################################################################"
		noi checkfile, file(${init_root}/genotypeqc/``function'_bim_short'/``function'_bim_short'-genotypeqc.log)
		import delim using  ${init_root}/genotypeqc/``function'_bim_short'/``function'_bim_short'-genotypeqc.log, clear
		gen threshold = substr(v1,23,50)
		gen parameter = "genotypeqc2report_"
		replace parameter = parameter + "plink" in 14
		replace parameter = parameter + "plink2" in 15
		replace parameter = parameter + "infile" in 17
		replace parameter = parameter + "insnp" in 18
		replace parameter = parameter + "inind" in 19
		replace parameter = parameter + "array" in 20
		replace parameter = parameter + "jaccard" in 21
		replace parameter = parameter + "inbuild" in 22
		replace parameter = parameter + "outfile" in 24
		replace parameter = parameter + "outsnp" in 25
		replace parameter = parameter + "outind" in 26
		replace parameter = parameter + "outbuild" in 27
		replace parameter = parameter + "refid" in 28
		replace parameter = parameter + "frqcompare" in 29
		replace parameter = parameter + "rounds" in 33
		replace parameter = parameter + "mac" in 34
		replace parameter = parameter + "geno" in 35
		replace parameter = parameter + "mind" in 36
		replace parameter = parameter + "hardy" in 37
		replace parameter = parameter + "het" in 38
		replace parameter = parameter + "kind" in 40
		replace parameter = parameter + "kinf" in 41
		replace parameter = parameter + "kins" in 42
		replace parameter = parameter + "kint" in 43
		replace parameter = parameter + "ancestry" in 44
		gen global = "global " + parameter + " " + threshold
		drop if parameter == "genotypeqc2report_"
		outsheet global using temp.do, non noq replace
		do temp.do
		erase temp.do
		noi di as text"#########################################################################"
		noi di as text""
		}
qui { // module 3 - define files
		noi di as text"#########################################################################"
		noi di as text"# SECTION - 3: define files"
		noi di as text"#########################################################################"
		foreach file in bim bed fam {
			noi checkfile, file(${init_root}/genotypeqc/``function'_bim_short'/``function'_bim_short'-qc-``function'_version'.`file')
			}
		noi checkfile, file(${init_root}/bim2frq_compare/``function'_bim_short'-mac5-hg19-refid/``function'_bim_short'-mac5-hg19-refid-``function'_bim2frq_short'_bim2frq_compare.gph)
		noi checkfile, file(${init_root}/bim2graphfrq/``function'_bim_short'-mac5-hg19-refid-bim2frq/``function'_bim_short'-mac5-hg19-refid-bim2frq_bim2graphfrq.gph)
		noi checkfile, file(${init_root}/bim2graphhet/``function'_bim_short'-mac5-hg19-refid-bim2frq/``function'_bim_short'-mac5-hg19-refid-bim2frq_bim2graphhet.gph)
		noi checkfile, file(${init_root}/bim2graphhwe/``function'_bim_short'-mac5-hg19-refid-bim2frq/``function'_bim_short'-mac5-hg19-refid-bim2frq_bim2graphhwe.gph)
		noi checkfile, file(${init_root}/bim2graphkin0/``function'_bim_short'-mac5-hg19-refid-bim2frq/``function'_bim_short'-mac5-hg19-refid-bim2frq_bim2graphkin0-scatter.gph)
		noi checkfile, file(${init_root}/bim2graphkin0/``function'_bim_short'-mac5-hg19-refid-bim2frq/``function'_bim_short'-mac5-hg19-refid-bim2frq_bim2graphkin0-hist.gph)
		noi checkfile, file(${init_root}/bim2graphmiss/``function'_bim_short'-mac5-hg19-refid-bim2frq/``function'_bim_short'-mac5-hg19-refid-bim2frq_bim2graphmiss-imiss.gph)
		noi checkfile, file(${init_root}/bim2graphmiss/``function'_bim_short'-mac5-hg19-refid-bim2frq/``function'_bim_short'-mac5-hg19-refid-bim2frq_bim2graphmiss-lmiss.gph)
		noi checkfile, file(${init_root}/bim2graphfrq/``function'_bim_short'-qc-``function'_version'/``function'_bim_short'-qc-``function'_version'_bim2graphfrq.gph)
		noi checkfile, file(${init_root}/bim2graphhet/``function'_bim_short'-qc-``function'_version'/``function'_bim_short'-qc-``function'_version'_bim2graphhet.gph)
		noi checkfile, file(${init_root}/bim2graphhwe/``function'_bim_short'-qc-``function'_version'/``function'_bim_short'-qc-``function'_version'_bim2graphhwe.gph)
		noi checkfile, file(${init_root}/bim2graphkin0/``function'_bim_short'-qc-``function'_version'/``function'_bim_short'-qc-``function'_version'_bim2graphkin0-scatter.gph)
		noi checkfile, file(${init_root}/bim2graphkin0/``function'_bim_short'-qc-``function'_version'/``function'_bim_short'-qc-``function'_version'_bim2graphkin0-hist.gph)
		noi checkfile, file(${init_root}/bim2graphmiss/``function'_bim_short'-qc-``function'_version'/``function'_bim_short'-qc-``function'_version'_bim2graphmiss-imiss.gph)
		noi checkfile, file(${init_root}/bim2graphmiss/``function'_bim_short'-qc-``function'_version'/``function'_bim_short'-qc-``function'_version'_bim2graphmiss-lmiss.gph)
		noi checkfile, file(${init_root}/bim2ancestry/``function'_bim_short'-qc-``function'_version'/``function'_bim_short'-qc-``function'_version'-bim2ancestry_pca.gph)
		noi checkfile, file(${init_root}/bim2ancestry/``function'_bim_short'-qc-``function'_version'/``function'_bim_short'-qc-``function'_version'-bim2ancestry_pca1-pca2.gph)
		noi checkfile, file(${init_root}/bim2ancestry/``function'_bim_short'-qc-``function'_version'/``function'_bim_short'-qc-``function'_version'-bim2ancestry_pca1-pca3.gph)
		noi checkfile, file(${init_root}/bim2ancestry/``function'_bim_short'-qc-``function'_version'/``function'_bim_short'-qc-``function'_version'-bim2ancestry_pca2-pca3.gph)
		noi di as text"#########################################################################"
		noi di as text""	
		}
qui { // module 4 - create report graphs
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 4: create report graphs"
	noi di as text"#########################################################################"
	qui { // snp-by-chromosome
		local already_created ${init_root}/genotypeqc2report/${genotypeqc2report_infile}/${genotypeqc2report_outfile}-snp-by-chromosome.png
		capture confirm file `already_created'
		if !_rc {
			noi di as text"# > .................. plot markers per chromosome"
			noi checkfile, file(`already_created')
			}
		else {		
			qui { // ${genotypeqc2report_infile}
				import delim  using ${init_root}/genotypeqc/${genotypeqc2report_infile}/${genotypeqc2report_infile}.bim, clear
				capture confirm variable v2
				if !_rc { 
					}
				else {
					rename v1 v
					split v,p(" ")
					drop v
					}
				rename (v1 v2 v4 v5 v6) (chr snp bp a1 a2)
				destring chr, replace
				sum chr
				noi di as text"# > .................. plot markers per chromosome for ${genotypeqc2report_infile}"
				hist chr,  xlabel(1(1)`r(max)') title("${genotypeqc2report_infile}") xtitle("Chromosome") caption("count based on ${genotypeqc2report_insnp} SNPs") discrete freq ylabel(#4,format(%9.0g)) nodraw saving(${genotypeqc2report_infile}-snp-by-chromosome.gph, replace)
				}
			qui { // ${genotypeqc2report_outfile}
				local already_created ${init_root}/bim2dta/${genotypeqc2report_outfile}/${genotypeqc2report_outfile}_bim.dta
				capture confirm file `already_created'
				if !_rc {
					noi checkfile, file(`already_created')
					}
				else {	
					bim2dta, bim(${init_root}/genotypeqc/${genotypeqc2report_infile}/${genotypeqc2report_outfile})
					noi checkfile, file(`already_created')
					}
				destring chr, replace
				sum chr
				noi di as text"# > .................. plot markers per chromosome for ${genotypeqc2report_outfile}"
				hist chr,  xlabel(1(1)`r(max)') title("${genotypeqc2report_outfile}") xtitle("Chromosome") caption("count based on ${genotypeqc2report_outsnp} SNPs") discrete freq ylabel(#4,format(%9.0g)) nodraw saving(${genotypeqc2report_outfile}-snp-by-chromosome.gph, replace)
				}
			noi di as text"# > .................. join plots of markers per chromosome"
			graph combine ${genotypeqc2report_infile}-snp-by-chromosome.gph ${genotypeqc2report_outfile}-snp-by-chromosome.gph, col(1) 
			graph export  ${genotypeqc2report_outfile}-snp-by-chromosome.png, as(png) replace width(8000) height(4000)
			erase ${genotypeqc2report_infile}-snp-by-chromosome.gph 
			erase ${genotypeqc2report_outfile}-snp-by-chromosome.gph
			window manage close graph
			checkfile, file(${genotypeqc2report_outfile}-snp-by-chromosome.png)
			}
		}
	qui { // combine quality control plots 
		foreach i in frq het hwe {
			local already_created ${init_root}/genotypeqc2report/${genotypeqc2report_infile}/${genotypeqc2report_outfile}_bim2graph`i'.png
			capture confirm file `already_created'
			if !_rc {
				noi di as text"# > .................. join bim2graph`i'"
				noi checkfile, file(`already_created')
				}
			else {		
				noi di as text"# > .................. join bim2graph`i'"
				graph combine ${init_root}/bim2graph`i'/${genotypeqc2report_infile}-mac5-hg19-refid-bim2frq/${genotypeqc2report_infile}-mac5-hg19-refid-bim2frq_bim2graph`i'.gph ///
				              ${init_root}/bim2graph`i'/${genotypeqc2report_outfile}/${genotypeqc2report_outfile}_bim2graph`i'.gph, col(1)  
				graph export  ${genotypeqc2report_outfile}_bim2graph`i'.png, as(png) replace width(8000) height(4000)
				window manage close graph
				checkfile, file(${genotypeqc2report_outfile}_bim2graph`i'.png)
				}
			}
		foreach i in hist scatter {
			local already_created ${init_root}/genotypeqc2report/${genotypeqc2report_infile}/${genotypeqc2report_outfile}_bim2graphkin0-`i'.png
			capture confirm file `already_created'
			if !_rc {
				noi di as text"# > .................. join bim2graphkin0-`i'"
				noi checkfile, file(`already_created')
				}
			else {		
				noi di as text"# > .................. join bim2graphkin0-`i'"
				graph combine ${init_root}/bim2graphkin0/${genotypeqc2report_infile}-mac5-hg19-refid-bim2frq/${genotypeqc2report_infile}-mac5-hg19-refid-bim2frq_bim2graphkin0-`i'.gph ///
				              ${init_root}/bim2graphkin0/${genotypeqc2report_outfile}/${genotypeqc2report_outfile}_bim2graphkin0-`i'.gph, col(1) 
				graph export  ${genotypeqc2report_outfile}_bim2graphkin0-`i'.png, as(png) replace width(8000) height(4000)
				window manage close graph
				checkfile, file(${genotypeqc2report_outfile}_bim2graphkin0-`i'.png)
				}
			}
		foreach i in imiss lmiss {
			local already_created ${init_root}/genotypeqc2report/${genotypeqc2report_infile}/${genotypeqc2report_outfile}_bim2graphmiss-`i'.png
			capture confirm file `already_created'
			if !_rc {
				noi di as text"# > .................. join bim2graphmiss-`i'"
				noi checkfile, file(`already_created')
				}
			else {		
				noi di as text"# > .................. join bim2graphmiss-`i'"
				graph combine ${init_root}/bim2graphmiss/${genotypeqc2report_infile}-mac5-hg19-refid-bim2frq/${genotypeqc2report_infile}-mac5-hg19-refid-bim2frq_bim2graphmiss-`i'.gph ///
				              ${init_root}/bim2graphmiss/${genotypeqc2report_outfile}/${genotypeqc2report_outfile}_bim2graphmiss-`i'.gph, col(1) 
				graph export  ${genotypeqc2report_outfile}_bim2graphmiss-`i'.png, as(png) replace width(8000) height(4000)
				window manage close graph
				checkfile, file(${genotypeqc2report_outfile}_bim2graphmiss-`i'.png)
				}
			}
		foreach i in bim2frq_compare {
			local already_created ${init_root}/bim2frq_compare/``function'_bim_short'-mac5-hg19-refid/``function'_bim_short'-mac5-hg19-refid-``function'_bim2frq_short'_`i'.png
			capture confirm file `already_created'
			if !_rc {
				noi di as text"# > .................. convert `i'"
				noi checkfile, file(`already_created')
				}
			else {		
				noi di as text"# > .................. convert `i'"
				graph combine ${init_root}/bim2frq_compare/``function'_bim_short'-mac5-hg19-refid/``function'_bim_short'-mac5-hg19-refid-``function'_bim2frq_short'_bim2frq_compare.gph
				graph export  ${genotypeqc2report_outfile}_`i'.png, as(png) replace width(8000) height(4000)
				window manage close graph
				checkfile, file(${genotypeqc2report_outfile}_`i'.png)
				}
			}
		foreach i in bim2ancestry {
			local already_created ${init_root}/genotypeqc2report/${genotypeqc2report_infile}/${genotypeqc2report_outfile}_`i'_pca.png
			capture confirm file `already_created'
			if !_rc {
				noi di as text"# > .................. convert `i'_pca"
				noi checkfile, file(`already_created')
				}
			else {		
				noi di as text"# > .................. convert `i'_pca"
				graph combine ${init_root}/`i'/${genotypeqc2report_outfile}/${genotypeqc2report_outfile}-`i'_pca.gph, title("test data: ${genotypeqc2report_outfile}""reference data: ${genotypeqc2report_ancestry}")
				graph export  ${genotypeqc2report_outfile}_`i'_pca.png, as(png) replace width(8000) height(4000)
				window manage close graph
				checkfile, file(${genotypeqc2report_outfile}_`i'_pca.png)
				}
			local already_created ${init_root}/genotypeqc2report/${genotypeqc2report_infile}/${genotypeqc2report_outfile}_`i'_pca1-pca2.png
			capture confirm file `already_created'
			if !_rc {
				noi di as text"# > .................. convert `i'_pca1-pca2"
				noi checkfile, file(`already_created')
				}
			else {		
				noi di as text"# > .................. convert `i'_pca1-pca2"
				graph combine ${init_root}/`i'/${genotypeqc2report_outfile}/${genotypeqc2report_outfile}-`i'_pca1-pca2.gph, title("test data: ${genotypeqc2report_outfile}""reference data: ${genotypeqc2report_ancestry}") 
				graph export  ${genotypeqc2report_outfile}_`i'_pca1-pca2.png, as(png) replace width(8000) height(4000)
				window manage close graph
				checkfile, file(${genotypeqc2report_outfile}_`i'_pca1-pca2.png)
				}
			local already_created ${init_root}/genotypeqc2report/${genotypeqc2report_infile}/${genotypeqc2report_outfile}_`i'_pca1-pca3.png
			capture confirm file `already_created'
			if !_rc {
				noi di as text"# > .................. convert `i'_pca1-pca3"
				noi checkfile, file(`already_created')
				}
			else {		
				noi di as text"# > .................. convert `i'_pca1-pca3"
				graph combine ${init_root}/`i'/${genotypeqc2report_outfile}/${genotypeqc2report_outfile}-`i'_pca1-pca3.gph, title("test data: ${genotypeqc2report_outfile}""reference data: ${genotypeqc2report_ancestry}") 
				graph export  ${genotypeqc2report_outfile}_`i'_pca1-pca3.png, as(png) replace width(8000) height(4000)
				window manage close graph
				checkfile, file(${genotypeqc2report_outfile}_`i'_pca1-pca3.png)
				}
			local already_created ${init_root}/genotypeqc2report/${genotypeqc2report_infile}/${genotypeqc2report_outfile}_`i'_pca2-pca3.png
			capture confirm file `already_created'
			if !_rc {
				noi di as text"# > .................. convert `i'_pca2-pca3"
				noi checkfile, file(`already_created')
				}
			else {		
				noi di as text"# > .................. convert `i'_pca2-pca3"
				graph combine ${init_root}/`i'/${genotypeqc2report_outfile}/${genotypeqc2report_outfile}-`i'_pca2-pca3.gph, title("test data: ${genotypeqc2report_outfile}""reference data: ${genotypeqc2report_ancestry}") 
				graph export  ${genotypeqc2report_outfile}_`i'_pca2-pca3.png, as(png) replace width(8000) height(4000)
				window manage close graph
				checkfile, file(${genotypeqc2report_outfile}_`i'_pca2-pca3.png)
				}
				
			}
		}
	noi di as text"#########################################################################"
	}
qui { // module 5 - create report
	noi di as text""
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 5: create report"
	noi di as text"#########################################################################"
	genotypeqc2report_subroutine
	noi checkfile, file(${init_root}/genotypeqc/${genotypeqc2report_infile}/${genotypeqc2report_outfile}-quality-control-report.docx)
	noi di as text"#########################################################################"
	noi di as text""
	}
qui cd `return'
restore

noi di as text"#########################################################################"
noi di as text"# Completed: $S_DATE $S_TIME"
noi di as text"#########################################################################"
end;
