*! 1.0.2 Richard Anney 30apr2019
* clean up-front functions
* change output name to include chr and loc
*! 1.0.3 Richard Anney 6feb202
* add bedtools wrapper to and clean up order code
*! 1.0.4 Richard Anney 6feb202
* remove bedtools to make it compatible with windows

program  graphgene
local    function graphgene
local    `function'_version v3
syntax  ,  range(string asis) generef(string asis) [label(string asis) png(string asis) biotype(string asis)]

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
	cd ${init_root}	
	!mkdir `function'
	cd `function'
	noi di as text"# > ..... project name "as result"``function'_generef_short'"
	!mkdir ``function'_generef_short'
	cd ``function'_generef_short'
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 2 - define files
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 2: define files"
	noi di as text"#########################################################################"
	noi checkfile, file(`generef'/Homo_sapiens.GRCh37.87.gtf_exon.dta)
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 3 - define region to plot
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3: define region to plot"
	noi di as text"#########################################################################"
	clear
	set obs 1
	gen `function'_range = "`range'"
	split `function'_range,p("chr"":""..")
	local chr 	= `function'_range2[1]
	local from 	= `function'_range3[1]
	local to 		= `function'_range4[1]
	noi di as text"# > .............. chr "as result"`chr'"
	noi di as text"# > .......... from bp "as result"`from'"
	noi di as text"# > ............ to bp "as result"`to'"
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 4 - define genes in region
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 4: define genes in region"
	noi di as text"#########################################################################"
	noi di as text"# > ......... define genes in region "
	local already_created ${init_root}/intersect_genes/``function'_generef_short'/``function'_generef_short'-with-chr`chr'_`from'_`to'_intersect_genes.dta
	capture confirm file `already_created'
	if !_rc {
		}
	else {
		intersect_genes, range(`range') generef(`generef') 
		}
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 5 - assign position to plot
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 5: assign where to plot"
	noi di as text"#########################################################################"
	use `already_created', clear
	keep chr start end label
	duplicates drop
	save `function'-chr`chr'_`from'_`to'-genes_in_region.dta,replace
	use  `already_created', clear
	replace symbol = "-" if biotype == "dummy"
	egen group = group (symbol start end)
	sum group
	global group `r(max)'
	noi di as text"# > elements in region "as result"`r(max)'"
	gen class = ""
	lab var class "transcript class for plotting"
	replace class = "protein_coding" if biotype == "protein_coding"
	foreach i in antisense lincRNA miRNA misc_RNA pseudogene snRNA snoRNA {
		replace class = "other(1)" if biotype == "`i'"
		}
	foreach i in _3prime_overlapping_ncrna IG_C_gene IG_C_pseudogene IG_D_gene IG_J_gene IG_J_pseudogene IG_V_gene IG_V_pseudogene Mt_rRNA Mt_tRNA TR_C_gene TR_D_gene TR_J_gene TR_J_pseudogene TR_V_gene  polymorphic_pseudogene processed_transcript protein_coding rRNA sense_intronic sense_overlapping {
		replace class = "other(2)" if biotype == "`i'"
		}
	qui { // crop boundaries
		gen exclude = ""
		qui { // transcripts within region
			replace exclude = "`range'" if chr == `chr' & inrange(txstart,`from',`to')
			replace exclude = "`range'" if chr == `chr' & inrange(txend,`from',`to')
			}
		qui { // transcripts the encompass the region
			replace exclude = "`range'" if chr == `chr' & txstart <`from' & txend >`to'
			}		
		drop if exclude == ""
		drop exclude
		merge m:1  chr start end label using `function'-chr`chr'_`from'_`to'-genes_in_region.dta
		drop _m
		}
	qui { // define plotting criteria
		gen label_plot = "`label'"
		lab var label_plot "label for plotting"
		if label_plot == "off" {
			replace label = ""
			noi di as text"# > . plot gene labels "as result"no"
			}
		else {
			noi di as text"# > . plot gene labels "as result"yes"
			}
		gen biotype_plot = "`biotype'"
		lab var biotype_plot "biotype for plotting"
		if biotype_plot == "all" {
			noi di as text"# > .... plot biotypes "as result"all"
			}
		else if biotype_plot == "protein_coding" {
			noi di as text"# > .... plot biotypes "as result"protein_coding only"
			drop if biotype != "protein_coding" & biotype != "dummy"
			}
		else {
			noi di as text"# > .... plot biotypes "as result"protein_coding only (default : set to all)"
			drop if biotype != "protein_coding" & biotype != "dummy"
			}
		}
	qui { // define plotting order
		gene_plotorder, from(`from') to(`to')
		}
	qui { // define label position
		gen pos = 1
		lab var pos "position of label"
		replace pos = 11  if start < `from' + 10000
		replace pos = 1 if start == `from'
		}
	qui { // expand txstart txend
	gen range_siz = `to'-`from'
	sum range_siz
	if `r(max)' < 10000 {	
			gen txstart_ = txstart
			gen txend_   = txend
			}
	else if `r(max)' < 100000 {
			gen txstart_ = txstart  - 100
			gen txend_   = txend   + 100
			}
	else if `r(max)' < 1000000 {
			gen txstart_ = txstart  - 200
			gen txend_   = txend   + 200
			}
	else  {
			gen txstart_ = txstart  - 400
			gen txend_   = txend   + 400
			}
	lab var txstart_ "plotting position for transcription start"
	lab var txend_ "plotting position for transcription end"
	drop range_siz
		}
	save ``function'_generef_short'-chr`chr'_`from'_`to'_`function'.dta, replace
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 6 - plot genes
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 6: plot genes"
	noi di as text"#########################################################################"	
	qui { // plot graph
		count
		sum order
		if `r(N)' > 1 {
			#delim;
			twoway rspike  start end 		order if class == "protein_coding", hor lcolor(green) lwidth(vvthin) 
			|| rspike  txstart_ txend_ 	order if class == "protein_coding", hor lcolor(green) lwidth(${gene_plotorder_e}) 
			|| rspike  start end 		order if class == "other(1)", hor lcolor(red) lwidth(vvthin) 
			|| rspike  txstart_ txend_ 	order if class == "other(1)", hor lcolor(red) lwidth(${gene_plotorder_e}) 
			|| rspike  start end 		order if class == "other(2)", hor lcolor(blue) lwidth(vvthin) 
			|| rspike  txstart_ txend_ 	order if class == "other(2)", hor lcolor(blue) lwidth(${gene_plotorder_e}) 						
			|| scatter order start if pos == 11  , msymbol(i) mlabel(label) mlabpos(9) mlabcolor(black) mlabsize(${gene_plotorder_f}) 
			|| scatter order end   if pos == 1   , msymbol(i) mlabel(label) mlabpos(3) mlabcolor(black) mlabsize(${gene_plotorder_f}) 
			legend(off) 
			ylab(1(1)`r(max)',labc(white) nogrid) ytitle("genes")
			xlab(,format(%10.0f) nogrid)
			xtitle("Chromosome `chr'") caption("source: ``function'_generef_short'", size(small))   
			graphregion(margin(zero)) nodraw
			saving(``function'_generef_short'-chr`chr'_`from'_`to'_`function'.gph, replace) 
			;
			#delim cr
			}
		else {
			noi di as text"# > .................. "as error"no transcripts in region - plotting as blank"
			tw scatteri 1 1, msymbol(i) ylab("") xlab("") ytitle("") xtitle("") yscale(off) xscale(off) plotregion(lpattern(blank)) nodraw saving(``function'_generef_short'-chr`chr'_`from'_`to'_`function'.gph, replace)   
			noi di as text"# > ......... saved to "as result"``function'_generef_short'-chr`chr'_`from'_`to'_`function'.gph"
			}
		}
	qui { // export to png
		clear 
		set obs 1
		gen png = "`png'"
		if png == "yes" {
			graph use ``function'_generef_short'-chr`chr'_`from'_`to'_`function'.gph
			noi di as text"# > ..... saved png to "as result"``function'_generef_short'-chr`chr'_`from'_`to'_`function'.png"
			graph export ``function'_generef_short'-chr`chr'_`from'_`to'_`function'.png, as(png) height(1000) width(3000) replace
			window manage close graph
			}
		else {
			}
		}
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 7 - write log file
	file open myfile using "``function'_generef_short'-chr`chr'_`from'_`to'_`function'.log", write replace
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
	file write myfile `"# > ..... project name ``function'_generef_short'-chr`chr'_`from'_`to'_`function'"' _n
	file write myfile `"# >  co-ordinates file ``function'_generef_short'"' _n
	file write myfile `"# > ..... genome build hg19"' _n
	file write myfile `"# > .............. chr chr`chr'"' _n
	file write myfile `"# > ............. from `from'"' _n
	file write myfile `"# > ............... to `to'"' _n
	file write myfile `"# > ..... saved graphs ``function'_generef_short'-chr`chr'_`from'_`to'_`function'.gph"' _n
	file write myfile `"# > ....... saved data ``function'_generef_short'-chr`chr'_`from'_`to'_`function'.dta"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > overview of `function'"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > 1. use reference co-ordinates"' _n
	file write myfile `"# > 2. select protein coding transcripts (adjustable in code)"' _n
	file write myfile `"# > 3. plot genes within defined boundaries"' _n
	file write myfile `"#########################################################################"' _n
	file close myfile	
	}
qui { // module 8 - clean folder
	files2dta, dir(`c(pwd)')
	split files,p("`function'")
	gen keep = .
	replace keep = 1 if files2 == ".png"
	replace keep = 1 if files2 == ".gph"
	replace keep = 1 if files2 == ".dta"
	replace keep = 1 if files2 == ".log"
	drop if keep == 1
	gen script = "erase " + files
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
