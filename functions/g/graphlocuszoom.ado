*! 1.0.2 Richard Anney 30apr2019
* clean up-front functions
* change output name to include chr and loc

program  graphlocuszoom
local    function graphlocuszoom
local    `function'_version v3
syntax , generef(string asis) ldref(string asis) recombref(string asis) [index(string asis) summaryqc(string asis) project(string asis) gwsp(real 7.3) range(string asis) png(string asis) label(string asis) biotype(string asis)]

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
	local `function'_generef_short Homo_sapiens.GRCh37.87.gtf_exon
	ref_path_to_short, ref(`recombref')
	local `function'_recombref_short ${ref_short}
	ref_path_to_short, ref(`ldref')
	local `function'_ldref_short ${ref_short}
	cd ${init_root}	
	!mkdir `function'
	cd `function'
	noi di as text""
	}
qui { // module 2 - define files
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 2: define files"
	noi di as text"#########################################################################"
	foreach i in bim bed fam {
		noi checkfile, file(`ldref'.`i')
		}
	noi checkfile, file(`recombref')
	noi checkfile, file($plink)
	checktabbed
	noi di as text"#########################################################################"
	noi di as text""
	}
qui { // module 3 - define project data	
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3: define project data	"
	noi di as text"#########################################################################"
	noi di as text"# > define project name"
	gen `function'_summaryqc = "`summaryqc'"
	if `function'_summaryqc != "" {
		drop `function'_summaryqc
		ref_path_to_short, ref(`summaryqc')
		noi di as text"# > .................. "as result"summaryqc specified"
		noi di as text"# > .................. "as result"project name based on summaryqc"
		local project ${ref_short}
		}
	else {
		noi di as text"# > .................. "as result"summaryqc not specified"
		noi di as text"# > .................. "as result"current data will be used - checking variables"
		qui { // check variables
			keep chr bp snp p
			foreach i in chr bp p {
				capture confirm variable `i'
				if !_rc {
					noi di as text"# > .................. "as result"`i' "as text"is present"
					}
				else {
					noi di as text"# > .................. "as error"`i' "as text"is absent"
					exit
					}
				capture confirm variable numeric `i'
				if !_rc {
					noi di as text"# > .................. "as result"`i' "as text"is numeric"
					}
				else {
					noi di as text"# > .......... convert "as result"`i' "as text"to numeric"
					destring `i', replace
					}
				}
			foreach i in snp {
				capture confirm variable `i'
				if !_rc {
					noi di as text"# > .................. "as result"`i' "as text"is present"
					}
				else {
					noi di as text"# > .................. "as error"`i' "as text"is absent"
					exit
					}
				}
			}
		gen `function'_project = "`project'"
		if `function'_project == "" {
			noi di as text"# > .................. "as result"project not specified"
			local project unknown
			}
		else {
			noi di as text"# > ..... project name "as result"`project'"
			}
		}	
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 4 - create folders
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 4: create folders	"
	noi di as text"#########################################################################"
	noi di as text"# > ..... project name "as result"`project'"
	!mkdir `project'
	cd `project'
	gen `function'_summaryqc = "`summaryqc'"
	if `function'_summaryqc == "" {		
		keep chr bp snp p 
		tostring chr, replace
		save `project'-summaryqc.dta, replace
		clear
		set obs 1
		gen `function'_summaryqc = "${init_root}/`function'/`project'/`project'"
		local summaryqc = `function'_summaryqc[1]
		}
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 5 - define region
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 5: define region"
	noi di as text"#########################################################################"
	clear 
	set obs 1
	gen `function'_range = "`range'"
	if `function'_range != ""   { // set range based on user input
		noi di as text"# > . range defined as "as result"`range'"
		split `function'_range,p("chr"":""..")
		local chr 	= `function'_range2[1]
		local from 	= `function'_range3[1]
		local to 		= `function'_range4[1]
		noi di as text"# > .............. chr "as result"chr`chr'"
		noi di as text"# > .......... from bp "as result"`from'"
		noi di as text"# > ............ to bp "as result"`to'"
		gen `function'_index = "`index'"
		if `function'_index != "" { // set index based on user input
			use if chr =="`chr'" using `summaryqc'-summaryqc.dta, clear
			destring bp, replace
			drop if bp < `from'
			drop if bp > `to'
			gen found = 0
			replace found = 1 if snp == "`index'"
			sum found
			if `r(max)' == 1        { // user input index is present
				noi di as text"# > ........ index SNP "as result"`index'"
				}
			else                    { // user input index is absent - use maxP alternative
				noi di as text"# > ........ index SNP "as result"`index' - absent"
				sort p
				local index = snp[1]
				noi di as text"# > ... maxP index SNP "as result"`index'"
				}
			}
		else                      { // index is absent - use maxP alternative
			use if chr =="`chr'" using `summaryqc'-summaryqc.dta, clear
			destring bp, replace
			drop if bp < `from'
			drop if bp > `to'
			sort p
			local index = snp[1]
			noi di as text"# > ... maxP index SNP "as result"`index'"
			}
		}
	else 	                      { // set range based index
		gen `function'_index = "`index'"
		if `function'_index != "" { // set index based on user input
			use snp chr bp if snp =="`index'" using `summaryqc'-summaryqc.dta, clear
			count
			if `r(N)' == 0          { // user input index is absent
				noi di as text"# > .................. "as error"`index' not found in data"
				exit
				}
			else                    { // user input index present
				noi di as text"# > ........ index SNP "as result"`index'"
				destring bp, replace
				replace bp = round(bp, 100000)
				gen from = bp - 500000
				gen to   = bp + 500000
				local chr 	= chr[1]
				local from 	= from[1]
				local to 		= to[1]
				local range chr`chr':`from'..`to'
				noi di as text"# > . range defined as "as result"`range'"
				noi di as text"# > .............. chr "as result"chr`chr'"
				noi di as text"# > .......... from bp "as result"`from'"
				noi di as text"# > ............ to bp "as result"`to'"				
				}
			}
		}
	noi di as text"#########################################################################"
	noi di as text""
	}
qui { // module 6 - define genes in region
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 6: define genes in region"
	noi di as text"#########################################################################"
	local already_created ${init_root}/graphgene/``function'_generef_short'/``function'_generef_short'-chr`chr'_`from'_`to'_graphgene.dta
	capture confirm file `already_created'
	if !_rc {
		local home `c(pwd)'
		cd ${init_root}/graphgene/``function'_generef_short'
		noi checkfile, file(``function'_generef_short'-chr`chr'_`from'_`to'_graphgene.dta)
		cd `home'
		}
	else {
		qui graphgene, range(`range') generef(`generef') label(`label') biotype(`biotype')
		local home `c(pwd)'
		cd ${init_root}/graphgene/``function'_generef_short'
		noi checkfile, file(``function'_generef_short'-chr`chr'_`from'_`to'_graphgene.dta)
		cd `home'
		}
	noi di as text"#########################################################################"
	noi di as text""
	}
qui { // module 7 - define recombination rate in region
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 7: define recombination rate in region"
	noi di as text"#########################################################################"
	local already_created `project'-chr`chr'_`from'_`to'-recombref.dta
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		}
	else {
		use `recombref', clear
		keep if chr == `chr'
		drop if bp < `from'
		drop if bp > `to'
		expand 3
		egen x = seq(),by(chr bp)
		replace bp = bp-1 if x == 1
		replace bp = bp+1 if x == 3
		replace recomb_rate = 0 if x != 2
		drop x
		gsort chr bp -recomb
		egen x = seq(),by(chr bp)
		keep if x == 1
		drop x
		save `project'-chr`chr'_`from'_`to'-recombref.dta, replace
		noi checkfile, file(`already_created')
		}
	noi di as text"#########################################################################"
	noi di as text""
	}	
qui { // module 8 - define LD for index in region 
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 8: define LD for `index' in region "
	noi di as text"#########################################################################"
	local already_created `project'-chr`chr'_`from'_`to'-ldref-`index'.dta
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		}
	else {
		use chr bp snp if chr == "`chr'" using `summaryqc'-summaryqc.dta, replace
		destring bp, replace
		drop if bp < `from'
		drop if bp > `to'
		outsheet snp using `project'-chr`chr'_`from'_`to'_`function'.extract, non noq replace
		if chr != "23" { // define ld in autosomes
			!$plink --bfile `ldref' ///
							--extract `project'-chr`chr'_`from'_`to'_`function'.extract ///
							--allow-no-sex ///
							--r2 --ld-window-kb 100000 --ld-window 99999 --ld-window-r2 0 ///
							--ld-snp `index' ///
							--out `project'-chr`chr'_`from'_`to'_`function'
			}
		else           { // define ld in x-chromosome
			!$plink --bfile `ldref' ///
							--extract `project'-chr`chr'_`from'_`to'_`function'.extract ///
							--allow-no-sex ///
							--r2 --ld-xchr 1 --ld-window-kb 100000 --ld-window 99999 --ld-window-r2 0 ///
							--ld-snp `index' ///
							--out `project'-chr`chr'_`from'_`to'_`function'
				}
		!$tabbed `project'-chr`chr'_`from'_`to'_`function'.ld
		import delim using `project'-chr`chr'_`from'_`to'_`function'.ld.tabbed, varnames(1) clear
		keep snp_b r2
		rename (snp_b r2) (snp rsquare)
		save `project'-chr`chr'_`from'_`to'-ldref-`index'.dta ,replace	
		noi checkfile, file(`already_created')
		}
	noi di as text"#########################################################################"
	noi di as text""
	}
qui { // module 9 - combine datasets
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 9: combine datasets"
	noi di as text"#########################################################################"
	local already_created `project'-chr`chr'_`from'_`to'-`index'-combined.dta
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		}
	else {
		use snp chr bp p if chr =="`chr'" using `summaryqc'-summaryqc.dta, clear
		destring chr, replace
		destring bp, replace
		drop if bp < `from'
		drop if bp > `to'
		merge 1:1 snp using `project'-chr`chr'_`from'_`to'-ldref-`index'.dta
		drop _m
		merge m:1 chr bp using `project'-chr`chr'_`from'_`to'-recombref.dta
		drop _m
		append using ${init_root}/graphgene/``function'_generef_short'/``function'_generef_short'-chr`chr'_`from'_`to'_graphgene.dta
		keep  bp snp p rsquare recomb_rate start end txstart_ txend_ class pos label biotype
		order bp snp p rsquare recomb_rate start end txstart_ txend_ class pos label biotype
		save `project'-chr`chr'_`from'_`to'-`index'-combined.dta ,replace	
		noi checkfile, file(`already_created')
		}
	noi di as text"#########################################################################"
	noi di as text""
	}
qui { // module 10 - prepare data for plotting
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 10: prepare data for plotting"
	noi di as text"#########################################################################"
	use `project'-chr`chr'_`from'_`to'-`index'-combined.dta, clear
	qui { // define -log10(p) by ld
		gen log10p = -log10(p)
		gen log10p_index = log10p if snp == "`index'"
		gen log10p0 = log10p if inrange(rsquare,0,.1)
		gen log10p1 = log10p if inrange(rsquare,.1,.2)
		gen log10p2 = log10p if inrange(rsquare,.2,.3)
		gen log10p3 = log10p if inrange(rsquare,.3,.4)
		gen log10p4 = log10p if inrange(rsquare,.4,.5)
		gen log10p5 = log10p if inrange(rsquare,.5,.6)
		gen log10p6 = log10p if inrange(rsquare,.6,.7)
		gen log10p7 = log10p if inrange(rsquare,.7,.8)
		gen log10p8 = log10p if inrange(rsquare,.8,.9)
		gen log10p9 = log10p if inrange(rsquare,.9,1)
		}
	qui { // define gene plotting criteria
		gen label_plot = "`label'"
		if label_plot[1] == "off" {
			replace label = ""
			noi di as text"# > . plot gene labels "as result"no"
			}
		else {
			noi di as text"# > . plot gene labels "as result"yes"
			}
		gen biotype_plot = "`biotype'"
		if biotype_plot[1] == "all" {
			noi di as text"# > .... plot biotypes "as result"all"
			}
		else if biotype_plot[1] == "protein_coding" {
			noi di as text"# > .... plot biotypes "as result"protein_coding only"
			drop if biotype != "protein_coding" & start != .
			}
		else  {
			noi di as text"# > .... plot biotypes "as result"all"
			}
		}
	qui { // define plotting order
		gene_plotorder, from(`from') to(`to')
		}
	qui { // create legend variable
		foreach num of num 0/9 {
			gen legend_p`num' = . 
			}
		}
	qui { // crop gene transcripts
		replace start = `from' if start <= `from' 
		replace end = `to' if end >= `to' & end ! = .
		replace txstart_ = `from' if txstart_ <= `from' 
		replace txend_ = `to' if txend_ >= `to' & txend_ ! = .
		}
	qui { // define color scheme
		colorpalette YlOrRd, nograph
			foreach color of num 1/9 {
				local color`color' "`r(p`color')'"
				}		
		*	colorscheme 9, palette(YlOrRd) 
		*	foreach color of num 1/9{
		*		local color`color' "`r(color`color')'"
		*		}
			local color0 gs10
			}
	qui { // create ylabel1 ytitle1
			gen ymin = 0
			sum log10p
			gen ymax = round(`r(max)' + 5,10)
			gen ygap = round(ymax / 5,1)
			local ymin = ymin[1]
			local ymax = ymax[1]
			local ygap = ygap[1]
			local ylabel1 ylabel(`ymin'(`ygap')`ymax', nogrid)
			local ytitle1	ytitle("-log10 P-value", size(small) justification(right))
			}
	qui { // create xlabel xtitle
			gen xgap = round((`to'-`from') /2, 100000)
			local xgap = xgap[1]
			local xlabel xlabel(`from'(`xgap')`to', labs(small) format(%10.0f) nogrid)
			local xtitle xtitle(" ""Chromosome `chr' (hg19)", size(small)) 
			}
	qui { // define rendering - legend plot 
			foreach num of num 1/9 {
				local legend_plot`num' sort msymbol(S) msize(small) mlwidth(vthin) mlc(gs0) mcolor("`color`num''") yaxis(1)
				}
			local legend_plot0 sort msymbol(S) msize(small) mlwidth(vthin) mlc(gs0) mcolor(`color0') yaxis(1)  
			}
	qui { // define rendering - index_snp
			local index_marker	sort mcolor("107 174 214") msymbol(D) mlwidth(vthin) msize(medium) mlc(gs0) yaxis(1)
			local index_point		sort mcolor(black)         msymbol(o) mlwidth(vthin) msize(small)  mlc(gs0) mlabel(snp) mlabpos(12) mlabcolor(black) mlabsize(vsmall) yaxis(1)  
			}
	qui { // define rendering - ld_snp
			foreach num of num 1/9 {
				local ld_marker`num' sort msymbol(O) mlwidth(vthin) mlc(gs0)  mcolor("`color`num''")  msize(vsmall) yaxis(1)
				}
			local ld_marker0       sort msymbol(O) mlwidth(vthin) mlc(gs10) mcolor(`color0')      msize(vsmall) yaxis(1)
			}
	qui { // update order to render underneath panel
			replace order = -order
			replace order =  order -2	
			sum order
			gen dummy_yaxis2 = .
			replace dummy_yaxis2 = `r(min)' in 1
			replace dummy_yaxis = `r(max)' in 2
			replace dummy_yaxis = ymax in 3
			}
	qui { // create ylabel2 ytitle2
			gen recomb_rate_adj = recomb_rate / ymax 
			gen recomb_label = _n - 1
			replace recomb_label =round((100/ymax)*recomb_label,1)
			gen N = _n -1
			foreach num in 0 20 40 60 80 100 {
				gen recomb_label`num' = .
				replace recomb_label`num' = N if recomb_label == `num'
				sum recomb_label`num'
				local recomb_label`num' = `r(max)'
				*noi di "`recomb_label`num''"
				drop recomb_label`num'
				}
			local ylabel2 ylabel(`recomb_label0'"0" `recomb_label20'"20" `recomb_label40'"40" `recomb_label60'"60" `recomb_label80'"80" `recomb_label100'"100", axis(2))   
			*noi di`"`ylabel2'"'
			local ytitle2	ytitle("Recombination Rate (cM/Mb)", size(small) axis(2) )
			}
	noi di as text"#########################################################################"
	noi di as text""
	}
qui { // module 11 - plot graph
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 11: plot graph"
	noi di as text"#########################################################################"
	#delimit ;
	tw scatter legend_p0    	bp		, `legend_plot0'
	||scatter  legend_p1    	bp		, `legend_plot1'
	||scatter  legend_p2    	bp		, `legend_plot2'
	||scatter  legend_p3    	bp		, `legend_plot3'
	||scatter  legend_p4    	bp		, `legend_plot4'
	||scatter  legend_p5    	bp		, `legend_plot5'
	||scatter  legend_p6    	bp		, `legend_plot6'
	||scatter  legend_p7    	bp		, `legend_plot7'
	||scatter  legend_p8    	bp		, `legend_plot8'
	||scatter  legend_p9    	bp		, `legend_plot9'
	||scatter  log10p0      	bp		, `ld_marker0'     
	||scatter  log10p1      	bp		, `ld_marker1'     
	||scatter  log10p2      	bp		, `ld_marker2'     
	||scatter  log10p3      	bp		, `ld_marker3'     
	||scatter  log10p4      	bp		, `ld_marker4'     
	||scatter  log10p5      	bp		, `ld_marker5'     
	||scatter  log10p6      	bp		, `ld_marker6'     
	||scatter  log10p7      	bp		, `ld_marker7'     
	||scatter  log10p8      	bp		, `ld_marker8'     
	||scatter  log10p9      	bp		, `ld_marker9'     
	||scatter  log10p_index 	bp 		, `index_marker' 
	||scatter  log10p_index 	bp 		, `index_point'
	||rspike  start         	end 		order if class == "protein_coding", hor lcolor(green) lwidth(vvthin) 
	||rspike  txstart_      	txend_ 	order if class == "protein_coding", hor lcolor(green) lwidth(${gene_plotorder_e}) 
	||rspike  start         	end 		order if class == "other(1)", 			hor lcolor(red) 	lwidth(vvthin) 
	||rspike  txstart_      	txend_  order if class == "other(1)", 			hor lcolor(red) 	lwidth(${gene_plotorder_e}) 
	||rspike  start        	 	end 		order if class == "other(2)", 			hor lcolor(blue) 	lwidth(vvthin) 
	||rspike  txstart_     	 	txend_ 	order if class == "other(2)", 			hor lcolor(blue) 	lwidth(${gene_plotorder_e}) 						
	||scatter order         	start 	if pos == 11  , msymbol(i) mlabel(label) mlabpos(9) mlabcolor(black) mlabsize(${gene_plotorder_f}) 
	||scatter order         	end   	if pos == 1   , msymbol(i) mlabel(label) mlabpos(3) mlabcolor(black) mlabsize(${gene_plotorder_f}) 
	||line    recomb_rate_adj bp		, sort lpattern(solid) lwidth(vthin) lcolor(blue) 
	||scatter dummy_yaxis   	bp		,  msymbol(i) yaxis(2) 
	`ylabel1'
	`ytitle1'
	`ylabel2'
	`ytitle2'
	`xlabel'
	`xtitle'
	legend(region(lc(black)) order(10 "0.9" 9 "0.8" 8 "0.7" 7 "0.6" 6 "0.5" 5 "0.4" 4 "0.3" 3 "0.2" 2 "0.1" 1 "0.0") size(vsmall) rowgap(zero) symp(3) textfirst ring(0) bm(tiny) pos(11) row(10) subtitle("rsquare",size(vsmall)))
	yline(-1, lpattern(solid) lwidth(vthin) lcolor(black))
	yline(5, lpattern(dash) lwidth(vthin) lcolor(orange))
	yline(7.2, lpattern(dash) lwidth(vthin) lcolor(red))
	title("`project'")
	ysize(5) xsize(10)
	fysize(150) fxsize(200)
	nodraw
	saving(`project'-chr`chr'_`from'_`to'_`function'.gph, replace)
	;
	#delimit cr
	noi di as text"# > ..... saved gph to "as result"`project'-chr`chr'_`from'_`to'_`function'.gph"
	clear 
	set obs 1
	gen png = "`png'"
	if png == "yes" {
		noi di as text"# > ..... saved png to "as result"`project'-chr`chr'_`from'_`to'_`function'.png"
		graph use `project'-chr`chr'_`from'_`to'_`function'.gph
		graph export  `project'-chr`chr'_`from'_`to'_`function'.png, as(png) height(2000) width(6000) replace
		window manage close graph
		}
	else {
		}
	noi di as text"#########################################################################"	
	noi di as text""
	}
qui { // module 10 - write log file
	file open myfile using "`project'-chr`chr'_`from'_`to'_`function'.log", write replace
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
	file write myfile `"# > ..... project name `project'"' _n
	file write myfile `"# > ..... ld reference ``function'_ldref_short'"' _n
	file write myfile `"# >  recombination ref ``function'_recombref_short'"' _n
	file write myfile `"# > ............ index `index'"' _n
	file write myfile `"# > ..... genome build hg19"' _n
	file write myfile `"# > .............. chr chr`chr'"' _n
	file write myfile `"# > ............. from `from'"' _n
	file write myfile `"# > ............... to `to'"' _n
	file write myfile `"# > ..... saved graphs `project'-chr`chr'_`from'_`to'_`function'.gph"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > overview of `function'"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > 1. identify index snp in data"' _n
	file write myfile `"# > 2. determine chromosome and region to plot"' _n
	file write myfile `"# > 3. define gene transcripts in region using graphgene"' _n
	file write myfile `"# > 4. calculate ld between index and snps in region using reference"' _n
	file write myfile `"# > 5. parse recombination rate for region"' _n
	file write myfile `"# > 6. combine data"' _n
	file write myfile `"# > 7. plot combined data"' _n
	file write myfile `"#########################################################################"' _n
	file close myfile	
	}
qui { // module 12 - clean folder
	files2dta, dir(${init_root}/`function'/`project'/)
	keep files
	split files, p("_`function'")
	drop if files2 == ".log"
	drop if files2 == ".gph"
	drop if files2 == ".png"
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

   
