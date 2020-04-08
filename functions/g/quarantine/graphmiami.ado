*! 1.0.2 Richard Anney 15jan2019
* - remove need for reference SNPs
*! 1.0.3 Richard Anney 30apr2019
* - update upfront functions

program  graphmiami
local    function graphmiami
local    `function'_version v4
syntax , summaryqc1(string asis) summaryqc2(string asis) range(string asis) generef(string asis) ldref(string asis) label(string asis) [png(string asis) biotype(string asis)]

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
	ref_path_to_short, ref(`summaryqc1')
	local `function'_summaryqc1_short ${ref_short}
	ref_path_to_short, ref(`summaryqc2')
	local `function'_summaryqc2_short ${ref_short}
	local `function'_generef_short Homo_sapiens.GRCh37.87.gtf_exon
	ref_path_to_short, ref(`ldref')
	local `function'_ldref_short ${ref_short}
	local project ``function'_summaryqc1_short'-``function'_summaryqc2_short'
	cd ${init_root}	
	!mkdir `function'
	cd `function'
	noi di as text"# > ..... project name "as result"`project'"
	!mkdir `project'
	cd `project'
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 2 - define files
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 2: define files"
	noi di as text"#########################################################################"
	noi di as text"# > ......... locating files to process"
	foreach i in bim bed fam {
		noi checkfile, file(`ldref'.`i')
		}
	noi checkfile, file(`summaryqc1'-summaryqc.dta)
	noi checkfile, file(`summaryqc2'-summaryqc.dta)
	noi checkfile, file(${plink})
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 3 - define region to plot
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3: define region to plot"
	noi di as text"#########################################################################"
	noi di as text"# > .......... GWAS #1 "as result"``function'_summaryqc1_short'"
	noi di as text"# > .......... GWAS #2 "as result"``function'_summaryqc2_short'"
	noi di as text"# > ......... plotting "
	clear
	set obs 1
	gen `function'_range = "`range'"
	split `function'_range,p("chr"":""..")
	local chr 	= `function'_range2[1]
	local from 	= `function'_range3[1]
	local to 		= `function'_range4[1]
	outsheet `function'_range2 `function'_range3 `function'_range4  using `function'-chr`chr'_`from'_`to'.bed, non noq replace
	noi di as text"# > .............. chr "as result"`chr'"
	noi di as text"# > .......... from bp "as result"`from'"
	noi di as text"# > ............ to bp "as result"`to'"
	noi di as text"#########################################################################"
	noi di as text""
	}
qui { // module 4 - plot genes
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 4: plot genes"
	noi di as text"#########################################################################"
	local already_created ${init_root}/graphgene/``function'_generef_short'/``function'_generef_short'-chr`chr'_`from'_`to'_graphgene.dta
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		}
	else {	
		graphgene, range(`range') generef(`generef') label(`label') biotype(all)
		noi checkfile, file(`already_created')
		}
	noi di as text"#########################################################################"
	noi di as text""
	}
qui { // module 5 - combine summaryqc
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 5: combine summaryqc"
	noi di as text"#########################################################################"
	local already_created `project'_chr`chr'_`from'_`to'_`function'.dta
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		use `already_created', clear
		outsheet snp using `project'_chr`chr'_`from'_`to'_`function'.extract, non noq replace
		}
	else{
		foreach i in 1 2 {
			local already_created ``function'_summaryqc`i'_short'_chr`chr'_`from'_`to'_`function'.dta
			capture confirm file `already_created'
			if !_rc {
				noi checkfile, file(`already_created')
				use `already_created', clear
				}
			else {
				use `summaryqc`i''-summaryqc.dta, clear
				keep snp chr bp p
				for var chr bp p: destring X,replace
				keep if chr == `chr'
				drop if bp < `from'
				drop if bp > `to'
				gen log10p_`i' = round(-log10(p),.01)
				keep snp chr bp log10p_`i'
				duplicates drop
				save ``function'_summaryqc`i'_short'_chr`chr'_`from'_`to'_`function'.dta,replace
				noi checkfile, file(`already_created')
				}
			}
		merge 1:1 snp using ``function'_summaryqc1_short'_chr`chr'_`from'_`to'_`function'.dta
		erase ``function'_summaryqc1_short'_chr`chr'_`from'_`to'_`function'.dta
		erase ``function'_summaryqc2_short'_chr`chr'_`from'_`to'_`function'.dta
		drop _m
		save `project'_chr`chr'_`from'_`to'_`function'.dta, replace
		noi checkfile, file(`project'_chr`chr'_`from'_`to'_`function'.dta)
		outsheet snp using `project'_chr`chr'_`from'_`to'_`function'.extract, non noq replace
		}
	noi di as text"#########################################################################"
	noi di as text""
	}
qui { // module 6 - calculate LD in reference
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 6: calculate LD"
	noi di as text"#########################################################################"
	foreach i in 1 2 {
		gsort -log10p_`i'
		gen index = "global index`i' " in 1
		replace index = index + snp in 1
		outsheet index in 1 using temp.do, non noq replace
		do temp.do
		erase temp.do
		noi di as text"# > .. define index as "as result"${index`i'}"
		drop index
		}
	noi di as text"# > .. define LD using "as result"``function'_ldref_short'"
	local already_created `project'_`function'-chr`chr'_`from'_`to'_r2.dta
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		use `already_created', clear
		}
	else {
		if chr != 23 {
			!$plink --bfile `ldref' --extract `project'_chr`chr'_`from'_`to'_`function'.extract --allow-no-sex --r2             --ld-window-kb 100000 --ld-window 99999 --ld-window-r2 0 --ld-snp ${index1} --out `project'_chr`chr'_`from'_`to'_`function'_1
			!$plink --bfile `ldref' --extract `project'_chr`chr'_`from'_`to'_`function'.extract --allow-no-sex --r2             --ld-window-kb 100000 --ld-window 99999 --ld-window-r2 0 --ld-snp ${index2} --out `project'_chr`chr'_`from'_`to'_`function'_2
			}
		else {
			!$plink --bfile `ldref' --extract `project'_chr`chr'_`from'_`to'_`function'.extract --allow-no-sex --r2 --ld-xchr 1 --ld-window-kb 100000 --ld-window 99999 --ld-window-r2 0 --ld-snp ${index1}  --out `project'_chr`chr'_`from'_`to'_`function'_1
			!$plink --bfile `ldref' --extract `project'_chr`chr'_`from'_`to'_`function'.extract --allow-no-sex --r2 --ld-xchr 1 --ld-window-kb 100000 --ld-window 99999 --ld-window-r2 0 --ld-snp ${index2}  --out `project'_chr`chr'_`from'_`to'_`function'_2
			}
		!$tabbed `project'_chr`chr'_`from'_`to'_`function'_1.ld
		import delim using `project'_chr`chr'_`from'_`to'_`function'_1.ld.tabbed, varnames(1) clear
		keep snp_b r2
		rename (snp_b r2) (snp rsquare1)
		merge 1:1 snp  using `project'_chr`chr'_`from'_`to'_`function'.dta
		keep snp chr bp log10p_1 log10p_2 rsquare1
		save `project'_`function'-chr`chr'_`from'_`to'_r2.dta,replace
		!$tabbed `project'_chr`chr'_`from'_`to'_`function'_2.ld
		import delim using `project'_chr`chr'_`from'_`to'_`function'_2.ld.tabbed, varnames(1) clear
		keep snp_b r2
		rename (snp_b r2) (snp rsquare2)
		merge 1:1 snp  using `project'_`function'-chr`chr'_`from'_`to'_r2.dta
		keep snp chr bp log10p_1 log10p_2 rsquare1 rsquare2
		save `project'_`function'-chr`chr'_`from'_`to'_r2.dta,replace		
		checkfile, file(`project'_`function'-chr`chr'_`from'_`to'_r2.dta)
		}
	noi di as text"#########################################################################"
	noi di as text""
	}
qui { // module 7 - preparing data for plot
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 7: preparing data for plot"
	noi di as text"#########################################################################"
	noi { // join gene data
		append using ${init_root}/graphgene/``function'_generef_short'/``function'_generef_short'-chr`chr'_`from'_`to'_graphgene.dta
		drop label_plot biotype_plot order
		qui { // define plotting criteria
			gen label_plot = "`label'"
			if label_plot == "off" {
				replace label = ""
				noi di as text"# > . plot gene labels "as result"no"
				}
			else {
				noi di as text"# > . plot gene labels "as result"yes"
				}
			gen biotype_plot = "`biotype'"
			if biotype_plot == "all" {
				noi di as text"# > .... plot biotypes "as result"all"
				}
			else if biotype_plot == "protein_coding" {
				noi di as text"# > .... plot biotypes "as result"protein_coding only"
				drop if biotype != "protein_coding" & start != .
				}
			else {
				noi di as text"# > .... plot biotypes "as result"protein_coding only (default : set to all)"
				drop if biotype != "protein_coding" & start != .
				}
			}
		qui { // define plotting order
			gene_plotorder, from(`from') to(`to')
			}
		qui { // define colors
			colorscheme 9, palette(Reds) 
			foreach color of num 1/9{
				local gwas1_color`color' "`r(color`color')'"
				}
			colorscheme 9, palette(Blues) 
			foreach color of num 1/9{
				local gwas2_color`color' "`r(color`color')'"
				}
			}
		qui { // update main gwas (limit to p < 1e-15)
			gen log10p_1_plot = log10p_1
			for var log10p_1_plot : replace X = 15  if X > 15 & X != .
			}		
		qui { // update reflection gwas (limit to p < 1e-15)
			sum order
			gen y = `r(max)'
			gen log10p_2_plot = log10p_2
			for var log10p_2_plot : replace X = -X 
			for var log10p_2_plot : replace X = -15  if X < -15 & X != .
			for var log10p_2_plot : replace X = X -y - 4
			}
		qui { // define order for gene plots
			replace order = -order - 2
			sort bp
			}
		qui { // define ylab
			sum log10p_2_plot if log10p_2 == 0
			gen ylab_0  = `r(min)'
			gen ylab_5  =  ylab_0 - 5
			gen ylab_10 =  ylab_0 - 10
			gen ylab_15 =  ylab_0 - 15
			gen yline   =  ylab_0 - 7.2
			for var ylab_0 - yline : tostring X,replace
			gen ylab = `"ylab("' + ylab_15 + `""15" "' + ylab_10 + `""10" "' + ylab_5 + `""5" "' + ylab_0 + `""0" 0"0" 5"5" 10"10" 15"15")"'
			local ylab 	= ylab[1]
			gen yline2 = "yline(7.2) yline(" + yline + ")"
			local yline 	= yline2[1]
			}
		qui { // define titles	
			gen dx_1 = "``function'_summaryqc1_short'" in 1 
			gen dx_2 = "``function'_summaryqc2_short'" in 1 
			gen dx_1p = 15
			gen dx_2p = -15 - y - 8
			}
		qui { // flag LD
			foreach i in 1 2 {
				gen `function'_p_`i'r0_1   = log10p_`i'_plot if inrange(rsquare`i',0.0,0.1)
				gen `function'_p_`i'r1_2   = log10p_`i'_plot if inrange(rsquare`i',0.1,0.2)
				gen `function'_p_`i'r2_3   = log10p_`i'_plot if inrange(rsquare`i',0.2,0.3)
				gen `function'_p_`i'r3_4   = log10p_`i'_plot if inrange(rsquare`i',0.3,0.4)
				gen `function'_p_`i'r4_5   = log10p_`i'_plot if inrange(rsquare`i',0.4,0.5)
				gen `function'_p_`i'r5_6   = log10p_`i'_plot if inrange(rsquare`i',0.5,0.6)
				gen `function'_p_`i'r6_7   = log10p_`i'_plot if inrange(rsquare`i',0.6,0.7)
				gen `function'_p_`i'r7_8   = log10p_`i'_plot if inrange(rsquare`i',0.7,0.8)
				gen `function'_p_`i'r8_9   = log10p_`i'_plot if inrange(rsquare`i',0.8,0.9)
				gen `function'_p_`i'r9_1   = log10p_`i'_plot if inrange(rsquare`i',0.9,1.0)
				}
			}
		qui { // define graph common commands
			local `function'_log10p_plot sort msymbol(O) mlwidth(vthin) 
			}
		qui { // plot graph
			#delimit ;
			tw scatter  `function'_p_1r0_1 bp, ``function'_log10p_plot' mcolor("`gwas1_color3'") msize(small)
			|| scatter  `function'_p_1r1_2 bp, ``function'_log10p_plot' mcolor("`gwas1_color3'") msize(small)
			|| scatter  `function'_p_1r2_3 bp, ``function'_log10p_plot' mcolor("`gwas1_color4'") msize(small)
			|| scatter  `function'_p_1r3_4 bp, ``function'_log10p_plot' mcolor("`gwas1_color4'") msize(small)
			|| scatter  `function'_p_1r4_5 bp, ``function'_log10p_plot' mcolor("`gwas1_color5'") msize(small)
			|| scatter  `function'_p_1r5_6 bp, ``function'_log10p_plot' mcolor("`gwas1_color6'") msize(small)
			|| scatter  `function'_p_1r6_7 bp, ``function'_log10p_plot' mcolor("`gwas1_color7'") msize(small)
			|| scatter  `function'_p_1r7_8 bp, ``function'_log10p_plot' mcolor("`gwas1_color8'") msize(small)
			|| scatter  `function'_p_1r8_9 bp, ``function'_log10p_plot' mcolor("`gwas1_color9'") msize(small)
			|| scatter  `function'_p_1r9_1 bp, ``function'_log10p_plot' mcolor("`gwas1_color9'") msize(small)
			|| scatter  `function'_p_2r0_1 bp, ``function'_log10p_plot' mcolor("`gwas2_color3'") msize(small)
			|| scatter  `function'_p_2r1_2 bp, ``function'_log10p_plot' mcolor("`gwas2_color3'") msize(small)
			|| scatter  `function'_p_2r2_3 bp, ``function'_log10p_plot' mcolor("`gwas2_color4'") msize(small)
			|| scatter  `function'_p_2r3_4 bp, ``function'_log10p_plot' mcolor("`gwas2_color4'") msize(small)
			|| scatter  `function'_p_2r4_5 bp, ``function'_log10p_plot' mcolor("`gwas2_color5'") msize(small)
			|| scatter  `function'_p_2r5_6 bp, ``function'_log10p_plot' mcolor("`gwas2_color6'") msize(small)
			|| scatter  `function'_p_2r6_7 bp, ``function'_log10p_plot' mcolor("`gwas2_color7'") msize(small)
			|| scatter  `function'_p_2r7_8 bp, ``function'_log10p_plot' mcolor("`gwas2_color8'") msize(small)
			|| scatter  `function'_p_2r8_9 bp, ``function'_log10p_plot' mcolor("`gwas2_color9'") msize(small)
			|| scatter  `function'_p_2r9_1 bp, ``function'_log10p_plot' mcolor("`gwas2_color9'") msize(small)
			|| scatter dx_1p bp, m(none) mlabel(dx_1) mlabpos(3) mlabcolor(black)   
			|| scatter dx_2p bp, m(none) mlabel(dx_2) mlabpos(3) mlabcolor(black)   
			|| rspike  start end 		order if class == "protein_coding", hor lcolor(green) lwidth(vvthin) 
			|| rspike  txstart_ txend_ 	order if class == "protein_coding", hor lcolor(green) lwidth(${gene_plotorder_e}) 
			|| rspike  start end 		order if class == "other(1)", hor lcolor(red) lwidth(vvthin) 
			|| rspike  txstart_ txend_ 	order if class == "other(1)", hor lcolor(red) lwidth(${gene_plotorder_e}) 
			|| rspike  start end 		order if class == "other(2)", hor lcolor(blue) lwidth(vvthin) 
			|| rspike  txstart_ txend_ 	order if class == "other(2)", hor lcolor(blue) lwidth(${gene_plotorder_e}) 						
			|| scatter order start if pos == 11  , msymbol(i) mlabel(label) mlabpos(9) mlabcolor(black) mlabsize(${gene_plotorder_f}) 
			|| scatter order end   if pos == 1   , msymbol(i) mlabel(label) mlabpos(3) mlabcolor(black) mlabsize(${gene_plotorder_f}) 
			`ylab'
			xlab(,format(%10.0f))
			legend(off) ytitle("-log10(P)") 
			xtitle("Chromosome `chr'") 
			`yline' nodraw
			saving(`project'_chr`chr'_`from'_`to'_`function'.gph, replace)
			;
			#delimit cr
			noi di as text"# > .... plot graph to "as result"`project'_chr`chr'_`from'_`to'_`function'.gph"
			}
		qui { // plot png 
			clear 
			set obs 1
			gen png = "`png'"
			if png == "yes" {
				graph use `project'_chr`chr'_`from'_`to'_`function'.gph
				noi di as text"# > ..... saved png to "as result"`project'_chr`chr'_`from'_`to'_`function'.png"
				graph export `project'_chr`chr'_`from'_`to'_`function'.png, as(png) height(2000) width(3000) replace
				window manage close graph
				}
			else {
				}
			}
		}
	noi di as text"#########################################################################"
	noi di as text""
	}
qui { // module 8 - write log file
	file open myfile using "`project'_chr`chr'_`from'_`to'_`function'.log", write replace
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
	file write myfile `"# > .......... gwas #1 ``function'_summaryqc1_short'"' _n
	file write myfile `"# > .......... gwas #2 ``function'_summaryqc2_short'"' _n
	file write myfile `"# >  co-ordinates file ``function'_generef_short'"' _n
	file write myfile `"# > ..... ld reference ``function'_ldref_short'"' _n
	file write myfile `"# > ..... genome build hg19"' _n
	file write myfile `"# > .............. chr chr`chr'"' _n
	file write myfile `"# > ............. from `from'"' _n
	file write myfile `"# > ............... to `to'"' _n
	file write myfile `"# > ..... saved graphs `project'_chr`chr'_`from'_`to'_`function'.gph"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > overview of `function'"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > 1. parse region from gwas summary files "' _n
	file write myfile `"# > 2. merge parsed files "' _n
	file write myfile `"# > 3. plot genes within region "' _n
	file write myfile `"# > 4. calculate LD in region "' _n
	file write myfile `"# > 5. convert gwas#2 to mirror "' _n
	file write myfile `"# > 6. limit max-p to 1e-15 "' _n
	file write myfile `"# > 7. plot region "' _n
	file write myfile `"#########################################################################"' _n
	file close myfile	
	}	
qui { // module 9 - clean folder
	files2dta, dir(${init_root}/`function'/`project'/)
	split files, p("_`function'")
	drop if files2 == ".log"
	drop if files2 == ".gph"
	drop if files2 == ".dta"
	drop if files2 == ".png"
	gen script = "erase " + files
	outsheet script using temp.do, non noq replace 
	do temp.do
	erase temp.do
	clear
	}
qui cd `return'
restore
qui di as text"#########################################################################"
qui di as text"# Completed: $S_DATE $S_TIME"
qui di as text"#########################################################################"
end;
	
