*! 1.0.2 Richard Anney 15apr2019
* clean up-front functions
*! 1.0.1 Richard Anney 6mar2019

program  bim2ancestry_keep
local    function bim2ancestry_keep
local    `function'_version v2
syntax , bim(string asis) ancestry(string asis) ref(string asis) [sd(real 4) pc(real 3)]

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
	noi bim_path_to_short, bim(`bim')
	local `function'_bim_short ${bim_short}
	cd ${init_root}	
	!mkdir `function'
	cd `function'
	noi di as text"# > ..... project name "as result"``function'_bim_short'"
	!mkdir ``function'_bim_short'
	cd ``function'_bim_short'
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 2 - define files
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 2: define files"
	qui { // check files
		checktabbed
		noi di as text"#########################################################################"
		noi di as text""	
		}
	}
qui { // module 3 - identify individuals similar to ancestry
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3: identify individuals similar to `ancestry'"
	noi di as text"#########################################################################"
	local already_created ${init_root}/bim2ancestry/``function'_bim_short'/``function'_bim_short'-bim2ancestry-merged_population.dta
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		use `already_created', clear
		}
	else {
		bim2ancestry, bim(`bim') ref(`ref')
		noi checkfile, file(`already_created')
		}
	count if pop == "TEST" 
	local `function'_Nin `r(N)'
	noi di as text"# > ..................... N in file "as result "`r(N)'"
	gen drop = .
	noi di as text"# > .......... ancestry to match to "as result "`ancestry'"
	noi di as text"# > .. deviation from ancestry mean "as result "< `sd' sd"
	noi di as text"# > ............... PCs to evaluate "as result "`pc'"
	gen pop2 = "`ancestry'"
	replace pop = "EAS" if pop == "CHB" & pop2 == "EAS"
	replace pop = "EAS" if pop == "JPT" & pop2 == "EAS"
	replace pop = "EAS" if pop == "CHS" & pop2 == "EAS"
	replace pop = "EAS" if pop == "CDX" & pop2 == "EAS"
	replace pop = "EAS" if pop == "KHV" & pop2 == "EAS"
	replace pop = "EUR" if pop == "CEU" & pop2 == "EUR"
	replace pop = "EUR" if pop == "TSI" & pop2 == "EUR"
	replace pop = "EUR" if pop == "FIN" & pop2 == "EUR"
	replace pop = "EUR" if pop == "GBR" & pop2 == "EUR"
	replace pop = "EUR" if pop == "IBS" & pop2 == "EUR"
	replace pop = "AFR" if pop == "YRI" & pop2 == "AFR"
	replace pop = "AFR" if pop == "LWK" & pop2 == "AFR"
	replace pop = "AFR" if pop == "GWD" & pop2 == "AFR"
	replace pop = "AFR" if pop == "MSL" & pop2 == "AFR"
	replace pop = "AFR" if pop == "ESN" & pop2 == "AFR"
	replace pop = "AFR" if pop == "ASW" & pop2 == "AFR" 
	replace pop = "AFR" if pop == "ACB" & pop2 == "AFR"
	replace pop = "AMR" if pop == "MXL" & pop2 == "AMR"
	replace pop = "AMR" if pop == "PUR" & pop2 == "AMR"
	replace pop = "AMR" if pop == "CLM" & pop2 == "AMR" 
	replace pop = "AMR" if pop == "PEL" & pop2 == "AMR"
	replace pop = "SAS" if pop == "GIH" & pop2 == "SAS" 
	replace pop = "SAS" if pop == "PJL" & pop2 == "SAS"
	replace pop = "SAS" if pop == "BEB" & pop2 == "SAS"
	replace pop = "SAS" if pop == "STU" & pop2 == "SAS"
	replace pop = "SAS" if pop == "ITU" & pop2 == "SAS"
	drop pop2
	foreach i of num 1 / `pc' { 
		sum pc`i' if pop == "`ancestry'"
		gen pc`i'_mean = `r(mean)'
		gen pc`i'_sd   = `r(sd)'
		gen pc`i'_l = pc`i'_mean - (`sd' * pc`i'_sd)
		gen pc`i'_u = pc`i'_mean + (`sd' * pc`i'_sd)
		replace drop = `i' if pc`i' > pc`i'_u & pop == "TEST" & drop == .
		replace drop = `i' if pc`i' < pc`i'_l & pop == "TEST" & drop == .
		count if drop == `i'
		noi di as text"# > .............. dropped from pc`i' "as result "`r(N)'"
		drop pc`i'_mean - pc`i'_u
		}
	noi di as text"#########################################################################"
	}
qui { // module 4 - plot graphs
	noi di as text""	
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 4: plot graphs"
	noi di as text"#########################################################################"
	qui { // ancestries
		/*
		CHB	Han Chinese in Beijing, China	EAS	
		JPT	Japanese in Tokyo, Japan	EAS	
		CHS	Southern Han Chinese	EAS	
		CDX	Chinese Dai in Xishuangbanna, China	EAS	
		KHV	Kinh in Ho Chi Minh City, Vietnam	EAS	
		CEU	Utah Residents (CEPH) with Northern and Western European Ancestry	EUR	
		TSI	Toscani in Italia	EUR	
		FIN	Finnish in Finland	EUR	
		GBR	British in England and Scotland	EUR	
		IBS	Iberian Population in Spain	EUR	
		YRI	Yoruba in Ibadan, Nigeria	AFR	
		LWK	Luhya in Webuye, Kenya	AFR	
		GWD	Gambian in Western Divisions in the Gambia	AFR	
		MSL	Mende in Sierra Leone	AFR	
		ESN	Esan in Nigeria	AFR	
		ASW	Americans of African Ancestry in SW USA	AFR	
		ACB	African Caribbeans in Barbados	AFR	
		MXL	Mexican Ancestry from Los Angeles USA	AMR	
		PUR	Puerto Ricans from Puerto Rico	AMR	
		CLM	Colombians from Medellin, Colombia	AMR	
		PEL	Peruvians from Lima, Peru	AMR	
		GIH	Gujarati Indian from Houston, Texas	SAS	
		PJL	Punjabi from Lahore, Pakistan	SAS	
		BEB	Bengali from Bangladesh	SAS	
		STU	Sri Lankan Tamil from the UK	SAS	
		ITU	Indian Telugu from the UK	SAS	
		*/
		}
	qui { // define ancestry colours
		local eur  mfcolor("red")
		colorscheme 9, palette(Reds)  
		foreach i of num 1/9 {
			local eur`i' mfcolor("`=r(color`i')'")
			}
		local afr  mfcolor("green")
		colorscheme 9, palette(Greens)  
		foreach i of num 1/9 {
			local afr`i' mfcolor("`=r(color`i')'")
			}
		local sas  mfcolor("purple")
		colorscheme 9, palette(Purples)  
		foreach i of num 1/9 {
			local sas`i' mfcolor("`=r(color`i')'")
			}
		local eas  mfcolor("blue")
		colorscheme 9, palette(Blues)  
		foreach i of num 1/9 {
			local eas`i' mfcolor("`=r(color`i')'")
			}
		local amr  mfcolor("yellow")
		colorscheme 9, palette(YlOrBr)  
		foreach i of num 1/9 {
			local amr`i' mfcolor("`=r(color`i')'")
			}
		local test mfcolor("yellow")
		}
	qui { // plot graph
		gen p2 = pop
		egen x =seq(),by(p2)
		replace p2 = "" if x !=1
		global format1 "msiz(small) msymbol(o) mlc(black) mlw(vvthin) mlabel(p2) mlabposition(1) mlabsize(half_tiny) mlabgap(3)"
		global format2 "msiz(small) msymbol(s) mlc(black) mlw(vvthin)"
		foreach i of num 1/3 { 
			foreach j of num  1/3 { 
				tw scatter pc`j' pc`i' if pop == "CHB" ,  $format1 `eas9' ///
				|| scatter pc`j' pc`i' if pop == "JPT" ,  $format1 `eas8' ///
				|| scatter pc`j' pc`i' if pop == "CHS" ,  $format1 `eas7' ///
				|| scatter pc`j' pc`i' if pop == "CDX" ,  $format1 `eas6' ///
				|| scatter pc`j' pc`i' if pop == "KHV" ,  $format1 `eas5' ///
				|| scatter pc`j' pc`i' if pop == "CEU" ,  $format1 `eur9' ///
				|| scatter pc`j' pc`i' if pop == "TSI" ,  $format1 `eur8' ///
				|| scatter pc`j' pc`i' if pop == "FIN" ,  $format1 `eur7' ///
				|| scatter pc`j' pc`i' if pop == "GBR" ,  $format1 `eur6' ///
				|| scatter pc`j' pc`i' if pop == "IBS" ,  $format1 `eur5' ///
				|| scatter pc`j' pc`i' if pop == "YRI" ,  $format1 `afr9' ///
				|| scatter pc`j' pc`i' if pop == "LWK" ,  $format1 `afr8' ///
				|| scatter pc`j' pc`i' if pop == "GWD" ,  $format1 `afr7' ///
				|| scatter pc`j' pc`i' if pop == "MSL" ,  $format1 `afr6' ///
				|| scatter pc`j' pc`i' if pop == "ESN" ,  $format1 `afr5' ///
				|| scatter pc`j' pc`i' if pop == "ASW" ,  $format1 `afr4' ///
				|| scatter pc`j' pc`i' if pop == "ACB" ,  $format1 `afr3' ///
				|| scatter pc`j' pc`i' if pop == "MXL" ,  $format1 `amr9' ///
				|| scatter pc`j' pc`i' if pop == "PUR" ,  $format1 `amr8' ///
				|| scatter pc`j' pc`i' if pop == "CLM" ,  $format1 `amr7' ///
				|| scatter pc`j' pc`i' if pop == "PEL" ,  $format1 `amr6' ///
				|| scatter pc`j' pc`i' if pop == "GIH" ,  $format1 `sas9' /// 
				|| scatter pc`j' pc`i' if pop == "PJL" ,  $format1 `sas8' /// 
				|| scatter pc`j' pc`i' if pop == "BEB" ,  $format1 `sas7' /// 
				|| scatter pc`j' pc`i' if pop == "STU" ,  $format1 `sas6' /// 
				|| scatter pc`j' pc`i' if pop == "ITU" ,  $format1 `sas5' /// 
				|| scatter pc`j' pc`i' if pop == "TEST" & drop == .,  $format2 `test' ///
				legend(off) saving(``function'_bim_short'-`function'-`ancestry'-keep_pca`j'-pca`i'.gph, replace) nodraw
				}
			}
		graph combine  ``function'_bim_short'-`function'-`ancestry'-keep_pca1-pca2.gph ``function'_bim_short'-`function'-`ancestry'-keep_pca1-pca3.gph ``function'_bim_short'-`function'-`ancestry'-keep_pca2-pca3.gph , col(4) title("") nodraw saving(``function'_bim_short'-`function'-`ancestry'-keep_pca.gph, replace)
		noi di as text"# > ...................... saved to "as result"``function'_bim_short'-`function'-`ancestry'-keep_pca.gph"
		}
	qui { // define keep
		keep if pop == "TEST"
		keep if drop == .
		count if pop == "TEST" 
		local `function'_Nout `r(N)'

		noi di as text"# > ............. N in cleaned file "as result "`r(N)'"
		noi di as text"# > ...................... saved to "as result"``function'_bim_short'-`function'.keep"
		outsheet fid iid using ``function'_bim_short'-`function'-`ancestry'.keep, non noq replace
		}
	}
qui { // module 5 - write log file
	file open myfile using "``function'_bim_short'-`function'-`ancestry'-keep.log", write replace
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
	file write myfile `"# > ..... project name ``function'_bim_short'"' _n
	file write myfile `"# > .... plink version ${plink}"' _n
	file write myfile `"# > ......... N before ``function'_Nin'"' _n
	file write myfile `"# > .......... N after ``function'_Nout'"' _n
	file write myfile `"# > ......... limit to `ancestry'"' _n
	file write myfile `"# > ........... within `sd' sd of pc mean"' _n
	file write myfile `"# > ..... evaluate top `pc' ancestry PCs"' _n
	file write myfile `"# > ..... saved graphs ``function'_bim_short'-`function'-`ancestry'-keep_pca.gph"' _n
	file write myfile `"# > ..... saved graphs ``function'_bim_short'-`function'-`ancestry'-keep_pca1-pca2.gph"' _n
	file write myfile `"# > ..... saved graphs ``function'_bim_short'-`function'-`ancestry'-keep_pca1-pca3.gph"' _n
	file write myfile `"# > ..... saved graphs ``function'_bim_short'-`function'-`ancestry'-keep_pca2-pca3.gph"' _n
	file write myfile `"# > ....... saved data ``function'_bim_short'-`function'-`ancestry'.keep"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > overview of `function'"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > 1. import bim2ancestry-merged_population.dta "' _n
	file write myfile `"# > 2. define mean and sd of each pc for ancestry group"' _n
	file write myfile `"# > 3. define test individual who are within sd range for the defined number of pcs"' _n
	file write myfile `"# > 4. keep test individuals that meet similarity"' _n
	file write myfile `"# > 5. plot against reference populations"' _n
	file write myfile `"#########################################################################"' _n
	file close myfile	
	}
qui { // module 8 - clean folder
	files2dta, dir(`c(pwd)')
	gen keep = .
	split files,p("kee")
	replace keep = 1 if files3 == "p_pca.gph"
	replace keep = 1 if files3 == "p_pca1-pca2.gph"
	replace keep = 1 if files3 == "p_pca1-pca3.gph"
	replace keep = 1 if files3 == "p_pca2-pca3.gph"
	replace keep = 1 if files3 == "p.log"
	replace keep = 1 if files3 == "p"
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


