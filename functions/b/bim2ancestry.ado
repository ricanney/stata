*! 1.0.3 Richard Anney 16apr2019
* clean up-front functions
* clean down stream
* work in local data folder
* create log file
*! 1.0.2 Richard Anney 21jan2019
* fix folder clean up 
* update to include bim name in output
* update to save files in bim2ancestry file and source
*! 1.0.1 Richard Anney 5dec2018
* updated to generate better pc plots - 4dec2018
* requires a file *.aims containing ancestry informative markers
* add routine to deal with missnps after flipping
* add routine to remove legend and label according to population in graph

program  bim2ancestry
local    function bim2ancestry
local    `function'_version v2
syntax , bim(string asis) ref(string asis) 

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
	ref_path_to_short, ref(`ref')
	local `function'_bim_short ${bim_short}
	local `function'_ref_short ${ref_short}
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
	noi di as text"#########################################################################"
	qui { // check files
		foreach file in bim bed fam { 
			noi checkfile, file(`bim'.`file')
			}
		foreach file in bim bed fam population aims { 
			noi checkfile, file(`ref'.`file')
			}
		checkfile, file(${plink})
		checktabbed
		noi di as text"#########################################################################"
		noi di as text""	
		}
	}
qui { // module 3 - identify overlapping aims snps
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3: extract ancestry informative markers from plink binaries"
	noi di as text"#########################################################################"
	qui { // extract aims SNPs from test and ref
		noi di as text"# > .................. extracting AIM SNPs "
		foreach set in ref bim {
			local already_created ${init_root}/bim2ancestry/``function'_bim_short'/``function'_`set'_short'-aims.bed
			capture confirm file `already_created'
			if !_rc {
				noi checkfile, file(`already_created')
				}
			else {	
				!$plink --bfile ``set''        --extract `ref'.aims --mac 5 --make-bed --out ``function'_`set'_short'-aims
				noi checkfile, file(`already_created')
				}
			}
		}
	qui { // restrict aims SNPs to rymk
		noi di as text"# > .................. limit SNPs to RYMK"
		foreach set in ref bim {
			local already_created ${init_root}/bim2rymk/``function'_`set'_short'-aims/``function'_`set'_short'-aims_rymk.dta
			capture confirm file `already_created'
			if !_rc {
				noi checkfile, file(`already_created')
				}
			else {	
				bim2rymk, bim(${init_root}/bim2ancestry/``function'_bim_short'/``function'_`set'_short'-aims)
				noi checkfile, file(`already_created')
				}	
			}
		}
	qui { // merge and remove incompatible SNPs
		noi di as text"# > .................. merge and remove incompatible SNPs"
		local already_created ${init_root}/bim2ancestry/``function'_bim_short'/``function'_bim_short'-aims.extract
		capture confirm file `already_created'
		if !_rc {
			noi checkfile, file(`already_created')
			}
		else {
			use ${init_root}/bim2rymk/``function'_bim_short'-aims/``function'_bim_short'-aims_rymk.dta, clear
			rename (a1 a2 gt) (bim_a1 bim_a2 bim_gt)
			merge 1:1 snp using ${init_root}/bim2rymk/``function'_ref_short'-aims/``function'_ref_short'-aims_rymk.dta
			rename (a1 a2 gt) (ref_a1 ref_a2 ref_gt)
			keep if _m == 3
			ta ref_gt bim_gt,m
			drop if ref_gt == ""
			drop if bim_gt == "" 
			drop if ref_gt == "R"  &  bim_gt == "M" 
			drop if ref_gt == "R"  &  bim_gt == "K" 
			drop if ref_gt == "Y"  &  bim_gt == "M" 
			drop if ref_gt == "Y"  &  bim_gt == "K" 
			drop if ref_gt == "M"  &  bim_gt == "R" 
			drop if ref_gt == "M"  &  bim_gt == "Y" 
			drop if ref_gt == "K"  &  bim_gt == "R" 
			drop if ref_gt == "K"  &  bim_gt == "Y" 
			ta ref_gt bim_gt,m
			outsheet using ``function'_bim_short'-aims.extract, replace non noq
			noi checkfile, file(`already_created')
			}
		}
	qui { // limit to common SNPs
		noi di as text"# > .................. extracting common SNPs"
		foreach set in ref bim {
			local already_created ${init_root}/bim2ancestry/``function'_bim_short'/``function'_`set'_short'-aims-overlap.bed
			capture confirm file `already_created'
			if !_rc {
				noi checkfile, file(`already_created')
				}
			else {			
				!$plink --bfile ``function'_`set'_short'-aims --make-bed --extract ``function'_bim_short'-aims.extract   --out ``function'_`set'_short'-aims-overlap
				noi checkfile, file(`already_created')
				}
			}
		}
	qui { // flip markers
		noi di as text"# > .................. identifying flip errors"
		local already_created ${init_root}/bim2ancestry/``function'_bim_short'/``function'_bim_short'-aims-merged.bed
		capture confirm file `already_created'
		if !_rc {
			noi di as text"# > .................. merging aims binaries"
			noi checkfile, file(`already_created')
			}
		else {	
			!$plink --bfile ``function'_bim_short'-aims-overlap --bmerge ``function'_ref_short'-aims-overlap --merge-equal-pos --allow-no-sex --make-bed --out ``function'_bim_short'-aims-merged
			capture confirm file ``function'_bim_short'-aims-merged-merge.missnp
			if !_rc {
				noi di as text"# > .................. fixing flip errors"
				!$plink --bfile ``function'_bim_short'-aims-overlap --flip ``function'_bim_short'-aims-merged-merge.missnp --allow-no-sex --make-bed --out ``function'_bim_short'-aims-overlap-flip
				noi di as text"# > .................. merging datasets"
				!$plink --bfile ``function'_bim_short'-aims-overlap-flip --bmerge ``function'_ref_short'-aims-overlap --merge-equal-pos --allow-no-sex --make-bed --out ``function'_bim_short'-aims-merged
				}
			else {
				noi di as text"# > .................. no flip errors / datasets merged"
				}
			noi di as text"# > .................. merging aims binaries"
			noi checkfile, file(`already_created')
			}
		}
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 4 - calculate eigenvectors for `function'-merged file
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 4: calculate eigenvecs for merged files"
	noi di as text"#########################################################################"
	local already_created ${init_root}/bim2eigenvec/``function'_bim_short'-aims-merged/``function'_bim_short'-aims-merged_eigenvec.dta
	capture confirm file `already_created'
	if !_rc {
			noi di as text"# > .................. calculate eigenvecs"
			noi checkfile, file(`already_created')
			}
	else {	
			noi di as text"# > .................. calculate eigenvecs"
			bim2eigenvec, bim(${init_root}/bim2ancestry/``function'_bim_short'/``function'_bim_short'-aims-merged)
			noi checkfile, file(`already_created')
			}
	noi di as text"#########################################################################"
	}
qui { // module 5 - define ancestries
	noi di as text""	
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 5: define ancestries"
	noi di as text"#########################################################################"
	local already_created ${init_root}/bim2ancestry/``function'_bim_short'/``function'_bim_short'-`function'-merged_population.dta
	capture confirm file `already_created'
	if !_rc {
		noi di as text"# > .................. merge with reference ancestry"
		use ``function'_bim_short'-`function'-merged_population.dta, clear
		count if pop == "TEST"
		local `function'_N `r(N)'
		noi di as text"# > .................. individuals in analysis = "as result"``function'_N'" 		
		noi checkfile, file(`already_created')
		}
	else {
		noi di as text"# > .................. merge with reference ancestry"
		import delim using "`ref'.population", clear case(lower) varnames(1)
		replace pop = strupper(pop)
		merge 1:1 fid iid using ${init_root}/bim2eigenvec/``function'_bim_short'-aims-merged/``function'_bim_short'-aims-merged_eigenvec.dta
		replace pop = "TEST" if pop == ""
		count if pop == "TEST"
		local `function'_N `r(N)'
		noi di as text"# > .................. individuals in analysis = "as result"``function'_N'"  
		drop _m
		lab var fid 			"Family Identifier"
		lab var iid				"Individual Identifier"
		lab var pop				"Population"
		foreach num of num 1/10 {
			lab var pc`num'		"Ancestry Principal Component `num'"
			}
		save ``function'_bim_short'-`function'-merged_population.dta, replace
		noi checkfile, file(`already_created')
		}
	noi di as text"#########################################################################"
	}
qui { // module 6 - plot graphs
	noi di as text""	
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 6: plot graphs"
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
		local already_created ${init_root}/bim2ancestry/``function'_bim_short'/``function'_bim_short'-`function'_pca.gph
		capture confirm file `already_created'
		if !_rc {
			noi di as text"# > ......... saved to "as result"``function'_bim_short'-`function'_pca.gph"
			}
		else {
			use 	``function'_bim_short'-`function'-merged_population.dta, clear
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
					|| scatter pc`j' pc`i' if pop == "TEST",  $format2 `test' ///
					legend(off) saving(``function'_bim_short'-`function'_pca`j'-pca`i'.gph, replace) nodraw
					}
				}
			graph combine  ``function'_bim_short'-`function'_pca1-pca2.gph ``function'_bim_short'-`function'_pca1-pca3.gph ``function'_bim_short'-`function'_pca2-pca3.gph , col(4) title("") nodraw saving(``function'_bim_short'-`function'_pca.gph, replace)
			noi di as text"# > ......... saved to "as result"``function'_bim_short'-`function'_pca.gph"
			}
		noi di as text"#########################################################################"
		noi di as text""
		}	
	}
qui { // module 7 - write log file
	file open myfile using "``function'_bim_short'-`function'.log", write replace
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
	file write myfile `"# > ......... N in fam ``function'_N'"' _n
	file write myfile `"# > ........ reference ``function'_ref_short'"' _n
	file write myfile `"# > ..... saved graphs ``function'_bim_short'-`function'_pca.gph"' _n
	file write myfile `"# > ..... saved graphs ``function'_bim_short'-`function'_pca1-pca2.gph"' _n
	file write myfile `"# > ..... saved graphs ``function'_bim_short'-`function'_pca1-pca3.gph"' _n
	file write myfile `"# > ..... saved graphs ``function'_bim_short'-`function'_pca2-pca3.gph"' _n
	file write myfile `"# > ....... saved data ``function'_bim_short'-`function'-merged_population.dta"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > overview of `function'"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > 1. limit to ancestry informative markers"' _n
	file write myfile `"# > 2. limit to RYMK SNPs"' _n
	file write myfile `"# > 3. remove incompatible SNPS (e.g. R/M or Y/K)"' _n
	file write myfile `"# > 4. limit to overlapping markers"' _n
	file write myfile `"# > 5. merge/flip/merge"' _n
	file write myfile `"# > 6. calculate eigenvectors (removing high ld regions)"' _n
	file write myfile `"# > 7. plot against reference populations"' _n
	file write myfile `"#########################################################################"' _n
	file close myfile	
	}
qui { // module 8 - clean folder
	files2dta, dir(`c(pwd)')
	gen keep = .
	replace keep = 1 if file == "``function'_bim_short'-`function'_pca.gph"
	replace keep = 1 if file == "``function'_bim_short'-`function'_pca1-pca2.gph"
	replace keep = 1 if file == "``function'_bim_short'-`function'_pca1-pca3.gph"
	replace keep = 1 if file == "``function'_bim_short'-`function'_pca2-pca3.gph"
	replace keep = 1 if file == "``function'_bim_short'-`function'-merged_population.dta"
	replace keep = 1 if file == "``function'_bim_short'-`function'.log"
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


