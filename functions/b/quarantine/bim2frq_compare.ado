*! 1.0.3 Richard Anney 17apr2019
* clean up-front functions
*! 1.0.2 Richard Anney 24jan2019
* - add folder structure 
* - update for compatability with bim2frq
*! 1.0.1 Richard Anney 4dec2018
* updated to take advantage of bim2frq var routine - 4dec2018
* - allow snps to be excluded by variance around reference frequency - 4dec2018
* - shuffle maf to af according to alphabetic order of alleles (renders a better tw scatter) - 4dec2018

program  bim2frq_compare
local    function bim2frq_compare
local    `function'_version v3
syntax , bim(string asis) ref(string asis) [sd(real 4)]

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
	foreach file in bim bed fam {
		noi checkfile, file(`bim'.`file')
		}
	noi checkfile, file(${init_root}/bim2frq/``function'_ref_short'/``function'_ref_short'_frq.dta)
	noi checkfile, file(${plink})
	checktabbed
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 3 - calculate frq by chromosome	
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 3: calculate allele frequency"
	noi di as text"#########################################################################"
	local already_created ${init_root}/bim2frq/``function'_bim_short'/``function'_bim_short'_frq.dta
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		}
	else {
		noi di as text"# > .................. create "as result "``function'_bim_short'_frq.dta"
		noi bim2frq, bim(`bim')
		noi checkfile, file(`already_created')
		}
	local already_created ${init_root}/bim2frq/``function'_ref_short'/``function'_ref_short'_frq.dta
	capture confirm file `already_created'
	if !_rc {
		noi checkfile, file(`already_created')
		}
	else {
		noi di as text"# > .................. create "as result "``function'_ref_short'_frq.dta"
		noi bim2frq, bim(`ref')
		noi checkfile, file(`already_created')
		}
	noi di as text"#########################################################################"
	noi di as text" "
	}
qui { // module 4 - merge frequency files	
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 4: merge frequency files"
	noi di as text"#########################################################################"	
	qui { // merge and rename variables
		noi di as text"# > .......... import "as result"``function'_bim_short'"
		use ${init_root}/bim2frq/``function'_bim_short'/``function'_bim_short'_frq.dta, clear
		count
		local `function'_bimN `r(N)'
		noi di as text"# > .... SNPs in file "as result"``function'_bimN'"
		for var a1 a2 maf gt var: rename X bim_X
		drop if bim_a1 == ""
		drop if bim_a2 == ""
		destring bim_maf, replace force
		noi di as text"# > ........... merge "as result"``function'_ref_short'"
		merge 1:1 snp using ${init_root}/bim2frq/``function'_ref_short'/``function'_ref_short'_frq.dta
		count if _m != 1
		local `function'_refN `r(N)'
		noi di as text"# > .... SNPs in file "as result"``function'_refN'"
		keep if _m == 3
		count
		local `function'_mergeN `r(N)'
		noi di as text"# > ... SNPs in merge "as result"``function'_mergeN'"
		for var a1 a2 maf gt var: rename X ref_X
		drop if ref_a1 == ""
		drop if ref_a2 == ""
		destring ref_maf, replace force
		drop _m
		}
	qui { // process matched SNPs
		qui { // drop incompatible markers
			noi di as text"# > ................. drop incompatible SNPs"
			gen comp = .
			replace comp = 1 if bim_gt == ref_gt
			replace comp = 1 if bim_gt == "Y"  & ref_gt == "R"
			replace comp = 1 if bim_gt == "R"  & ref_gt == "Y"
			replace comp = 1 if bim_gt == "K"  & ref_gt == "M"
			replace comp = 1 if bim_gt == "M"  & ref_gt == "K"
			replace comp = 1 if bim_gt == "ID" & ref_gt == "ID"
			count if comp == .
			local `function'_incompN `r(N)'
			noi di as text"# > ............ drop "as result"``function'_incompN'"as text" incompatible SNPs"
			drop if comp == .
			drop comp
			}
		qui { // drop ambiguous markers
			noi di as text"# > ................. drop ambiguous SNPs"
			gen ambig = .
			replace ambig = 1 if bim_gt == "S"
			replace ambig = 1 if bim_gt == "W"
			count if ambig == 1
			local `function'_ambigN `r(N)'
			noi di as text"# > ............ drop "as result"``function'_ambigN'"as text" ambiguous SNPs"
			drop if ambig == 1
			drop ambig
			}
		qui { // fix strand - identify markers to flip
			noi di as text"# > ................. match strands"
			gen flip = .
			replace flip = 1 if bim_gt == "Y" & ref_gt == "R"
			replace flip = 1 if bim_gt == "R" & ref_gt == "Y"
			replace flip = 1 if bim_gt == "K" & ref_gt == "M"
			replace flip = 1 if bim_gt == "M" & ref_gt == "K"
			foreach allele in a1 a2 {
				gen `allele' = bim_`allele'
				replace `allele' = "A" if bim_`allele' == "T" & flip == 1
				replace `allele' = "C" if bim_`allele' == "G" & flip == 1
				replace `allele' = "G" if bim_`allele' == "C" & flip == 1
				replace `allele' = "T" if bim_`allele' == "A" & flip == 1
				}
			foreach allele in a1 a2 {
				replace  bim_`allele' = `allele'
				drop `allele'
				}
			drop flip
			}
		qui { // map to common a1
			noi di as text"# > ................. map to common reference allele"
			drop bim_gt ref_gt
			rename (bim_a1 bim_a2) (a1 a2)
			recodegenotype
			rename  (a1 a2 recodegenotype_gt ) (bim_a1 bim_a2 gt)
			gen a1 = ""
			replace a1 = "A" if gt == "R"
			replace a1 = "C" if gt == "Y"
			replace a1 = "A" if gt == "M"
			replace a1 = "G" if gt == "K"
			}
		qui { // map freq to common a1
			noi di as text"# > ................. map to frq to reference allele"
			gen frq = bim_maf
			replace frq = 1-bim_maf if a1 == bim_a2
			gen maf = ref_maf
			replace maf = 1-maf if a1 == ref_a2
			}
		qui { // upper and lower bounds from maf var 
			noi di as text"# > ................. calculate bounds from maf variance"
			gen sd = sqrt(ref_var)
			gen lower = maf - (`sd' * sd)
			gen upper = maf + (`sd' * sd)
			keep snp a1 maf frq lower upper 
			rename frq bim_maf
			}
		qui { // identify snp with frq out of range
			noi di as text"# > ................. define out-of bounds SNPs"
			gen oob = .
			replace oob = 1 if bim_maf > upper
			replace oob = 1 if bim_maf < lower
			count if oob == 1
			local `function'_oobN `r(N)'
			noi di as text"# > ............ flag "as result"``function'_oobN'"as text" out-of-bounds SNPs"
			outsheet snp if oob == 1 using ``function'_bim_short'-``function'_ref_short'_`function'.exclude, non noq replace
			noi di as text"# > ....... saving as "as result"``function'_bim_short'-``function'_ref_short'_`function'.exclude"
			}
		}
	noi di as text"#########################################################################"
	noi di as text" "
	}
qui { // module 5 - plot allele frequencies	
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 5: plot allele frequencies	"
	noi di as text"#########################################################################"		
	noi di as text"# > ................. prepare data to plot"
	count
	local `function'_totN `r(N)'
	qui { // speed up graphing - remove duplicate points
		for var bim_maf maf: replace X = round(X,.001)
		egen x = seq(),by(bim_maf maf oob)
		drop if x!=1 & oob ==.
		drop x
		}
	qui { // plot uncorrected 	
		tw scatter maf bim_maf , mlc(black) mfc(blue) mlw(vvthin) m(o) || scatter maf bim_maf if oob == 1,  legend(off) mlc(black) mfc(yellow) mlw(vvthin) m(o) xtitle("frq (``function'_bim_short')") ytitle("frq (``function'_ref_short')")  saving(``function'_bim_short'-1.gph,replace) nodraw
		}
	qui { // plot corrected 	
		tw scatter maf bim_maf if oob == . , mlc(black) mfc(blue) mlw(vvthin) m(o) xtitle("frq (``function'_bim_short')") ytitle("frq (``function'_ref_short')")  saving(``function'_bim_short'-2.gph,replace) nodraw
		}
	qui { // plot combined
		graph combine ``function'_bim_short'-1.gph ``function'_bim_short'-2.gph, ycommon nodraw title("test data: ``function'_bim_short'""reference data: ``function'_ref_short'")
		graph save ``function'_bim_short'-``function'_ref_short'_`function'.gph, replace
		noi di as text"# > . saving graph as "as result"``function'_bim_short'-``function'_ref_short'_`function'.gph"
		}
	noi di as text"#########################################################################"
	noi di as text" "
	}
qui { // module 4 - write log file
	file open myfile using "``function'_bim_short'-``function'_ref_short'_`function'.log", write replace
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
	file write myfile `"# > ... reference SNPs ``function'_ref_short'"' _n
	file write myfile `"# > ...... SNPs in bim ``function'_bimN'"' _n
	file write myfile `"# > ...... SNPs in ref ``function'_refN'"' _n
	file write myfile `"# > .. SNPs post merge ``function'_mergeN'"' _n
	file write myfile `"# > ............. drop ``function'_incompN' incompatible SNPs"' _n
	file write myfile `"# > ............. drop ``function'_ambigN' ambiguous SNPs"' _n
	file write myfile `"# > ... SNPs to assess ``function'_totN'"' _n
	file write myfile `"# > . out-of-bounds as `sd' sd from the ref mean (estimated from variance of maf in ref)"' _n
	file write myfile `"# > ............. flag ``function'_oobN' out-of-bound SNPs"' _n
	file write myfile `"# > ....... saved data ``function'_bim_short'-``function'_ref_short'_`function'.exclude"' _n
	file write myfile `"# > ...... saved graph ``function'_bim_short'-``function'_ref_short'_`function'.gph"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > overview of `function'"' _n
	file write myfile `"#########################################################################"' _n
	file write myfile `"# > 1. create/use and _frq.dta files from bim and ref"' _n
	file write myfile `"# > 2. merge bim and ref"' _n
	file write myfile `"# > 3. drop incompatible SNPs"' _n
	file write myfile `"# > 4. drop ambiguous SNPs"' _n
	file write myfile `"# > 5. define frequency boundaries from variance in referece data"' _n
	file write myfile `"# > 6. flag out of bounds SNPs"' _n
	file write myfile `"# > 7. plot frequencies"' _n
	file write myfile `"#########################################################################"' _n
	file close myfile	
	}
qui { // module 5 - clean folder
	files2dta, dir(`c(pwd)')
	gen keep = .
	replace keep = 1 if file == "``function'_bim_short'-``function'_ref_short'_`function'.exclude"
	replace keep = 1 if file == "``function'_bim_short'-``function'_ref_short'_`function'.gph"
	replace keep = 1 if file == "``function'_bim_short'-``function'_ref_short'_`function'.log"
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
