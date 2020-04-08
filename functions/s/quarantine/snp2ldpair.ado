*! 1.0.1 Richard Anney 7feb2020

program  snp2ldpair
local    function snp2ldpair
local    `function'_version v1

syntax ,  ldref(string asis) [snp(string asis) snp_list(string asis)]

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
	ref_path_to_short, ref(`ldref')
	local `function'_ldref_short ${ref_short}
	cd ${init_root}	
	!mkdir `function'
	cd `function'
	!mkdir ``function'_ldref_short'
	cd ``function'_ldref_short'
	noi di as text""	
	}
qui { // module 2 - define files
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 2: define files"
	noi di as text"#########################################################################"
	noi checkfile, file(${plink})
	noi di as text"# > ......... locating files to process"
	qui { // cheack ldref
		foreach file in bim bed fam {
			noi checkfile, file(`ldref'.`file')
			}
		}	
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 3 - define LD partners
	noi di as text"#########################################################################"
	noi di as text"# SECTION - 4: define LD partners for INDEX SNPs"
	noi di as text"#########################################################################"
	noi di as text"# > define LD pairs at r^2 > 0.8 in ``function'_ldref_short'"
	clear
	set obs 1
	gen SNP = "`snp'"
	if SNP != "" {
		noi di as text"# > ....... SNP to map "as result"`snp'"
		local already_created ${init_root}/`function'/``function'_ldref_short'/``function'_ldref_short'-`snp'-`function'.dta
		capture confirm file `already_created'
		if !_rc {
			noi checkfile, file(`already_created')
			}
		else {
			outsheet using ``function'_ldref_short'-`snp'.snp-in, non noq replace
			noi di as text"# > mapping LD-pairs using plink"
			!plink2 --bfile `ldref' --r2 in-phase --ld-window 100000 --ld-window-kb 100000 --ld-snp-list ``function'_ldref_short'-`snp'.snp-in --ld-window-r2 0.8 --out ``function'_ldref_short'-`snp'.snp-out
			!$tabbed ``function'_ldref_short'-`snp'.snp-out.ld
			import delim using ``function'_ldref_short'-`snp'.snp-out.ld.tabbed, clear varnames(1)
			rename (snp_a snp_b) (SNP PROXY_SNP)
			for var chr_b bp_b: tostring X, replace
			gen PROXY_POSITION = "chr" + chr_b + ":" + bp_b
			rename r2 R2
			gen SNP_A1 = substr(phase,1,1)
			gen PROXY_A1 = substr(phase,2,1)
			gen SNP_A2 = substr(phase,4,1)
			gen PROXY_A2 = substr(phase,5,1)
			keep  SNP SNP_A1 SNP_A2 PROXY_SNP PROXY_A1 PROXY_A2 PROXY_POSITION R2
			order SNP SNP_A1 SNP_A2 PROXY_SNP PROXY_A1 PROXY_A2 PROXY_POSITION R2
			gen a1 = SNP_A1
			gen a2 = SNP_A2
			recodegenotype
			rename recode SNP_GT
			drop a1 a2
			save ``function'_ldref_short'-`snp'-`function'.dta, replace
			file open myfile using "``function'_ldref_short'-`snp'-`function'.log", write replace
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
			file write myfile `"# > ....... reference ``function'_ldref_short'"' _n
			file write myfile `"# > ...... saved data ``function'_ldref_short'-`snp'-`function'.dta"' _n
			file write myfile `"#########################################################################"' _n
			file write myfile `"# > overview of `function'"' _n
			file write myfile `"#########################################################################"' _n
			file write myfile `"# > 1. identify SNPs"' _n
			file write myfile `"# > 2. extend LD to capture proxy_SNPs"' _n
			file write myfile `"#########################################################################"' _n
			file close myfile	
			}
		}
	else {
		clear
		set obs 1
		gen SNP_LIST = "`snp_list'"
		if SNP_LIST != "" {
			ref_path_to_short, ref(`snp_list')
			local `function'_sl_short ${ref_short}
			noi di as text"# > .. SNP list to map "as result"`snp_list'"
			local already_created ${init_root}/`function'/``function'_ldref_short'/``function'_ldref_short'-``function'_snp_list_short'-`function'.dta
			capture confirm file `already_created'
			if !_rc {
				noi checkfile, file(`already_created')
				}
			else {
				noi di as text"# > mapping LD-pairs using plink"
				!$plink  --bfile `ldref' --r2 in-phase --ld-window 100000 --ld-window-kb 100000 --ld-snp-list `snp_list' --ld-window-r2 0.8 --out ``function'_ldref_short'-``function'_sl_short'.snp-out
				!$tabbed ``function'_ldref_short'-``function'_sl_short'.snp-out.ld
				import delim using ``function'_ldref_short'-``function'_sl_short'.snp-out.ld.tabbed, clear varnames(1)
				rename (snp_a snp_b) (SNP PROXY_SNP)
				for var chr_b bp_b: tostring X, replace
				gen PROXY_POSITION = "chr" + chr_b + ":" + bp_b
				rename r2 R2
				gen SNP_A1 = substr(phase,1,1)
				gen PROXY_A1 = substr(phase,2,1)
				gen SNP_A2 = substr(phase,4,1)
				gen PROXY_A2 = substr(phase,5,1)
				keep  SNP SNP_A1 SNP_A2 PROXY_SNP PROXY_A1 PROXY_A2 PROXY_POSITION R2
				order SNP SNP_A1 SNP_A2 PROXY_SNP PROXY_A1 PROXY_A2 PROXY_POSITION R2
				gen a1 = SNP_A1
				gen a2 = SNP_A2
				recodegenotype
				rename recode SNP_GT
				drop a1 a2
				save ``function'_ldref_short'-``function'_sl_short'-`function'.dta, replace
				file open myfile using "``function'_ldref_short'-``function'_sl_short'-`function'.log", write replace
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
				file write myfile `"# > ....... reference ``function'_ldref_short'"' _n
				file write myfile `"# > ...... saved data ``function'_ldref_short'-``function'_sl_short'-`function'.dta"' _n
				file write myfile `"#########################################################################"' _n
				file write myfile `"# > overview of `function'"' _n
				file write myfile `"#########################################################################"' _n
				file write myfile `"# > 1. identify SNPs"' _n
				file write myfile `"# > 2. extend LD to capture proxy_SNPs"' _n
				file write myfile `"#########################################################################"' _n
				file close myfile	
				}
			}	
		}		
	noi di as text"#########################################################################"
	noi di as text""	
	}
qui { // module 5 - clean folder
	files2dta, dir(`c(pwd)')
	gen drop = .
	split files, p("ldpair")
	drop if files2 == ".log"
	drop if files2 == ".dta"
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

   
