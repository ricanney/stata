/*
*program*
 snp2build

*description* 
 a command to check genome build from snp chr bp variables

*syntax*
 snp2build,  ref(-reference-)
 
 -reference- download bim2build.dta from github.com/ricanney
*/

program snp2build
syntax , ref(string asis)
noi di as text" "
noi di as text"#########################################################################"
noi di as text"# snp2build"
noi di as text"#########################################################################"
noi di as text"# Started: $S_DATE $S_TIME"
noi di as text"#########################################################################"
preserve
qui { // 1 - introduction
	noi di as text"# > snp2build .............. checking build against (ref) "as result"`ref'"
	noi checkfile, file(`ref')
	}
qui { // 2 - limit sample size
	drop if _n > 50000
	}
qui { // 3 - measure overlap
	keep snp chr bp
	}
qui { // 4 - merge against reference 
		duplicates drop snp, force
		merge 1:1 snp using `ref'
		keep if _m == 3
		tostring bp, replace
		compress
		}
qui { // 5 - measure overlap with reference
		foreach i in 17 18 19 { 
			gen hg`i'_0 = .
			gen hg`i'_1 = .
			replace hg`i'_0 = 1 if bp == hg`i'_chromStart 
			replace hg`i'_1 = 1 if bp == hg`i'_chromEnd 
			}		
		count
		gen all = `r(N)'
		foreach i in 17 18 19 { 
			foreach j in 0 1 { 
				sum hg`i'_`j'
				gen phg`i'_`j' = r(N) / all
				}
			}
		keep in 1
		keep phg17_0 - phg19_1
		xpose, clear v
		rename v1 percentMatched
		rename _v build
		replace build = "hg17 +0" if build == "phg17_0"
		replace build = "hg17 +1" if build == "phg17_1"
		replace build = "hg18 +0" if build == "phg18_0"
		replace build = "hg18 +1" if build == "phg18_1"
		replace build = "hg19 +0" if build == "phg19_0"
		replace build = "hg19 +1" if build == "phg19_1"
		gsort -p
		gen MostLikely = "+++" in 1
		replace MostLikely = "++" if p > 0.9 & MostLikely == ""
		replace MostLikely = "+" if p > 0.8 & MostLikely == ""
		}
qui { // 6 - report findings 
		noi di as text"# > snp2build ................. build overlap reported in " as result"snp2build.report"
		outsheet using snp2build.report, replace noq	
		graph hbar percentMatched , over(build,sort(percentMatched) lab(labs(large))) title("Percentage Match Genome Build") yline(.9, lcol(red))  
		graph export snp2build.eps, as(eps) replace
		!convert snp2build.eps snp2build.png
		noi di as text"# > snp2build .................. build overlap plotted to " as result"snp2build.png"
		keep in 1 
		tostring per, replace force
		gen a = ""
		replace a = "global snp2build " + build
		outsheet a using _tmp.do, non noq replace
		do _tmp.do
		noi di as text"# > snp2build ....................... build identified as " as result"${snp2build}"
		erase _tmp.do
		erase snp2build.eps
		erase snp2build.png
		erase snp2build.report
		}
noi di as text"#########################################################################"
noi di as text"# Completed: $S_DATE $S_TIME"
noi di as text"#########################################################################"
restore
end;	
