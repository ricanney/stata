* update code for windows and check
global init_root H:/software/stata/data

global bim1 			H:/databank/2-genotypes/restricted/neonate-ahe/neonate-ahe-qc-v9
global bim2 			H:/databank/2-genotypes/restricted/controls-combined-wtccc1/controls-combined-wtccc1 
global bim2build  H:/databank/3-packages/bim2build/bim2build.dta
global bim2refid  ${init_root}/create_hrc/HRC.r1-1.GRCh37.wgs.mac5/all-hrc-1.1-chrall-mac5
global generef	  ${init_root}/create_generef/Homo_sapiens.GRCh37.87.gtf/
global recombref  ${init_root}/create_recombref/genetic_map_HapMapII_GRCh37/genetic_map_HapMapII_GRCh37-recombination-rate.dta 
global ref		 	  ${init_root}/create_1000genomes/ftp.1000genomes.ebi.ac.uk/all-1000g-phase3-chrall-mac5-v2
global ldref      ${init_root}/create_1000genomes/ftp.1000genomes.ebi.ac.uk/gbr-1000g-phase3-chrall-mac5-v2
global hrcref     ${init_root}/create_hrc/HRC.r1-1.GRCh37.wgs.mac5/all-hrc-1.1-chrall-mac5
global summaryqc1 ${init_root}/summaryqc/demontis-walters-martin-2017-eur-adhd/demontis-walters-martin-2017-eur-adhd
global summaryqc2 ${init_root}/summaryqc/ripke-neale-corvin-2014-scz/ripke-neale-corvin-2014-scz
global summaryqc3 ${init_root}/summaryqc/grove-ripke-als-2017-eur-asd/grove-ripke-als-2017-eur-asd
global summaryqc4 ${init_root}/summaryqc/anney-klei-pinto-2010-autism-cstr/anney-klei-pinto-2010-autism-cstr

qui { // (+) bim_path_to_short
	discard
	bim_path_to_short, bim(${bim}) 
	noi di "${bim_short}"
	qui { // dependencies
			*"none"
			}
	}
qui { // (+) bim2ancestry
	discard
	bim2ancestry, bim(${bim}) ref(${ref})
	graph use "H:\software\stata\data\bim2ancestry\neonate-ahe-qc-v9\neonate-ahe-qc-v9-bim2ancestry_pca.gph"
	graph export "H:\software\stata\data\bim2ancestry\neonate-ahe-qc-v9\neonate-ahe-qc-v9-bim2ancestry_pca.png", as(png) height(1000) width(3000) replace
	graph export "H:\software\stata\functions\b\images\example-bim2ancestry_pca.png", as(png) height(1000) width(3000) replace
	qui { // dependencies
		*(+)bim_path_to_short
		*(+)bim2rymk
		*(+)bim2eigenvec
		*(+)checkfile
		*(+)checktabbed
		*(+)files2dta
		*(+)ref_path_to_short
		}
	}
qui { // (+) bim2count
	discard
	bim2count, bim(${bim})
	noi di "SNPs = ${bim2count_snp}"
	noi di "SNPs = ${bim2count_ind}"
	qui { // dependencies
		*(+)bim_path_to_short
		*(+)checkfile
		}
	}
qui { // (+) bim2dta
	discard
	noi bim2dta, bim(${bim})
	qui { // dependencies
		*(+)bim_path_to_short
		*(+)bim2count
		*(+)checkfile
		*(+)checkloc_name
		*(+)files2dta
		*(+)recodegenotype
		}
	}
qui { // (+) bim2eigenvec
	discard
	noi bim2eigenvec, bim(${bim}) 
	* adjust = of prinicpal components using pc() flag
	qui { // dependencies
		*(+)bim_path_to_short
		*(+)bim2count
		*(+)bim2ldexclude
		*(+)checkfile
		*(+)checktabbed
		*(+)files2dta
		}
	}
qui { // (+) bim2frq
	discard
	bim2frq, bim(${bim})
	qui { // dependencies
		*(+)bim_path_to_short
		*(+)checkfile
		*(+)checktabbed
		*(+)files2dta
		*(+)recodegenotype
		}
	}
qui { // (+) bim2ldexclude
	discard
	noi bim2ldexclude, bim(${bim}) 
	qui { // dependencies
		*(+)bim_path_to_short
		*(+)checkfile
		*(+)files2dta
		}	
	}
qui { // (+) bim2rymk
	discard
	bim2rymk, bim(${bim})
	qui { // dependencies
		*(+)bim_path_to_short
		*(+)bim2count
		*(+)bim2dta
		*(+)checkfile
		*(+)files2dta
		}
	}
qui { // (+) checkfile
	discard
	checkfile, file(${ref}.bim)
	qui { // dependencies
			*"none"
			}
		}
qui { // (-) checkfolder
	qui { // dependencies
			}
		}
qui { // (+) checkloc_name
	discard
	use H:\software\stata\data\bim2dta\neonate-ahe-qc-v9\neonate-ahe-qc-v9_bim.dta, clear
	drop loc_name
	checkloc_name
	qui { // dependencies
		*(+)recodegenotype
			}
		}
qui { // (+) checktabbed
	discard 
	noi checktabbed
	qui { // dependencies
		*(+)checkfile
		}
	}
qui { // (+) create_generef
	discard
	create_generef, version(87)
	qui { // dependencies
		*(+)checkfile
		}
	}
qui { // (+) create_recombref
	discard
	create_recombref	
	qui { // dependencies
		*(+)checkfile
		}	
	}

qui { // (+) create_1000genomes
	discard
	noi create_1000genomes
	qui { // dependencies
		*(+)bim2dta
		*(+)checkfile
		*(+)checkloc_name
		*(+)files2dta
		*(+)recodegenotypes
		}
	}
qui { // (-) create_hapmap3
}
qui { // (+) create_hrc
	discard
	create_hrc
	qui { // dependencies
		*(+)bim2dta
		*(+)checkfile
		*(+)checkloc_name
		*(+)files2dta
		*(+)recodegenotypes
		}
	}
qui { // (+) fam2dta
	discard
	noi fam2dta, fam(${bim})
	qui { // dependencies
		*(+)bim_path_to_short
		*(+)bim2count
		*(+)checkfile
		*(+)files2dta
		}
	}
qui { // (+) files2dta
		discard
		files2dta, dir(H:\software\stata\functions\f)
		qui { // dependencies
			*"none"
			}
		}
qui { // (+) graphgene
	discard
	graphgene, range(chr20:1000000..3000000) generef(${generef}) label(on) png(yes) biotype(all)
	qui { // dependencies
		*(+)checkfile
		*(+)files2dta
		*(+)gene_plotorder
		*(+)intersect_genes
		}
	}

qui { // (+) graphlocuszoom
	discard
	graphlocuszoom, index(rs4680) generef(${generef}) ldref(${ldref}) recombref(${recombref}) summaryqc(${summaryqc1}) label(on) png(yes) biotype(all) 
	graphlocuszoom, index(rs10942171) generef(${generef}) ldref(${ldref}) recombref(${recombref}) summaryqc(${summaryqc1}) label(on) png(yes) biotype(all) 
	graphlocuszoom, index(rs1233578) generef(${generef}) ldref(${ldref}) recombref(${recombref}) summaryqc(${summaryqc2}) label(on) png(yes) biotype(all) 
	graphlocuszoom, index(rs1233578) generef(${generef}) ldref(${ldref}) recombref(${recombref}) summaryqc(${summaryqc2}) label(on) png(yes) range(chr6:26000000..30000000) biotype(protein_coding) 
	graphlocuszoom, index(rs2224274) generef(${generef}) ldref(${ldref}) recombref(${recombref}) summaryqc(${summaryqc3}) label(on) png(yes)  
	graphlocuszoom,  generef(${generef}) ldref(${ldref}) recombref(${recombref}) summaryqc(${summaryqc4}) label(on) png(yes) range(chr20:14000000..16000000) biotype(all)  
	graphlocuszoom,  generef(${generef}) ldref(${ldref}) recombref(${recombref}) summaryqc(${summaryqc3}) label(on) png(yes) range(chr20:14000000..16000000) biotype(all)  

	discard
	graphlocuszoom,  generef(${generef}) ldref(${ldref}) recombref(${recombref}) summaryqc(${summaryqc3}) png(yes) index(rs2224274)
	graphlocuszoom,  generef(${generef}) ldref(${ldref}) recombref(${recombref}) summaryqc(${summaryqc3}) png(yes) range(chr20:14000000..16000000)
	
	use ${summaryqc3}-summaryqc.dta, clear
	keep chr snp bp p
	graphlocuszoom, generef(${generef}) ldref(${ldref}) recombref(${recombref}) project(PGC-iPSYCH-ASD-GWAS-2019) png(yes) index(rs2224274)
	qui { // dependencies
		*(+)checkfile
		*(+)checktabbed
		*(+)files2dta
		*(+)gene_plotorder
		*(+)graphgene
		*(+)ref_path_to_short
		}
	}
qui { // (-) graphmanhattan
}
qui { // (+) intersect_genes
	discard
	intersect_genes, range(chr20:1000000..3000000) generef(${generef}) 
	clear
	set obs 2
	gen intersect_genes_range = ""
	replace intersect_genes_range = "chr20:1000000..3000000" in 1
	replace intersect_genes_range = "chr20:2000000..4000000" in 2
	save ${init_root}/intersect_genes_range.dta, replace
	discard
	intersect_genes, range_list(${init_root}/intersect_genes_range.dta) generef(${generef}) 	
	qui { // dependencies
		*(+)checkfile
		*(+)ref_path_to_short__*.dta
		}
	}		
qui { // (+) recodegenotype
	discard
	import delim using ${ref}.bim, clear
	rename (v5 v6) (a1 a2)
	recodegenotype
	qui { // dependencies
			*"none"
			}
		}
qui { // (+) ref_path_to_short
	discard
	ref_path_to_short, ref(${ref}) 
	noi di "${ref_short}"
	qui { // dependencies
			*"none"
			}	
	}
qui { // (+) summaryqc
	discard
	qui { // pre-process dataset
		use H:\software\stata\data\bim2dta\alspac-ext-b2661-2-1000g-qc-v9-25yr-subset\alspac-ext-b2661-2-1000g-qc-v9-25yr-subset_bim.dta, clear
		keep snp a1 a2
		rename (a1 a2) (b1 b2)
		merge 1:1 snp using alspac-ext-b2661-2-1000g-qc-v9-25yr-subset-s25BAARS.dta
		gen a2 = b2
		replace a2 = b1 if a1 == a2
		drop b1 b2 _m
		gen n = 3122
		tostring bp, replace
		tostring chr, replace
		gen or = exp(beta)
		gen l95  = or - (1.96 * se)
		gen u95  = or + (1.96 * se)
		drop if chr == "23"
		drop if p == .
		}
	discard
	noi summaryqc, input(alspac-ext-b2661-2-1000g-qc-v9-25yr-subset-s25BAARS) output(cardiff-alspac-genetics-2020-s25BAARS) id_ref(${hrcref}) build_ref(${bim2build})
	qui { // dependencies
		*(+)checkfile
		*(+)files2dta
		*(-)graphmanhattan
		*(+)ref_path_to_short
		*(+)snp2build
		*(+)snp2refid
		}
	}
qui { // (+) snp2build
	use H:\software\stata\data\summaryqc\cardiff-alspac-genetics-2020-s25BAARS\cardiff-alspac-genetics-2020-s25BAARS-summaryqc.dta, clear
	snp2build, ref($bim2build)
	qui { // dependencies
		*(+)checkfile
		}
	}
qui { // (-) snp2refid
	discard
	use H:\software\stata\data\summaryqc\cardiff-alspac-genetics-2020-s25BAARS\cardiff-alspac-genetics-2020-s25BAARS-summaryqc.dta, clear		
	snp2refid, ref($bim2refid)
	qui { // dependencies
		*(+)bim2dta
		*(+)checkfile
		*(+)checkloc_name
		*(+)ref_path_to_short
		}
	}
	x
