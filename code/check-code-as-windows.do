* update code for windows and check

global bim1 			H:/databank/2-genotypes/restricted/neonate-ahe/neonate-ahe-qc-v9
global bim2 			H:/databank/2-genotypes/restricted/controls-combined-wtccc1/controls-combined-wtccc1 
global generef	  H:/software/stata/data/create_generef/Homo_sapiens.GRCh37.87.gtf/
global recombref  H:/software/stata/data/create_recombref/genetic_map_HapMapII_GRCh37/genetic_map_HapMapII_GRCh37-recombination-rate.dta 
global ref		 	  H:/software/stata/data/create_1000genomes/ftp.1000genomes.ebi.ac.uk/all-1000g-phase3-chrall-mac5-v2
global ldref      H:/software/stata/data/create_1000genomes/ftp.1000genomes.ebi.ac.uk/
global summaryqc1 H:/software/stata/data/summaryqc/demontis-walters-martin-2017-eur-adhd/demontis-walters-martin-2017-eur-adhd
global summaryqc2 H:/software/stata/data/summaryqc/ripke-neale-corvin-2014-scz/ripke-neale-corvin-2014-scz
global summaryqc3 H:/software/stata/data/summaryqc/grove-ripke-als-2017-eur-asd/grove-ripke-als-2017-eur-asd
global summaryqc4 H:/software/stata/data/summaryqc/anney-klei-pinto-2010-autism-cstr/anney-klei-pinto-2010-autism-cstr

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
*qui { // (-) create_hapmap3
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

qui { // (-) graphlocuszoom
	discard
	graphlocuszoom, index(rs4680) generef(${generef}) ldref(${ldref}) recombref(${recombref}) summaryqc(${summaryqc1}) label(on) png(yes) biotype(all) 
	graphlocuszoom, index(rs10942171) generef(${generef}) ldref(${ldref}) recombref(${recombref}) summaryqc(${summaryqc1}) label(on) png(yes) biotype(all) 
	graphlocuszoom, index(rs1233578) generef(${generef}) ldref(${ldref}) recombref(${recombref}) summaryqc(${summaryqc2}) label(on) png(yes) biotype(all) 
	graphlocuszoom, index(rs1233578) generef(${generef}) ldref(${ldref}) recombref(${recombref}) summaryqc(${summaryqc2}) label(on) png(yes) range(chr6:26000000..30000000) biotype(protein_coding) 
	graphlocuszoom, index(rs2224274) generef(${generef}) ldref(${ldref}) recombref(${recombref}) summaryqc(${summaryqc3}) label(on) png(yes)  
	graphlocuszoom,  generef(${generef}) ldref(${ldref}) recombref(${recombref}) summaryqc(${summaryqc4}) label(on) png(yes) range(chr20:14000000..16000000) biotype(all)  
	graphlocuszoom,  generef(${generef}) ldref(${ldref}) recombref(${recombref}) summaryqc(${summaryqc3}) label(on) png(yes) range(chr20:14000000..16000000) biotype(all)  

	discard
	use ${summaryqc3}-summaryqc.dta, clear
	keep chr snp bp p
	graphlocuszoom, index(rs2224274)                generef(${generef}) ldref(${ldref}) recombref(${recombref}) project(PGC-iPSYCH ASD GWAS 2019)




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
qui { // (-) summaryqc
	