program genotypeqc2report_subroutine
syntax
end
mata:
dh = _docx_new()
// set default document font, size, color, and orientation
_docx_set_color(dh, "000000")
_docx_set_font(dh, "Consolas")
_docx_set_size(dh, 18)

// TITLE PAGE
_docx_paragraph_new(dh, "")
_docx_paragraph_add_text(dh, "# GENOTYPE ARRAY QUALITY CONTROL REPORT")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# =================================================================================")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# This report was created via the function -genotypeqc2report- ")  
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# AUTHOR ..................... Richard Anney")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# REPORT BUGS TO ............. anneyr@cardiff.ac.uk")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# SOURCE AVAILABLE FROM ...... https://github.com/ricanney")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# LICENSE .................... CC BY-NC (CC Attribution-NonCommercial)")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh,`"# ............................ all functions considered beta "caveat emptor""')
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# =================================================================================")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# STARTED .................... $S_DATE $S_TIME")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# USERNAME ................... `c(username)'")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# OPERATING SYSTEM ........... `c(os)'")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# =================================================================================")
_docx_paragraph_add_linebreak(dh) 
_docx_paragraph_add_text(dh, "# > INPUT BINARIES ........... ${genotypeqc2report_infile}")
_docx_paragraph_add_linebreak(dh) 
_docx_paragraph_add_text(dh, "# > INPUT BINARIES ARRAY ..... ${genotypeqc2report_array}; JACCARD = ${genotypeqc2report_jaccard}")
_docx_paragraph_add_linebreak(dh) 
_docx_paragraph_add_text(dh, "# > INPUT BINARIES BUILD ..... ${genotypeqc2report_inbuild}")
_docx_paragraph_add_linebreak(dh) 
_docx_paragraph_add_text(dh, "# > INPUT BINARIES N(SNP) .... ${genotypeqc2report_insnp}")
_docx_paragraph_add_linebreak(dh) 
_docx_paragraph_add_text(dh, "# > INPUT BINARIES N(IND) .... ${genotypeqc2report_inind}")
_docx_paragraph_add_linebreak(dh) 
_docx_paragraph_add_text(dh, "# =================================================================================")
_docx_paragraph_add_linebreak(dh) 
_docx_paragraph_add_text(dh, "# > OUTPUT BINARIES .......... ${genotypeqc2report_outfile}")
_docx_paragraph_add_linebreak(dh) 
_docx_paragraph_add_text(dh, "# > OUTPUT BINARIES BUILD .... ${genotypeqc2report_outbuild}")
_docx_paragraph_add_linebreak(dh) 
_docx_paragraph_add_text(dh, "# > OUTPUT BINARIES N(SNP) ... ${genotypeqc2report_outsnp}")
_docx_paragraph_add_linebreak(dh) 
_docx_paragraph_add_text(dh, "# > OUTPUT BINARIES N(IND) ... ${genotypeqc2report_outind}")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# =================================================================================")
_docx_paragraph_add_linebreak(dh) 
_docx_paragraph_add_text(dh, "# > REFERENCE FOR NAMING ..... ${genotypeqc2report_refid}")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# > REFERENCE FOR FREQUENCY .. ${genotypeqc2report_frqcompare}")
_docx_paragraph_add_linebreak(dh) 
_docx_paragraph_add_text(dh, "# > REFERENCE FOR ANCESTRY ... ${genotypeqc2report_ancestry}")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# =================================================================================")
_docx_paragraph_add_linebreak(dh) 
_docx_paragraph_add_text(dh, "# > ROUNDS OF QC ............. ${genotypeqc2report_rounds}")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# > --mac (threshold) ........ ${genotypeqc2report_mac}")
_docx_paragraph_add_linebreak(dh) 
_docx_paragraph_add_text(dh, "# > --geno (threshold) ....... ${genotypeqc2report_geno}")
_docx_paragraph_add_linebreak(dh) 
_docx_paragraph_add_text(dh, "# > --mind (threshold) ....... ${genotypeqc2report_mind}")
_docx_paragraph_add_linebreak(dh) 
_docx_paragraph_add_text(dh, "# > --hardy (threshold) ...... ${genotypeqc2report_hardy}")
_docx_paragraph_add_linebreak(dh) 
_docx_paragraph_add_text(dh, "# > --het (threshold) ........ ${genotypeqc2report_het}")
_docx_paragraph_add_linebreak(dh) 
_docx_paragraph_add_text(dh, "# > kinship .................. dup = ${genotypeqc2report_kind}; 1st = ${genotypeqc2report_kinf}; 2nd = ${genotypeqc2report_kins}; 3rd = ${genotypeqc2report_kint}")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# =================================================================================")

_docx_add_pagebreak(dh,)
_docx_paragraph_add_text(dh, "# =================================================================================")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# GENOME BUILD")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# =================================================================================")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# > INPUT BINARIES BUILD ..... ${genotypeqc2report_inbuild}")
_docx_paragraph_add_linebreak(dh) 
_docx_paragraph_add_text(dh, "# > OUTPUT BINARIES BUILD .... ${genotypeqc2report_outbuild}")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# =================================================================================")
_docx_paragraph_add_linebreak(dh) 
_docx_paragraph_add_text(dh, "This genome build of the input array is assessed using an rsid comparison using the")
_docx_paragraph_add_linebreak(dh) 
_docx_paragraph_add_text(dh, "function -bim2build-. As part of the -genotypeqc- routine all output binaries are  ")
_docx_paragraph_add_linebreak(dh) 
_docx_paragraph_add_text(dh, "converted to hg19+1.")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# =================================================================================")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# =================================================================================")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# GENOTYPE ARRAY")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# =================================================================================")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# > INPUT BINARIES ARRAY ..... ${genotypeqc2report_array}; JACCARD = ${genotypeqc2report_jaccard}")
_docx_paragraph_add_linebreak(dh) 
_docx_paragraph_add_text(dh, "# =================================================================================")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "The genotyping array platform is assessed using a snpid comparison using the       ")
_docx_paragraph_add_linebreak(dh) 
_docx_paragraph_add_text(dh, "function -bim2array-. Using snpid (over rsid) allows identification of non rsid     ")
_docx_paragraph_add_linebreak(dh) 
_docx_paragraph_add_text(dh, "arrays such as affy. The program assesses the genotype-array against a panel of    ")
_docx_paragraph_add_linebreak(dh) 
_docx_paragraph_add_text(dh, "known arrays, most of which were sourced and adapted from the resources made       ")
_docx_paragraph_add_linebreak(dh) 
_docx_paragraph_add_text(dh, "available at (http://www.well.ox.ac.uk/~wrayner/strand/). Importantly, these are   ")
_docx_paragraph_add_linebreak(dh) 
_docx_paragraph_add_text(dh, "only the best estimates based on overlap coefficients.")
_docx_paragraph_add_text(dh, "# =================================================================================")
_docx_paragraph_add_linebreak(dh)
_docx_add_pagebreak(dh,)
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# =================================================================================")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# COMPARISON OF ALLELE FREQUENCIES")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# =================================================================================")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# > REFERENCE FOR FREQUENCY .. ${genotypeqc2report_frqcompare}")
_docx_paragraph_add_linebreak(dh) 
_docx_paragraph_add_text(dh, "# =================================================================================")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "The -genotypeqc- pipeline uses the function -bim2frq_compare- to compare allele    ")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "frequencies against a reference set of genotypes. We use SNP concordance to prune  ")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "markers with deviation from the reference. SNPs with > 4 sd from the expected      ")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "frequency are excluded from the dataset. Deviation in allele frequency can result  ")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "from incorrect ancestry matching as well errors in marker annotation.")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# =================================================================================")
_docx_paragraph_add_linebreak(dh)
_docx_image_add(dh,"${init_root}/genotypeqc2report/${genotypeqc2report_infile}/${genotypeqc2report_outfile}_bim2frq_compare.png")

_docx_add_pagebreak(dh,)
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# =================================================================================")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# COMPARISON OF ANCESTRIES")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# =================================================================================")
_docx_paragraph_add_linebreak(dh) 
_docx_paragraph_add_text(dh, "# > REFERENCE FOR ANCESTRY ... ${genotypeqc2report_ancestry}")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# =================================================================================")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "Ancestry comparison is performed against a a reference panel using the function    ")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "-bim2ancestry-. These plots are created from the *.eigenvec files created using the")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "--pca flag in plink.                                                               ")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# =================================================================================")
_docx_paragraph_add_linebreak(dh)
_docx_image_add(dh,"${init_root}/genotypeqc2report/${genotypeqc2report_infile}/${genotypeqc2report_outfile}_bim2ancestry_pca1-pca2.png")
_docx_image_add(dh,"${init_root}/genotypeqc2report/${genotypeqc2report_infile}/${genotypeqc2report_outfile}_bim2ancestry_pca1-pca3.png")
_docx_image_add(dh,"${init_root}/genotypeqc2report/${genotypeqc2report_infile}/${genotypeqc2report_outfile}_bim2ancestry_pca2-pca3.png")
_docx_add_pagebreak(dh,)
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# =================================================================================")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# CHROMOSOME DISTRIBUTIONS")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# =================================================================================")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "Plotting the distribution of markers by chromosome is a useful sanity-check on your")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "data to make sure all chromosomes were included in the final report.")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# =================================================================================")
_docx_paragraph_add_linebreak(dh)
_docx_image_add(dh,"${init_root}/genotypeqc2report/${genotypeqc2report_infile}/${genotypeqc2report_outfile}-snp-by-chromosome.png")

_docx_add_pagebreak(dh,)
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# =================================================================================")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# KINSHIP DISTRIBUTIONS")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# =================================================================================")
_docx_paragraph_add_linebreak(dh) 
_docx_paragraph_add_text(dh, "# > kinship .................. dup = ${genotypeqc2report_kind}; 1st = ${genotypeqc2report_kinf}; 2nd = ${genotypeqc2reportkin_s}; 3rd = ${genotypeqc2report_kint}")
_docx_paragraph_add_linebreak(dh) 
_docx_paragraph_add_text(dh, "# > relPairs ................. ${genotypeqc2report_outfile}_bim2graphkin0.relPairs")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# =================================================================================")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "These plots are created from the *.kin0 files created using the --make-king-table  ")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "flag in plink2. The following plots show relatedness as defined by this kinship    ")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "matrix. Approximation of kinship thresholds (duplicates, 2nd and 3rd degree        ")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "relatives) are given by horizontal red lines. As part of the quality control       ")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "pipeline it was decided not to remove -related- individuals at this stage. The     ")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "exception are individuals with excessive relatedness (indicative of genotyping     ")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "error). The *.relPairs files can be used to define related pairs.                  ")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "A subset of unrelated individuals can be generated using the function              ")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "-bim2unrelated-.                                                                   ")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# =================================================================================")
_docx_paragraph_add_linebreak(dh)
_docx_image_add(dh,"${init_root}/genotypeqc2report/${genotypeqc2report_infile}/${genotypeqc2report_outfile}_bim2graphkin0-hist.png")
_docx_paragraph_add_linebreak(dh)
_docx_image_add(dh,"${init_root}/genotypeqc2report/${genotypeqc2report_infile}/${genotypeqc2report_outfile}_bim2graphkin0-scatter.png")

_docx_add_pagebreak(dh,)
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# =================================================================================")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# HETEROZYGOSITY DISTRIBUTIONS")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# =================================================================================")
_docx_paragraph_add_linebreak(dh) 
_docx_paragraph_add_text(dh, "# > --het (threshold) ........ ${genotypeqc2report_het}")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# =================================================================================")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "These plots are created from the *.het files created using the --het flag in plink.")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "The --het flag computes observed and expected autosomal homozygous genotype counts ")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "for each sample, and reports method-of-moments F coefficient estimates i.e.")
_docx_paragraph_add_linebreak(dh) 
_docx_paragraph_add_text(dh, " ")
_docx_paragraph_add_linebreak(dh) 
_docx_paragraph_add_text(dh, "([observed hom. count] - [expected count])")
_docx_paragraph_add_linebreak(dh) 
_docx_paragraph_add_text(dh, "------------------------------------------")
_docx_paragraph_add_linebreak(dh) 
_docx_paragraph_add_text(dh, "([total observations] - [expected count])")
_docx_paragraph_add_linebreak(dh) 
_docx_paragraph_add_text(dh, " ")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "Excessive heterozygosity can be indicative of mixed DNA samples; excessive         ")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "homozygosity can be indicative of poor DNA quality and allele drop-out. In this    ")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "pipeline, we exclude individuals whose rate of heterozygosity deviates from the    ")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "population (of genotyping array) mean. Importantly, heterozygosity rates differ on ")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "genotyping arrays and within individuals of differing ancestries. Where there are  ")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "mixed ancestry samples, you may consider relaxing the threshold and excluding only ")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "those individuals with extreme deviations from the sample mean.")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# =================================================================================")
_docx_paragraph_add_linebreak(dh)
_docx_image_add(dh,"${init_root}/genotypeqc2report/${genotypeqc2report_infile}/${genotypeqc2report_outfile}_bim2graphhet.png")
_docx_add_pagebreak(dh,)
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# =================================================================================")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# HARDY WEINBERG DISTRIBUTIONS")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# =================================================================================")
_docx_paragraph_add_linebreak(dh) 
_docx_paragraph_add_text(dh, "# > --hardy (threshold) ...... ${genotypeqc2report_hardy}")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# =================================================================================")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "These plots are created using the --hardy flag in plink. The --hardy flag writes a ")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "list of genotype-counts and Hardy-Weinberg equilibrium exact test statistics to    ")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "plink.hwe. Importantly, we apply less stringent thresholds than in other routines  ")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "as case samples, by design, may incur HW bias due to disease selection.")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# =================================================================================")
_docx_paragraph_add_linebreak(dh)
_docx_image_add(dh,"${init_root}/genotypeqc2report/${genotypeqc2report_infile}/${genotypeqc2report_outfile}_bim2graphhwe.png")

_docx_add_pagebreak(dh,)
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# =================================================================================")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# MISSINGNESS BY INDIVIDUAL")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# =================================================================================")
_docx_paragraph_add_linebreak(dh) 
_docx_paragraph_add_text(dh, "# > --mind (threshold) ....... ${genotypeqc2report_mind}")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# =================================================================================")
_docx_paragraph_add_linebreak(dh)	
_docx_paragraph_add_text(dh, "These plots are created using the --missing flag in plink. The --missing flag      ")
_docx_paragraph_add_linebreak(dh)	
_docx_paragraph_add_text(dh, "produces sample-based and variant-based missing data reports.")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# =================================================================================")
_docx_paragraph_add_linebreak(dh)
_docx_image_add(dh,"${init_root}/genotypeqc2report/${genotypeqc2report_infile}/${genotypeqc2report_outfile}_bim2graphmiss-imiss.png")

_docx_add_pagebreak(dh,)
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# =================================================================================")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# MISSINGNESS BY SNP")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# =================================================================================")
_docx_paragraph_add_linebreak(dh) 
_docx_paragraph_add_text(dh, "# > --geno (threshold) ....... ${genotypeqc2report_geno}")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# =================================================================================")
_docx_paragraph_add_linebreak(dh)	
_docx_paragraph_add_text(dh, "These plots are created using the --missing flag in plink. The --missing flag      ")
_docx_paragraph_add_linebreak(dh)	
_docx_paragraph_add_text(dh, "produces sample-based and variant-based missing data reports. As part of the       ")
_docx_paragraph_add_linebreak(dh)	
_docx_paragraph_add_text(dh, "-genotypeqc- pipeline we apply two --geno parameters during the QC; firstly, a     ")
_docx_paragraph_add_linebreak(dh)	
_docx_paragraph_add_text(dh, "pre-clean to remove excessive missing then completing the pipeline with a stricter ")
_docx_paragraph_add_linebreak(dh)	
_docx_paragraph_add_text(dh, "limit to give a stringent QC.")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# =================================================================================")
_docx_paragraph_add_linebreak(dh)
_docx_image_add(dh,"${init_root}/genotypeqc2report/${genotypeqc2report_infile}/${genotypeqc2report_outfile}_bim2graphmiss-lmiss.png")

_docx_add_pagebreak(dh,)
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# =================================================================================")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# ALLELE FREQUENCY")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# =================================================================================")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# > --mac (threshold) ........ ${genotypeqc2report_mac}")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# =================================================================================")
_docx_paragraph_add_linebreak(dh) 
_docx_paragraph_add_text(dh, "These plots are created using the --freq counts flag. We do not prune on allele    ")
_docx_paragraph_add_linebreak(dh) 
_docx_paragraph_add_text(dh, "frequency, using a standard minor allele count of 5 as an alternative.")
_docx_paragraph_add_linebreak(dh)
_docx_paragraph_add_text(dh, "# =================================================================================")
_docx_paragraph_add_linebreak(dh)
_docx_image_add(dh,"${init_root}/genotypeqc2report/${genotypeqc2report_infile}/${genotypeqc2report_outfile}_bim2graphfrq.png")
_docx_paragraph_add_linebreak(dh)

// SAVING DOCX
_docx_save(dh, "${init_root}/genotypeqc/${genotypeqc2report_infile}/${genotypeqc2report_outfile}-quality-control-report.docx", 1)
_docx_close(dh)
end
