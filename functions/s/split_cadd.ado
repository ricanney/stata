program  split_cadd
syntax

split cadd, p(" ")
split cadd1,p("-")
replace cadd6 = subinstr(cadd6,")","",.)
replace cadd6 = subinstr(cadd6,";","",.)
capture confirm variable cadd7
if !_rc {
	split cadd7,p("/"")")
	rename(cadd11 cadd12 cadd2 cadd3 cadd6 cadd71 cadd72) (cadd_snp cadd_a1 cadd_consequence cadd_consscore cadd_phred cadd_polyphencat cadd_siftcat)
	drop cadd1 cadd4 cadd5 cadd7
	}
else {
	gen cadd71 = ""
	gen cadd72 = ""
	rename(cadd11 cadd12 cadd2 cadd3 cadd6 cadd71 cadd72) (cadd_snp cadd_a1 cadd_consequence cadd_consscore cadd_phred cadd_polyphencat cadd_siftcat)
	drop cadd1 cadd4 cadd5 
	}
for var cadd_consscore cadd_phred: destring X, replace

end;

