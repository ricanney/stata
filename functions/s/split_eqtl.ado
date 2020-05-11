program  split_eqtl
syntax

split eqtl, p(" ")
split eqtl1,p("-")
replace eqtl4 = subinstr(eqtl4,"(","",.)
replace eqtl4 = subinstr(eqtl4,")","",.)
replace eqtl4 = subinstr(eqtl4,";","",.)
split eqtl5,p("(")
replace eqtl52 = subinstr(eqtl52,")","",.)
rename(eqtl11 eqtl12 eqtl2 eqtl3 eqtl4 eqtl51 eqtl52 eqtl8) (eqtl_snp eqtl_a1 eqtl_symbol eqtl_tissue eqtl_super_tissue eqtl_beta eqtl_se eqtl_p)
drop eqtl1 eqtl5-eqtl7
for var eqtl_beta eqtl_se eqtl_p: destring X, replace

end;

