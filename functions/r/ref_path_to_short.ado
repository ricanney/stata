*! 1.0.1 Richard Anney 15apr2019
*! 1.0.2 Richard Anney 30mar20202
* removed the creation of temp.do - *use global ref_short 	= ref1[1]* instead

program ref_path_to_short
syntax , ref(string asis)

qui { // define short ref name
	preserve
	clear
	set obs 1
	gen ref = ""
	replace ref =  "`ref'" in 1
	replace ref = subinstr(ref,"\","%%",.)
	replace ref = subinstr(ref,"/","%%",.)
	split ref, p("%%")
	gen ref999 = ""
	for var ref1-ref999: replace ref1 = X if X !=""
	keep ref1
	compress
		global ref_short 	= ref1[1]
	restore
	}
end; 

