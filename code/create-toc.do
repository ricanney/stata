cd H:/software/stata/data
dir2dta, dir(H:/software/stata/functions)
sort folder
sxpose, clear
gen _var99 = " "
aorder
gen folder = " "
for var _var1 - _var99: replace folder = folder + " " + X
replace folder = subinstr(folder,"  "," ",.)
keep folder
local folder = folder[1]
di "`folder'"
foreach i in `folder' {
	clear
	set obs 2
	gen v1 = ""
	replace v1 = "v 1" in 1
	replace v1 = "d Programs by Richard Anney, Cardiff University" in 2
	gen obs = _n
	replace obs = obs -2
	save H:/software/stata/functions/`i'/stata-toc.dta, replace
	files2dta, dir(H:/software/stata/functions/`i')
	split files,p(".")
	keep if files2 == "ado"
	keep files1
	sort files1
	gen obs = _n
	expand 10
	sort files
	egen seq= seq(),by(files)
	sum obs
	gen script = ""
	replace script = "import delim using H:/software/stata/functions/`i'/" + files1 + `".pkg, delim("$$") varnames(nonames) clear"' if seq == 1
	replace script = "keep if _n == 1" if seq == 2 
	replace script = "gen obs = 1" if seq == 3
	replace script = "append using H:/software/stata/functions/`i'/stata-toc.dta" if seq == 4
	replace script = "save H:/software/stata/functions/`i'/stata-toc.dta, replace" if seq == 5
	drop if script == ""
	outsheet script using create-toc-temp.do, non noq replace
	do create-toc-temp.do
	erase create-toc-temp.do
	sort obs v1
	replace v1 = "$$" + v1 if obs == 1
	replace v1 = subinstr(v1,"$$d","p",.)
	replace v1 = subinstr(v1,"ppd","p",.)
	outsheet v1 using H:/software/stata/functions/`i'/stata.toc, non noq replace
	erase H:/software/stata/functions/`i'/stata-toc.dta
	files2dta, dir(H:/software/stata/functions/`i')
	split files,p(".")
	keep if files2 == "ado"
	gen script = "net install " + files1 + ", from(https://raw.github.com/ricanney/stata/master/functions/`i') replace"
	keep script
	save H:\software\stata\code\net-install-`i'.dta, replace
	save H:\software\stata\code\net-install.dta, replace
	}
foreach i in `folder' {
	use H:\software\stata\code\net-install-`i'.dta, clear
	append using H:\software\stata\code\net-install.dta
	save H:\software\stata\code\net-install.dta, replace
	erase H:\software\stata\code\net-install-`i'.dta
	}
duplicates drop
sort script
outsheet script using H:\software\stata\code\net-install.do, non noq replace
erase H:\software\stata\code\net-install.dta
