program  gene_plotorder
local    function gene_plotorder
local    `function'_version v1
syntax, from(string asis) to(string asis)

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

qui { // define plot order
	save `function'_temp.dta,replace
	keep start end label 
	drop if start == .
	gen size = `to' - `from'
	gen buffer = size / 20
	drop if start == .
	egen group = group(start end label)
	egen x = seq(),by(group)
	keep if x == 1
	sum buffer
	gen buffer_start = round(start - `r(max)', 1000) - 10000
	gen buffer_end   = round(end + `r(max)', 1000) + 10000
	foreach i in start end {
	gen string_`i' = buffer_`i'
		tostring string_`i', replace
		}
	sum group
	local group_n = `r(max)'
	gen order = .
	foreach round of num 1 / `group_n' {
		egen y = seq(),by(order) 
		replace order = `round' if order == . & y == 1
		drop y
		foreach group of num 1 / `group_n' {
			sum buffer_end if order == `round'
			if `r(N)' !=0 {
				replace order = `round' if buffer_start > `r(max)' & order == . & group == `group'
				}
			}
		}
	keep order label start end
	sum order
	if `r(max)' > 100 {
		global `function'_f half_tiny
		global `function'_e *2
		}
	else if `r(max)' > 50 {
		global `function'_f half_tiny
		global `function'_e *2.5
		}
	else if `r(max)' > 20 {
		global `function'_f half_tiny
		global `function'_e *3
		}
	else {
		global `function'_f tiny
		global `function'_e *4
		}
	merge 1:m  label start end using `function'_temp.dta
	erase `function'_temp.dta
	drop _m
	}		
end;
