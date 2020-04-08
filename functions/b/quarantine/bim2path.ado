*! 1.0.2 Richard Anney 24apr2019
* clean up-front functions
*! 1.0.1 Richard Anney 10dec2018
* - fix end line in windows / unix - 10dec2018

program  bim2path
local    function bim2path
local    `function'_version v2
syntax , bim(string asis) 

preserve
qui { // define bim_path
	global bim_path `bim'
	}
qui { // define bim_folder
	clear
	set obs 1
	gen bim = ""
	replace bim =  "`bim'" in 1
	replace bim = subinstr(bim,"\","%%",.)
	replace bim = subinstr(bim,"/","%%",.)
	split bim, p("%%")
	drop bim
	sxpose, clear
	drop if _var1 == ""
	gen obs = _n
	gsort -obs
	replace obs = 0 in 1
	gsort obs
	drop if obs == 0
	keep _v
	rename _ x
	sxpose, clear
	gen _var99 = ""
	gen folder = "/"
	for var _var1 - _var99: replace folder = folder + X + "/"
	replace folder = subinstr(folder,"//","/",.)
	replace folder = substr(folder,1,length(folder)-1)
	keep folder
	gen script = "global bim_folder " + folder
	outsheet script using temp.do, non noq replace
	do temp.do
	erase temp.do	
	}


qui { // define bim_file
	bim_path_to_short, bim(`bim')
	global bim_file $bim_short
	}
di as text"# > ........ full path $bim_path"
di as text"# > ........... folder $bim_folder"
di as text"# > ............. file $bim_file"
restore
end;
