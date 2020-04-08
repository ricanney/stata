*! 1.0.0 Richard Anney 05nov2018
program delete_temp_dir
syntax 

qui { // 
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


	noi di as text"# > delete_temp_dir .............. temporary directory is "as result"${temp_dir}"
	noi di as text"# > delete_temp_dir ................ current directory is "as result"`c(pwd)'"
	clear
	set obs 1
	gen os = "`c(os)'"
	if os == "Unix" { 
		cd ..
		!rm -r ${temp_dir}
		}
	else if os == "Windows" {
		cd ..
		!rd /s /q "${temp_dir}" 
		}
	noi di as text"# > delete_temp_dir ..................................... directory deleted"
	noi di as text"# > delete_temp_dir ............... current directory now "as result"`c(pwd)'"
	noi di as text"#########################################################################"
	noi di as text"# Completed: $S_DATE $S_TIME"
	noi di as text"#########################################################################"
	}
end;
	
