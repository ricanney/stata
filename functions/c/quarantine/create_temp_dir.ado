*! 1.0.0 Richard Anney 05nov2018
program create_temp_dir
syntax 

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

preserve

qui { // 1 - introduction
	noi di as text"# > create_temp_dir ....... create temporary directory in "as result"`c(pwd)'"
	}
qui { // 2 - set seed
	clear
	set obs 52
	gen folderRandom = ""
	qui { // random seed from - https://www.stata.com/statalist/archive/2008-10/msg00032.html
			tokenize "`c(current_date)'" ,parse(" ")
			local seed_1 "`1'"
			tokenize "`c(current_time)'" ,parse(":")
			local seed_2 "`1'`3'`5'"
			local seed_final "`seed_1'`seed_2'"
			set seed `seed_final'
			}
	}
qui { // 3 - generate folder name
	foreach num of num 1/10 {
			gen a = uniform()
			gen b = "a"
			replace b ="b" in 2
			replace b ="c" in 3
			replace b ="d" in 4
			replace b ="e" in 5
			replace b ="f" in 6
			replace b ="g" in 7
			replace b ="h" in 8
			replace b ="i" in 9
			replace b ="j" in 10
			replace b ="k" in 11
			replace b ="l" in 12
			replace b ="m" in 13
			replace b ="n" in 14
			replace b ="o" in 15
			replace b ="p" in 16
			replace b ="q" in 17
			replace b ="r" in 18
			replace b ="s" in 19
			replace b ="t" in 20
			replace b ="u" in 21
			replace b ="v" in 22
			replace b ="w" in 23
			replace b ="x" in 24
			replace b ="y" in 25
			replace b ="z" in 26
			replace b ="A" in 27
			replace b ="B" in 28
			replace b ="C" in 29
			replace b ="D" in 30
			replace b ="E" in 31
			replace b ="F" in 32
			replace b ="G" in 33
			replace b ="H" in 34
			replace b ="I" in 35
			replace b ="J" in 36
			replace b ="K" in 37
			replace b ="L" in 38
			replace b ="M" in 39
			replace b ="N" in 40
			replace b ="O" in 41
			replace b ="P" in 42
			replace b ="Q" in 43
			replace b ="R" in 44
			replace b ="S" in 45
			replace b ="T" in 46
			replace b ="U" in 47
			replace b ="V" in 48
			replace b ="W" in 49
			replace b ="X" in 50
			replace b ="Y" in 51
			replace b ="Z" in 52
			sort a
			replace folderRandom = folderRandom + b[1]
			drop a b
			}
	gen a = "global temp_dir  " + folderRandom
	outsheet a using _x.do, non noq replace
	do _x.do
	erase _x.do
	}
qui { // 4 - create folder
	di as text"# >> creating random folder"
	!mkdir ${temp_dir}
	cd ${temp_dir}
	}
noi di as text"# > create_temp_dir ..................................... directory created"
noi di as text"# > create_temp_dir ....................... new directory " as result"`c(pwd)'"
noi di as text"# > create_temp_dir ...................... \${temp_dir} is " as result"${temp_dir}"
noi di as text"#########################################################################"
noi di as text"# Completed: $S_DATE $S_TIME"
noi di as text"#########################################################################"
restore
end;
	
