*! 1.0.1 Richard Anney 23jan2019

program checkfile
syntax , file(string asis) 
qui {	
	capture confirm file "`file'"
	if _rc==0 {
		noi di as text"# > .......... located" as result" `file'"
		}
	else {
		noi di as text"# > .... " as error "cannot locate" as result" `file'"
		exit
		}
	}
end;
	
