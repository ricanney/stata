*! 1.0.1 Richard Anney 23jan2019

program checkfolder
syntax , folder(string asis) 
	
qui { 
	capture confirm folder "`folder'"
	if _rc==0 {
		noi di as text"# > "as input"checkfolder "as text"................................. located" as result" `folder'"
		}
	else {
		noi di as text"# > "as input"checkfolder "as text"........................... " as error "cannot locate" as result" `folder'"
		exit
		}
end;
	
