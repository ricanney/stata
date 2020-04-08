*! 1.0.1 Richard Anney 15apr2019
program do_temp_script
syntax 

outsheet script using temp_script.do, non noq replace
do temp_script.do
erase temp_script.do 
end; 
