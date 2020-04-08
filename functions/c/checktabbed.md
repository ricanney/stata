# checkfile
__Author:__ Richard Anney
__Edited:__ 31 Mar 2020
__Current Version:__ v2

__Overview__
This function runs a perl script ``tabbed.pl``. ``tabbed.pl`` converts white-space delimited files to tab-delimited. 

__Dependencies__
``checkfile``
Note that your profile.do file must contain the following ``macro`` defining the location of the ``tabbed.pl`` file

```
global tabbed perl H:/software/perl/code/tabbed.pl
```

__tabbed.pl__
```
open (IN, $ARGV[0]);
open (OUT, ">$ARGV[0].tabbed");
while (<IN>) {
	$_ =~ s/\s+/\t/g;
	$_ =~ s/^\t//g;
	print OUT "$_\n";
  }
```

__Syntax__
```
checktabbed
```




