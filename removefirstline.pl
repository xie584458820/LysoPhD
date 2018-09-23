
use strict;
use warnings;

open(Inputfile,"$ARGV[0]")or die "failed to open:$!";
open(Outputfile,">$ARGV[1]") or die "can't open file!";
readline Inputfile;
while(<Inputfile>){
	print Outputfile $_;
	}
	close Outputfile;