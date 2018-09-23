use strict;
use warnings;

open(DATA,"$ARGV[0]") or die "failed to open:$!";
open(OUT,">$ARGV[1]");
our $num = 0;

while(<DATA>){
	chomp;
	if (/^>/){
		if ($num != 0){
			print OUT "\n".$_."\n";
		}else{
			print OUT $_."\n";
			$num += 1;
		}
	}elsif(/^(\w+)\W*$/){
		print OUT $1;
	}
}

close(DATA);
close(OUT);