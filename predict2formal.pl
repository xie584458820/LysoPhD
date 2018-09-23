
use strict;
use warnings;

open(DATA, "$ARGV[0]") or die "Cannot open:$!"; 
open(OUT, ">$ARGV[1]") or die "Cannot open:$!"; 
our $title;

while(<DATA>){
    chomp;
    if (/^>(\S+).*$/){
  	    $title = $1;
    }
	  if (/^orf/){
	  	  my @array = split;
	  	  print OUT "$array[0] $title $array[1] $array[2] $array[3] $array[4]\n";
	  }
}