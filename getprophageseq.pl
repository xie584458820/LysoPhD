use strict;
use warnings;

open(DATA,"$ARGV[1]") or die "failed to open:$!";
open(OUT,">$ARGV[2]");
our $title;
our $flag = 1;

while(<DATA>){
  chomp;
  $title = $_;
  our @temp = split(/_/,$title);
  open(SEQ,"$ARGV[0]") or die "failed to open:$!";
  while(<SEQ>){
  	chomp;
  	$flag += 1;
  	if (/^>$temp[1].*$/){
  		$flag = 0;
  		print OUT ">".$title."\n";
  	}
  	if ($flag == 1){
  		print OUT substr($_, 1000-$temp[3], $temp[3]+100000)."\n";
  		close(SEQ);
  		last;
  	}
  }
}
close(OUT);