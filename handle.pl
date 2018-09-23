use strict;
use warnings;

open("DATA","<$ARGV[0]");
open("OUT",">$ARGV[1]");
our $title;
our @temp;
our @temp1;
our $line;
our $flag=$ARGV[2];
our $splitnum;
our $split0;
our $jointemp;
our $joinall;
while(<DATA>){
	chomp;
	$title=$_;
	if(/>/){
		@temp=split(/\s+/,$_);
		$split0=$temp[0];
		#$splitnum=@temp;
		#for(my $i=0;$i<$splitnum;$i++){
		#	}
		$joinall=$split0."&".$flag;
		#print "$line\n";
		print OUT $joinall."\n";
		}
	else{
		print OUT $_."\n"; 
		}
	}
	
	close(DATA);
	close(OUT);