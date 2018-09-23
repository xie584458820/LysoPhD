use strict;
use warnings;
my $num1=0;
my $num2=0;
my @temp;
my $iden;
my $idennum;
my @temp1;
my $qlen;
my $alignlen;
my $percent;

open(DATA,"<$ARGV[0]");
open(OUT,">$ARGV[1]");
while(<DATA>){
	chomp;
	if(/^Query_1/){
		$num1++;
		@temp=split(/\s+/,$_);
		$iden=$temp[1];
		$qlen=$temp[2];
		$alignlen=$temp[3];
		$percent=$alignlen/$qlen;
		@temp1=split(/\./,$iden);
		$idennum=$temp1[0];
		if($idennum>=97){
			print OUT $_."\n";
			$num2++;
			}
		}
	}
print "$num1\n";
print "$num2\n";
	close DATA;
	close OUT;