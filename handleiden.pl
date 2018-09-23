use strict;
use warnings;
my $num1=0;
my $num2=0;
my @temp;
my $iden;
my $idennum;
my @temp1;
my $slen;
my $alignlen;
my $percent;
my $flag=$ARGV[2];
open(DATA,"<$ARGV[0]");
open(OUT,">$ARGV[1]");
if($flag==0){
while(<DATA>){
	chomp;
	if(/^ref.fa/){
		$num1++;
		@temp=split(/\s+/,$_);
		$iden=$temp[1];
		$slen=$temp[2];
		$alignlen=$temp[3];
		$percent=$alignlen/$slen;
		@temp1=split(/\./,$iden);
		$idennum=$temp1[0];
		if($idennum>=99 and (($slen<100 and $percent>0.95)or($slen>100 and $percent>0.93))){
			print OUT $_."\n";
			$num2++;
			}
		}
	}
}
else{
while(<DATA>){
	chomp;
	if(/^ref.fa/){
		$num1++;
		@temp=split(/\s+/,$_);
		$iden=$temp[1];
		$slen=$temp[2];
		$alignlen=$temp[3];
		$percent=$alignlen/$slen;
		@temp1=split(/\./,$iden);
		$idennum=$temp1[0];
		if($idennum>=99){
			print OUT $_."\n";
			$num2++;
			}
		}
	}	
	
	}


	close DATA;
	close OUT;