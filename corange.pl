#usage: perl corange.pl $catPossible $coPossible

use strict;
use warnings;

our @start;
our @end;
our @sortEnd;

open(DATA,"$ARGV[0]") or die "failed to open:$!";
open(OUT,">$ARGV[1]");
our $title;
our $flag = 1;
our $num=0;
our $count=0;
our $i=0;
our $seqName;
our $ctg;
our $tempstart;
our $tempend;
our $minstart;
our $maxend;



my $filename=$ARGV[0];
print "$filename\n";

if(-s $filename){
	print "111\n";
while(<DATA>){
  chomp;
  $title = $_;
  our @temp = split(/_/,$title);
	$seqName=$temp[0];
	$ctg=$temp[1];
	$start[$num]=$temp[2];
	$end[$num]=$temp[3];
		
	$num++;
  }
  our @sortStart = sort { $a<=>$b }@start;
  our $TotalNum = @sortStart;
  for($count=0;$count<$TotalNum;$count++){
  	for($i=0;$i<$TotalNum;$i++){
  		if($sortStart[$count]==$start[$i]){
  			$sortEnd[$count]=$end[$i];
  			}
  		}
  	}
  	print "@sortStart\n";
   	print "@sortEnd\n";
   $tempstart=$sortStart[0];
   $tempend=$sortEnd[0];
   $minstart=$sortStart[0];
   $maxend=$sortEnd[0];
  
  if($TotalNum==1){
  	print OUT $seqName."_".$ctg."_".$minstart."_".$maxend."_all\n";
  	}
 else{
  for($count=1;$count<$TotalNum;$count++){
  	$tempstart=&min($minstart,$sortStart[$count]);
  	$tempend=&max($maxend,$sortEnd[$count]);
  if ($sortStart[$count]-$sortEnd[$count-1]>18000 or $tempend-$tempstart>65000 ){
		print OUT $seqName."_".$ctg."_".$minstart."_".$maxend."_all\n";
	  $minstart=$sortStart[$count];
	  $maxend=$sortEnd[$count];
 	}
 	else{
 		$minstart=$tempstart;
 		$maxend=$tempend;
 		}
  	}
  	print OUT $seqName."_".$ctg."_".$minstart."_".$maxend."_all\n";  	
  }
}
else{
	print "There is no prophage seq!";
	} 
 
close(OUT);

sub min{
	my $a=$_[0];
	my $b=$_[1];
	if($a<$b){
		return $a;
		}
	else{
		return $b;
		}
	}
	
sub max{
	my $a=$_[0];
	my $b=$_[1];
	if($a>$b){
		return $a;
		}
	else{
		return $b;
		}
	}