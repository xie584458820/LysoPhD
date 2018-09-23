use strict;
use warnings;
use List::Util qw/sum/;
use List::Util qw/max min/;
use Data::Dumper qw(Dumper);
use 5.010;

open(DATA,"$ARGV[0]") or die "failed to open:$!";
open(OUT1,">$ARGV[1]");

our $wholePath = "$ARGV[0]";
our $fileName = "";
#get the file name without path
#print "the wholePath is:".$wholePath."\n";
if ($wholePath =~ m!^(.*)\/([^\/]*)$!){
	$fileName = $2;
} else {
	$fileName = $wholePath;
} 
#print "the filename is:".$fileName."\n";
#-----------------------------------
our ($seqName) = ($fileName =~ /^(.*).txt.title$/);
our @groupOrf = 00000;
our @eachline;
our @seperateORF;
our $ctg = 00000;
our $ctgTime = 0;

while(<DATA>){
  chomp;
  @eachline = split(/\s+/,$_);

  if ($eachline[2] ne $ctg){
  #record the maximum and minimum of the normal distribution ctgs
	if ($ctgTime>5){
	#--------record each group in each contig--------------
	  @groupOrf = (@groupOrf,$ctg);
	  &OneDimDBSCAN(@groupOrf);
	}

	$ctg = $eachline[2];  #record the last line
	@groupOrf=();
	$ctgTime = 0;
	@seperateORF = split(/orf/,$eachline[1]);
	$groupOrf[$ctgTime]=$seperateORF[1];
    $ctgTime=1;
  }
  else{
	@seperateORF = split(/orf/,$eachline[1]);
	$groupOrf[$ctgTime]=$seperateORF[1];
	$ctgTime +=1;	
  } 
}
if ($ctgTime>5){
	#--------record each group in each contig--------------
	  @groupOrf = (@groupOrf,$ctg);
	  &OneDimDBSCAN(@groupOrf);
}
close(DATA);
close(OUT1);

sub OneDimDBSCAN{
  our $Eps = 5;
  our $ctg = pop(@_);
  our @sortData = sort @_;
  our $TotalNum = @sortData;
  our $count;
  our $start;
  our $end;
  our @sepORF;
  our $subGroup = 1;
  our $startGroup = 0;
  for ($count=1;$count<$TotalNum;$count++){
	if ($sortData[$count]-$sortData[$count-1]<=$Eps){
		$subGroup +=1;
	}
  else{
	  if ($subGroup>=5){
	#-----------for one group in a contig, record its start and end-----------
		open(DATA1,"$ARGV[0]") or die "failed to open:$!";
		my $start =0;
		my $end =0;
		while(<DATA1>){
			chomp;
			@eachline = split(/\s+/,$_);
			@sepORF = split(/orf/,$eachline[1]);
			if ($sepORF[1] == $sortData[$startGroup] and $eachline[2] eq $ctg){
				$start = $eachline[3];
			}
			if ($sepORF[1] == $sortData[$count-1] and $eachline[2] eq $ctg){
				$end = $eachline[4];
			}
		}
		print OUT1 $seqName."_".$ctg."_".$start."_".$end."_all\n";
		close (DATA1);
		}
    #--------------------------------------------------------------------------
		$subGroup = 1;
		$startGroup = $count;
	}
  }
  if ($subGroup>=5){
	#-----------for one group in a contig, record its start and end-----------
		open(DATA1,"$ARGV[0]") or die "failed to open:$!";
		my $start =0;
		my $end =0;
		while(<DATA1>){
			chomp;
			@eachline = split(/\s+/,$_);
			@sepORF = split(/orf/,$eachline[1]);
			if ($sepORF[1] == $sortData[$startGroup] and $eachline[2] eq $ctg){
				$start = $eachline[3];
			}
			if ($sepORF[1] == $sortData[$count-1] and $eachline[2] eq $ctg){
				$end = $eachline[4];
			}
		}
		print OUT1 $seqName."_".$ctg."_".$start."_".$end."_all\n";
		close (DATA1);
		}
   $subGroup = 1;
   $startGroup = $count;
}