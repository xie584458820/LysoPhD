
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
our ($seqName) = ($fileName =~ /^(.*).table$/);
our @groupOrf = 00000;
our @eachline;
our @seperateORF;
our $ctg = 00000;
our $ctgTime = 0;
our $count=0; 
our $flag1=0;

while(<DATA>){
  chomp;
  @eachline = split(/\s+/,$_); 
  $ctg=$eachline[1];

if($ctgTime==0){
	@seperateORF = split(/orf/,$eachline[0]);
	$flag1=0;
		}
if(/lysin/i)
{	print "Find lysin!\n";
	@seperateORF = split(/orf/,$eachline[0]);
	$groupOrf[$ctgTime]=$seperateORF[1];
	$ctgTime +=1;	}
if(/capsid/i)
{	print "Find capsid!\n";
	@seperateORF = split(/orf/,$eachline[0]);
	$groupOrf[$ctgTime]=$seperateORF[1];
	$ctgTime +=1;	}
if(/terminase/i)
{	print "Find terminase!\n";
	@seperateORF = split(/orf/,$eachline[0]);
	$groupOrf[$ctgTime]=$seperateORF[1];
	$ctgTime +=1;	}
if(/tail/i)
{	print "Find tail!\n";
	@seperateORF = split(/orf/,$eachline[0]);
	$groupOrf[$ctgTime]=$seperateORF[1];
	$ctgTime +=1;	}
if(/holin/i)
{	print "Find holin!\n";
	@seperateORF = split(/orf/,$eachline[0]);
	$groupOrf[$ctgTime]=$seperateORF[1];
	$ctgTime +=1;	} 
if(/fiber/i)
{	print "Find holin!\n";
	@seperateORF = split(/orf/,$eachline[0]);
	$groupOrf[$ctgTime]=$seperateORF[1];
	$ctgTime +=1;	}
} 
@groupOrf = (@groupOrf,$ctg);
print "The groupOrf is: ".@groupOrf."\n";
&OneDimDBSCAN(@groupOrf);
close(DATA);
close(OUT1);


sub OneDimDBSCAN{
  our $Eps = 25000;
  our $ctg = pop(@_);
  our @sortData = sort @_;
  our $TotalNum = @sortData;
  our $count;
  our @sepORF;
  our $subGroup = 1;
  our $startGroup = 0;
  our $start=0;
  our $end=0;
  our @eachline2;
  our @sepORF2;
  our @orfposition1;
  our @orfposition2;
  print "The TotalNum is: ".$TotalNum."\n";
  if($TotalNum==1){
  	print "TotalNum==1\n";
  	open(DATA1,"$ARGV[0]") or die "failed to open:$!";
  	while(<DATA1>){
			chomp;
			@eachline = split(/\s+/,$_);
			@sepORF = split(/orf/,$eachline[0]);
    if ($sepORF[1] == $sortData[0] and $eachline[1] eq $ctg){
  	if($eachline[2]<$eachline[3]){$start = $eachline[2];$end=$eachline[3];} 
		else{$start = $eachline[3];$end = $eachline[2];}
	}
}	 

 if($start==0 and $end==0){		} 
		else{	
		print OUT1 $seqName."_".$ctg."_".$start."_".$end."_all\n";
    }
   close(DATA1); 
  }
  else{
  	
  	for ($count=0;$count<$TotalNum;$count++){
    open(DATA2,"$ARGV[0]") or die "failed to open:$!";
    while(<DATA2>){
    	chomp;
    	@eachline2 = split(/\s+/,$_);
    	@sepORF2 = split(/orf/,$eachline2[0]);
		if($sepORF2[1] == $sortData[$count] and $eachline2[1] eq $ctg){
			$orfposition1[$count]=$eachline2[2];
			$orfposition2[$count]=$eachline2[3];
			} 	
    	}
    	close(DATA2);
    }
  print "TotalNum>1\n";
  for ($count=1;$count<$TotalNum;$count++){
	if ($orfposition1[$count]-$orfposition1[$count-1]<=$Eps){
		print "sortData[count]-sortData[count-1]<=Eps\n";
		print "The sortData[count] is: ".$sortData[$count]."\n";
		print "The sortData[count-1] is: ".$sortData[$count-1]."\n";
		$subGroup +=1;
 	}
  else{
	#-----------for one group in a contig, record its start and end-----------
	  print "sortData[count]-sortData[count-1]<=Eps\n";
		open(DATA1,"$ARGV[0]") or die "failed to open:$!";
		while(<DATA1>){
			chomp;
			@eachline = split(/\s+/,$_);
			@sepORF = split(/orf/,$eachline[0]);
			if ($sepORF[1] == $sortData[$startGroup] and $eachline[1] eq $ctg){
			  if($eachline[2]<$eachline[3]){$start = $eachline[2]; }  
			  else{$start = $eachline[3];}
			  }
			if ($sepORF[1] == $sortData[$count-1] and $eachline[1] eq $ctg){
				if($eachline[2]<$eachline[3]){$end = $eachline[3]; }
				else{$end=$eachline[2];}
			}
		} 
		if($start==0 and $end==0)
{		}
  	else{
		print OUT1 $seqName."_".$ctg."_".$start."_".$end."_all\n";
	}
		close (DATA1); 
   #--------------------------------------------------------------------------
		$subGroup = 1;
		$startGroup = $count;
	}
  }
  
    if($startGroup==$TotalNum-1){
	#-----------for one group in a contig, record its start and end-----------
	  print "startGroup==TotalNum-1\n";
		open(DATA1,"$ARGV[0]") or die "failed to open:$!";
		while(<DATA1>){
			chomp;
			@eachline = split(/\s+/,$_);
			@sepORF = split(/orf/,$eachline[0]);
			print "The startGroup is:".$startGroup."\n";
			print "The sepORF is:".$sepORF[1]."\n";
			print "The sortData[startGroup] is:".$sortData[$startGroup]."\n";
			if ($sepORF[1] == $sortData[$startGroup] and $eachline[1] eq $ctg){
				print "The last orf is:".$sepORF[1]."\n";
				if($eachline[2]<$eachline[3])
				{$start = $eachline[2];$end = $eachline[3]}
				else
				{$start = $eachline[3];$end = $eachline[2];}

		}
	} 
	if($start==0 and $end==0)
{		}
	else{ 
		print OUT1 $seqName."_".$ctg."_".$start."_".$end."_all\n";
	}
		close(DATA1);
}
	if($startGroup<$TotalNum-1){
		print "startGroup<TotalNum-1\n";
		open(DATA1,"$ARGV[0]") or die "failed to open:$!";
		while(<DATA1>){
			chomp;
			@eachline = split(/\s+/,$_);
			@sepORF = split(/orf/,$eachline[0]);
			if ($sepORF[1] == $sortData[$startGroup] and $eachline[1] eq $ctg){
			  if($eachline[2]<$eachline[3]){$start = $eachline[2];
 } 
			  else{$start = $eachline[3];}
			  }
			if ($sepORF[1] == $sortData[$TotalNum-1] and $eachline[1] eq $ctg){
				if($eachline[2]<$eachline[3]){$end = $eachline[3]; }
				else{$end=$eachline[2];}
			}
		} 
		if($start==0 and $end==0)
{		}
    else
 {
 		print OUT1 $seqName."_".$ctg."_".$start."_".$end."_all\n";
	}
		close(DATA1); 		
		
		}
}
  

	
	 
	
	
	
	
}