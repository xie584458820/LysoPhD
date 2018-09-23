use strict;
use warnings;
use POSIX;


open(DATA,"$ARGV[1]") or die "failed to open:$!";
open(OUT,">$ARGV[2]");
our $title;
our $flag = 1;
our $flag1 = 1;
our $string;
our $contiglen;
our $count=1;
our $half;
our $halftemp;
our @temp1;
our @temp2;
our $ctg;
our $coPossible=$ARGV[1];

if(-s $coPossible){
open(SEQ1,"$ARGV[0]") or die "failed to open:$!";
while(<SEQ1>){
	chomp;
	$title=$_;
	$flag1+=1;
	if(/^>/){
		@temp1 = split(/\s+/,$title);
		@temp2 = split(/>/,$temp1[0]);
		$ctg=$temp2[1];
		$flag1=0;
		}
	if($flag1==1){
		$string=$_;
		}	
	}
$contiglen=length($string);
print "contiglen is: $contiglen\n";
if($contiglen<100000){
print OUT ">allcontig_".$ctg."_".$ctg."_0_".$contiglen."_all\n";
print OUT substr($string, 0, length($string))."\n";	
	}
else{
while(<DATA>){
  chomp;
  $title = $_;
  our @temp = split(/_/,$title);
  if($temp[2]==0 and $temp[3]==0)
  {
  	print "The coPossible is empty!\n";
  	}
  else{
  	print "Get coPossible!\n";
  open(SEQ,"$ARGV[0]") or die "failed to open:$!";
  while(<SEQ>){
  	chomp;
  	$flag += 1;
  	if (/^>$temp[1].*$/){
  		$flag = 0;
  		print OUT ">partcontig_".$title."\n"; 
  	}
  	if ($flag == 1){
  		
  		$halftemp=($temp[3]-$temp[2])/2;
  		$half=floor($halftemp);
  		if($temp[2]<=45000-$half){
  			print OUT substr($_, 0, 90000)."\n";
  		print "555\n";
  			}
  		if($temp[2]>45000-$half)
  		{
  		print "666\n";	
  		print "The half is $half\n";
  		print "The temp[2] is: $temp[2]\n";
  		print OUT substr($_, $temp[2]+$half-45000, 90000)."\n";
  	   }
  		
  		close(SEQ);
  		last;
  	}
  }
}
  
}
}
}

close(OUT);