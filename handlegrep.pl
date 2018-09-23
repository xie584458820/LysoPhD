
use strict;
use warnings;

my $num1=0;
my $num2=0;
my @temp;
my $id;
my $id1;
my $id2;

my $left1;
my $right1;
my $left2;
my $right2;

my @temp1;
my @temp2;

my $readlen;
my $alignlen;
my $percent;

my $line;
my $hitnum=0;
my @hit;

my $i=0;
my $j=0;
open(GRETEMP,"<$ARGV[0]");
open(TEMP,">$ARGV[2]");
open(GRE,">$ARGV[3]");

while(<GRETEMP>){
chomp;
$left1=0;
$left2=0;
$line=$_;
@temp=split(/\s+/,$_);
$id=$temp[0];
$id1=$id."&1";
$id2=$id."&2";
open(CON,"<$ARGV[1]");
while(<CON>){
	chomp;
	$percent=0;
	$readlen=0;
	$alignlen=0;
	
	if(/$id1/){
		@temp1=split(/\s+/,$_);
		$readlen=$temp1[2];
		$alignlen=$temp1[3];
		$percent=$alignlen/$readlen;
		if($percent>0.95){
			$left1=$temp1[6];
			}
		
			
		}
	if(/$id2/){
		@temp2=split(/\s+/,$_);
		$readlen=$temp2[2];
		$alignlen=$temp2[3];
		$percent=$alignlen/$readlen;
		if($percent>0.95){
			$left2=$temp2[6];
			}			
		}
	}

if(($left1!=0) and ($left2!=0)and ((  ($left1-$left2<2300) and ($left1-$left2>0) ) or ( ($left2-$left1<2300) and ($left2-$left1>0) ))){
	print "left1:$left1  left2:$left2\n";
	}
else{
	print TEMP "$line\n";
	}

close(CON);
	}
close(TEMP);
open(TEMP,"<$ARGV[2]");

while(<TEMP>){
	chomp;
	@temp=split(/\s+/,$_);
		$hit[$hitnum][0]=$temp[0];
	 	$hit[$hitnum][1]=$temp[1];
	 	$hit[$hitnum][2]=$temp[2];
	 	$hit[$hitnum][3]=$temp[3];
	 	$hitnum++;
	}
my $add=0;

		for($i=0;$i<$hitnum;$i++){
			for($j=0;$j<$hitnum;$j++){
			if($hit[$i][0] eq $hit[$j][0]){
			$add=$hit[$i][1]+$hit[$j][1];
			if($add==3){
				print GRE "$hit[$i][0]\t$hit[$i][1]\t$hit[$i][2]\t$hit[$i][3]\n";
				}
				

				
				}

				}
			}
	
	
	
close(GRETEMP);
close(TEMP);
close(GRE);	
	
	