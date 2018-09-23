use strict;
use warnings;
use Cwd;

our $title;
our $flag = 1;
our $count=1;
our $ifallcontig;
our $inputseqpath=$ARGV[0];
our $ctgphgnumpath;
our $possiblepath=$ARGV[2];
our $matrixpath=$ARGV[3];
our $matrixfile;
our $ctg;
our $seqpath;
our $ctgphgnum;
our $dir=getcwd;
our $contigpath;

our ($path_programm) = ($0=~/^(.*)\/extractphage.*$/); 

if(! -d $ARGV[1]){
	system("mkdir -p $ARGV[1]");
	}
if(! -d $ARGV[2]){
	system("mkdir -p $ARGV[2]");
	}
if(! -d $ARGV[3]){
	system("mkdir -p $ARGV[3]");
	}


if(-s $inputseqpath){
open(SEQ,$ARGV[0]) or die "failed to open:$!";
while(<SEQ>){
	chomp;
	$title=$_;
		$flag += 1;
		
		if(/^>/){
		our @temp = split(/_/,$title);
		print "temp: @temp\n";
	  $ctg=$temp[2];
		if(/^>allcontig/){
			$flag = 0;
			$ifallcontig=1;
			$ctgphgnumpath=$ARGV[1]."/".$ctg.".txt";
			$seqpath=$possiblepath."/".$ctg.".fa";  
			$matrixfile=$matrixpath."/".$ctg.".txt"; 
			
			print "The ctgphgnumpath is: $ctgphgnumpath\n";
			
			print "The seqname is: $seqpath\n";
			
			open(NUM,">$ctgphgnumpath") or die "failed to open:$!";
			
			print NUM ">$ctg"."_0\n";
			close(NUM);
			}
		else{
		 $flag = 0;
		 $ifallcontig=0;
		 $seqpath=$possiblepath."/".$ctg."_$count".".fa"; 
		 $matrixfile=$matrixpath."/".$ctg."_".$count.".txt"; 
		 print "The seqname is: $seqpath\n";
			}
		$count++;
		}
		$contigpath="$dir/contig/$ctg";
		if($flag == 1){
			if($ifallcontig==1){
		open(POS,">$seqpath"); 
		print POS $_;
		close(POS);		
		system("$path_programm/repeatfind $seqpath 13 $matrixfile $contigpath");
	}
	else{
		open(POS,">$seqpath"); 
		print POS $_;
		close(POS);		
		system("$path_programm/repeatfind $seqpath 13 $matrixfile $contigpath");
		
		}
												 }
			else{
				
				}
			
		
		}
		
		$count-=1;
		if($ifallcontig==0){
			print "111";
			$ctgphgnumpath=$ARGV[1]."/".$ctg.".txt";
			open(NUM,">$ctgphgnumpath") or die "failed to open:$!";
			print NUM "$ctg"."_".$count."\n";
			close(NUM);
			}
}