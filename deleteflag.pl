use strict;
use warnings;


my $flagposition=$ARGV[2];
my $taglengthtest=$ARGV[3];
my $differ;
my $taglength;
my $forward;
my $later;
my $jointwo;

open(DATA,"<$ARGV[0]");
open(OUT,">$ARGV[1]");

while(<DATA>){
	if(/>/){
	#print "line    :$_\n";
	$taglength=length($_);
	$differ=$taglength-$taglengthtest;
	$forward=substr($_,0,$flagposition+$differ);
	#print "forward :$forward\n";
	$later=substr($_,$flagposition+$differ+1);
	#print "later :$later\n";
	$jointwo=$forward.$later;
  $jointwo=~s/\s//g;
  #print "jointwo :$jointwo\n";
  print OUT "$jointwo\n"; 
}
else{
	print OUT "$_";
	}
	}