
use strict;
use warnings;

open(DATA,"$ARGV[0]") or die "failed to open:$!";
open(OUT1,">$ARGV[1]");
open(OUT2,">$ARGV[1].title");
our $num = 0;
our $numS = 0;
our $i = 0;
our @data;
our $flag;

while(<DATA>){
  chomp;
  $data[$num] = $_;
  if (/^Query= /){
  	$numS = $num;
	  $flag = 0;
  }
	if (/^ Identities = (\w+).*\((\w+)%.*Gaps.*$/){
	  $flag += 1;
	  if ($2>30 and $flag == 1){
	  	print OUT2 $data[$numS]."\n";
	  	for ($i=$numS;$i<=$num;$i++){
	  	  print OUT1 $data[$i]."\n";
	  	}
	  print OUT1 "\n"."\n"."\n"."\n";
	  }
	}
	$num += 1;
}
close(DATA);
close(OUT1);
close(OUT2);