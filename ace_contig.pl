use strict;
use warnings;


my $contignum=$ARGV[2];
my $flag=0;
my $path=$ARGV[1];
my $outpath=$ARGV[3];
my $tempfile=$path."/$contignum".".ace";
my @temp;
my @temp1;
my @temp2;
my $position;
my @contig;
my $fm_ctg;
my $to_ctg;
my $seq;

my $startnum;
my $endnum;

my @start;
my @endend;
my @self;


@contig=split /contig0+/,$contignum;
my $ctg=$contig[1];
print "The contig is $ctg\n";
open (DATA,"$ARGV[0]") or die "failed to open $ARGV[0]:$!";
open (TEMP,">$tempfile");
open (OUT,">$outpath");
while(<DATA>){
	chomp;
	if(/^CO/){
		$flag=0;
		}
	if(/$contignum/){
		$flag=1;
		}
	if($flag==1){
		print TEMP $_."\n";
		}
	}

close(TEMP);
open (RUN,"<$tempfile");

while(<RUN>){
	$seq = $_;
	  if (/^AF/){
  	@temp = split / /, $seq;
  	$position=$temp[3];
  	@temp1 = split /\./, $temp[1];
		if($temp1[-1]=~/fm/){
			@temp2=split /fm/, $temp1[-1];
			$fm_ctg=$temp2[1];
			if($fm_ctg != $ctg){
				if($temp[3]<1000){ push @start,$seq; }
				if($temp[3]>10000){ push @endend,$seq; }
				}
			else{
				push @self,$seq;
				}			
			}	
			
			
		if($temp1[-1]=~/to/){
			@temp2=split /to/, $temp1[-1];
								
			$to_ctg=$temp2[1];
			if($to_ctg != $ctg){
				if($temp[3]<1000){ push @start,$seq; }
				if($temp[3]>10000){push @endend,$seq; }
				}
				else{
					push @self,$seq;
					}
			}
  }
  }
	
if((@start>2) and (@endend>2)){
	print OUT "OnBacterialGenome\n";
	}

print OUT "\nCyclinzingReads:\n";
foreach my $k(@self){
	chomp $k;
		print OUT "$k\n";
	}
print OUT "\n\nStart:\n";
foreach my $i(@start){
	chomp $i;
		print OUT "$i\n";
	}

print OUT "\n\nEnd:\n";
foreach my $j(@endend){
	chomp $j;
		print OUT "$j\n";
	}

	
close(RUN);
close(DATA);
