#!/usr/bin/perl  -w
use strict;
die "perl $0 <input cds file> <output pep file> <translation table No(1 .. 12)>\n" unless (@ARGV==3);

my $in_path = shift @ARGV;
my $out_path = shift @ARGV;
my $transl_table = shift @ARGV;

open IN,"$in_path" or die "Cannot open input file: $!\n";
open OUT,">$out_path" or die "Cannot open output file: $!\n";

my ($i,$j,$k);
my ($in,$seq,$out);
my ($aas,$starts,$base1,$base2,$base3);
my $temp;
my $flag=0;
my @mod;
my (%hash,%hashStarts);

##You can get the translation table form http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c
#The Standard Code
if($transl_table == 1)
{
	$aas    = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
	$starts = "---M---------------M---------------M----------------------------";
	$base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG";
	$base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG";
	$base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG";
}

#The Vertebrate Mitochondrial Code
elsif($transl_table == 2){
	$aas    = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG";
	$starts = "--------------------------------MMMM---------------M------------";
	$base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG";
	$base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG";
	$base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG";
}

#The Yeast Mitochondrial Code
elsif($transl_table == 3){
	$aas    = "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
	$starts = "----------------------------------MM----------------------------";
	$base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG";
	$base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG";
	$base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG";
}

#The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
elsif($transl_table == 4){
	$aas    = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
	$starts = "--MM---------------M------------MMMM---------------M------------";
	$base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG";
	$base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG";
	$base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG";
}

#he Invertebrate Mitochondrial Code
elsif($transl_table == 5){
	$aas    = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG";
	$starts = "---M----------------------------MMMM---------------M------------";
	$base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG";
	$base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG";
	$base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG";
}

#he Ciliate, Dasycladacean and Hexamita Nuclear Code
elsif($transl_table == 6){
	$aas    = "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
	$starts = "-----------------------------------M----------------------------";
	$base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG";
	$base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG";
	$base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG";
}

##Tables 7 and 8 have been deleted
#The Echinoderm and Flatworm Mitochondrial Code
elsif($transl_table == 9){
	$aas    = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG";
	$starts = "-----------------------------------M---------------M------------";
	$base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG";
	$base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG";
	$base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG";  
}

#The Euplotid Nuclear Code
elsif($transl_table == 10){
	$aas    = "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
	$starts = "-----------------------------------M----------------------------";
	$base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG";
	$base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG";
	$base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG";
}

#The Bacterial, Archaeal and Plant Plastid Code
elsif($transl_table == 11){
	$aas    = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
	$starts = "---M---------------M------------MMMM---------------M------------";
	$base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG";
	$base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG";
	$base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG";
}

#The Alternative Yeast Nuclear Code
elsif($transl_table == 12){
	$aas    = "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
	$starts = "-------------------M---------------M----------------------------";
	$base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG";
	$base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG";
	$base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG";
}

##Table 13-23 were not be added

foreach $i (1..(length $aas)) 
{
	$temp =substr($base1,($i-1),1);
	$temp.=substr($base2,($i-1),1);
	$temp.=substr($base3,($i-1),1);
	$hash{$temp}=substr($aas,($i-1),1);
}

foreach $i (1..(length $starts)) 
{
    if(substr($starts,($i-1),1) eq "M")
    {
	$temp =substr($base1,($i-1),1);
	$temp.=substr($base2,($i-1),1);
	$temp.=substr($base3,($i-1),1);
	$hashStarts{$temp}=substr($starts,($i-1),1);
    }
};

$/ = ">";
<IN>;
while(my $line = <IN>){	
    chomp $line;
    my @a = split(/\n/,$line);
    my $head = shift @a;
    print OUT ">$head\n";
    my $seq = join "",@a;
    $seq = uc($seq);

    if($head =~ /Lack\s*5\'\-end/ || $head =~ /Lack\s*both\s*end/){
	for($i=0;$i<length($seq);$i+=3)	{
	    if(exists $hash{substr($seq,$i,3)}){
		$out.=$hash{substr($seq,$i,3)};
	    }else{
		$out.="X";
	    }
	}
    }else{
	if(exists $hashStarts{substr($seq,0,3)}){
	    $out.=$hashStarts{substr($seq,0,3)};
	    for($i=3;$i<length($seq);$i+=3){
		if(exists $hash{substr($seq,$i,3)}){
		    $out.=$hash{substr($seq,$i,3)};
		}else{
		    $out.="X";
		}
	    }			
	}else{
	    for($i=0;$i<length($seq);$i+=3){
		if(exists $hash{substr($seq,$i,3)}){
		    $out.=$hash{substr($seq,$i,3)};
		}else{
		    $out.="X";
		}
	    }
	}

    };
    my @b = split(//,$out);
    if($b[-1] eq "*"){
	pop @b;
    };
    $out = join "",@b;
    print OUT "$out\n";
    $out = "";
}
