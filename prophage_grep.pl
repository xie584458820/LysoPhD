use strict;                                                                                         
use warnings;
use Cwd;
use Bio::SeqIO;
use Bio::Seq;
use threads;
print "111\n";
#use Bio::AlignIO;
   #open ctg file
our $dir=getcwd();
our $thread;
our $numpath=$ARGV[0];
our $possiblepath=$ARGV[1];
our $matrixpath=$ARGV[2];
our $readspath=$ARGV[3];
our $N=$ARGV[4];
our $phagepath=$ARGV[5];
our $name=$ARGV[6];
our $readsnum;
our $stringtemp;
our @reads1;
our @reads2;
our $title;
our @grepreads1;
our $grepreadsnum1=0;
our @grepreads2;
our $grepreadsnum2=0;
our @grepreads3;
our $grepreadsnum3=0;
our @grepreads4;
our $grepreadsnum4=0;
our $ingenome;

our ($path_programm) = ($0=~/^(.*)\/prophage_grep.*$/); 


system("makeblastdb -in $readspath/cat2.fa -dbtype nucl -parse_seqids -out $readspath/read");


print "The current pwd is: $dir\n";
print "The outputfile name is: $name\n";

our $i=0;
our $j=0;
our $findnum=0;
our $phagecount=0;
our $iffind=0;
our $vertedfind=0;
our $findlast=0;

my $k=0;
my $h=0;


sub startposition{
	my ($bigseq,$smallseq)=@_;
	my $len;
	while($bigseq=~m/$smallseq/g){
		$len=length($smallseq);
		my $end=pos($bigseq);
		my $start=$end-$len+1;
		return($start);
		}
	}

sub endposition{
	my ($bigseq,$smallseq)=@_;
	my $len;
	while($bigseq=~m/$smallseq/g){
		$len=length($smallseq);
		my $end=pos($bigseq);
		return($end);
		}
	}

	
open(NUM,$numpath) or die "failed to open:$!";
while(<NUM>){
	print "222\n";
	print $_;
	my $realctg;
	my $ctg;
	my $phgnum;
	my $line=$_;
	if(/contig/){				
	my @temp = split(/_/,$_);
	$ctg=$temp[0];
	$phgnum=$temp[1];
	print "contig: $ctg\n"; 
	if($phgnum==0){
		if(scalar(threads->list())<50){
		threads->new(\&allcontigthread,$ctg);
		}
			foreach $thread(threads->list(threads::all)){
		if($thread->is_joinable()){
			$thread->join();
			print scalar(threads->list()),"\t$j\t",localtime(time),"\n";
			
			}
	}
	
								}
	 else{
	 	for($j=1;$j<=$phgnum;$j++){
				if(scalar(threads->list())<50){
		threads->new(\&partcontigthread,$j,$ctg);
		}
			foreach $thread(threads->list(threads::all)){
		if($thread->is_joinable())	{
			$thread->join();
			print scalar(threads->list()),"\t$j\t",localtime(time),"\n";
			
																}
			}

	 														}
			}
	}	
	}

close(NUM);

	foreach $thread(threads->list(threads::all)){
		$thread->join();
		print scalar(threads->list()),"\t$j\t",localtime(time),"\n";
		}



sub allcontigthread(){
		
		my ($ctg)=@_;
		my @ctgtemp=split(/>/,$ctg);
		my $realctg=$ctgtemp[1];
		my $section=$realctg;
		print "This is allcontig\n";
		print "The realctg is $realctg\n";
		my $seqpath=$possiblepath."/".$realctg.".fa";
		my $matrixnamepath=$matrixpath."/".$realctg.".txt";
		my $phagenamepath=$phagepath."/".$realctg;
		my $tempdir=$possiblepath."/".$realctg."temp";
		print "seqpath is $seqpath\n";
		print "matrixnamepath is $matrixnamepath\n";
		open(DATA1,"<$seqpath");
		my @text=<DATA1>;
		chomp(@text);
		my $string=join " ",@text;
		my $seqlen=length($string);
		close(DATA1);
		$iffind=0;
		&allcontiggrep($string,$seqlen,$phagenamepath,$realctg,$section,$seqpath);
		if($iffind==0){
		&partcontiggrep($string,$seqlen,$phagenamepath,$matrixnamepath,$realctg,$section,$seqpath);
			}
	}

sub partcontigthread(){
		my ($j,$ctg)=@_;
		my $section=$ctg."_".$j;
		my $seqpath=$possiblepath."/".$ctg."_".$j.".fa";
	  print "This is partcontig\n";
	  my $matrixnamepath=$matrixpath."/".$ctg."_".$j.".txt";
		my $phagenamepath=$phagepath."/".$ctg."_".$j;
		my $tempdir=$possiblepath."/".$ctg."temp";
		print "seqpath is $seqpath\n";
		print "matrixnamepath is $matrixnamepath\n";
		open(DATA1,"<$seqpath");
		my @text=<DATA1>;
		chomp(@text);
		my $string=join " ",@text;
		my $seqlen=length($string);
		close(DATA1);
	 	&partcontiggrep($string,$seqlen,$phagenamepath,$matrixnamepath,$ctg,$section,$seqpath);
	 		}






sub allcontiggrep{
	 my ($string,$seqlen,$phagenamepath,$ctg,$section)=@_;
	 $stringtemp=$string;
	 $vertedfind=0;
	 $findlast=0;
	 print "Verted greping...\n";
	 my $tempseq=substr($string,0,length($string));
	 my $checkseq=substr($string,length($string)-10000,9999);
	 &grepall("all",$tempseq,$stringtemp,$phagenamepath,$ctg,$section,0,0,15);
	 if($findlast==1) {$iffind=1;}   
	}
	
	
	

sub partcontiggrep{
	my ($string,$seqlen,$phagenamepath,$matrixnamepath,$ctg,$section)=@_;
	$stringtemp=$string;
	open(DATA4,"<$matrixnamepath");
	system("mkdir $dir/result");
	system("rm -rf $dir/result/$section") if(-e "$dir/result/$section");	
while(<DATA4>){
	chomp;
	$title=$_;
	my @temp=split(/\s+/,$title);
	my $repeatlen=$temp[0];
		print "repeatlen: $repeatlen\n";
	my $start1=$temp[1];
	my $start2=$temp[2];
	my $mysta;
	my $myend;
	my $mystarevcom;
	my $myendrevcom;
	my $tt;
	print "@temp\n";
	$vertedfind=0;
	$findlast=0;
	
	my $repeat1=substr($string,$start1,$repeatlen);
	my $repeat2=substr($string,$start2,$repeatlen);
	if($temp[3]){

	 my $tempseq1=substr($string,$start1,$start2-$start1); 
	 my $checkseqfor1=substr($string,$start1-5000,10000);
	 my $checkseqback1=substr($string,$start2-5000,10000);
	 my $checkseq1=$checkseqfor1.$checkseqback1;
	 print "temp[3]=1\n";
	 &grep1(1,$tempseq1,$checkseq1,$stringtemp,$phagenamepath,$ctg,$section,$start1,$start2,$repeatlen,$repeat1,$repeat2);
	 if($findlast==1) {last;}
		
		}
		
		
if($temp[3]==0){
	 our $tempseq5=substr($string,$start1,$start2-$start1); 
	 my $checkseqfor2=substr($string,$start1-5000,10000);
	 my $checkseqback2=substr($string,$start2-5000,10000);
	 my $checkseq2=$checkseqfor2.$checkseqback2;
	 print "Verted greping...\n";
	 &grep1(0,$tempseq5,$checkseq2,$stringtemp,$phagenamepath,$ctg,$section,$start1,$start2,$repeatlen,$repeat1,$repeat2);
	 if($findlast==1) {last;}
				}
	
	}
	
	
	
	if($findnum==0){
		}
	
	close(DATA4);
	}
   

sub revcom{
	my ($a)=@_;
	my $recom=reverse $a;
  $recom=~tr/ACGTacgt/TGCAtgca/;
	return ($recom);
	}

sub grepall{
	 my ($case,$tempseq1,$stringtemp,$phagenamepath,$ctg,$section,$start1,$start2,$repeatlen)=@_;
	 my $i=0;
	 my $j=0;
	 my $position1;
	 my $position2;
	 my $dist;
	 my $find=0;
	 my $read1temp;
	 my $read2temp;
	 my $read1temprevcom;
	 my $read2temprevcom;
	 my @temp;
	 my @temp1;
	 my @tempid1;
	 my @temp2;
	 my @tempid2;
	 my $id;
	 my @idtemp;
	 my $flag; 
	 my @hit;
	 my @cohit;
	 my $hitnum=0;
	 my $cohitnum=0;
	 my @align1;
	 my @align2;
	 my $alignnum1=0;
	 my $alignnum2=0;
	 my $refseq;
	 my $p1;
	 my $p2;
	 my $t;
	 my $forrev;
	 my $hitnumtemp;
	 my $left;
	 my $right;
	 my $leftcondition;
	 my $rightcondition;
	 my $longleftcondition;
	 my $longrightcondition;
	 my $loopnum=0;
	 my $islong=0;
	 my $allphage;
	 my $size;
	 my @firstlinearray;
	 
	 
	 print "$ctg\n";
   system("mkdir $dir/result/$section");
	 open(RESULT,">$dir/result/$section/result.txt");
	 print RESULT "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	 my $mysta= substr($tempseq1,0,$N);     
	 my $myend= substr($tempseq1,length($tempseq1)-$N,$N);  
	 my $len=length($tempseq1);
	 $refseq=$myend.$mysta;
	 print "The N is $N\n";
	 print "The length is $len\n";

	 print RESULT "If case $case, mysta is: $mysta myend is: $myend\n";
	 
	 open(REF,">$dir/result/$section/ref.fa");	
	 open(COHIT,">$dir/result/$section/cohit.fa");	
	 print REF ">ref.fa\n";
	 print REF "$refseq";
	 system("blastn -query $dir/result/$section/ref.fa -db $readspath/read -out $dir/result/$section/ref.table -num_alignments 30000 -num_threads 10 -outfmt \"7 qseqid pident slen length mismatch gapopen qstart qend sstart send sseqid\"");
	 system("blastn -query $dir/result/$section/ref.fa -db $readspath/read -out $dir/result/$section/ref_all.table -num_alignments 30000 -num_threads 10");
	 system("perl $path_programm/handleiden.pl $dir/result/$section/ref.table $dir/result/$section/ref99.table 1");
	 close(REF);
	 

	 open(TABLE,"<$dir/result/$section/ref99.table");
	 open(GRE,">>$dir/result/$section/grep.txt");
	 while(<TABLE>){
	 	chomp;
	 	print RESULT "The line is: $_ ";
	 	@temp=split(/\s+/,$_);
	 	$p1=$temp[6];
	 	$p2=$temp[7];

	 	if($p1>$p2){
	 		$t=$p2;
	 		$p2=$p1;
	 		$p1=$t;
	 		}
	 	@idtemp=split(/&/,$temp[10]);
	 	$id=$idtemp[0];
	 	$forrev=$idtemp[1];
	 	$hit[$hitnum][0]=$id;
	 	$hit[$hitnum][1]=$forrev;
	 	$hit[$hitnum][2]=$p1;
	 	$hit[$hitnum][3]=$p2;
	 	
	 	print RESULT "!!!  $hit[$hitnum][0]  $hit[$hitnum][2]  $hit[$hitnum][3]\n";
		$hitnum++;
	 	}
	 	
		print RESULT "hitnum is $hitnum\n";
		print RESULT "...........\n";
		for($i=0;$i<$hitnum;$i++){
			print RESULT "µÚ$iÂÖ\n";
			for($j=$i+1;$j<$hitnum;$j++){

				
			if($hit[$i][0] eq $hit[$j][0]){
				
				$left=&min($hit[$i][2],$hit[$j][2]);
				$right=&max($hit[$i][3],$hit[$j][3]);
				$hit[$i][2]=$left;
				$hit[$j][2]=$left;
				$hit[$i][3]=$right;
				$hit[$j][3]=$right;
				

				
				}

				}
			}
		for($i=0;$i<$hitnum;$i++){
			print COHIT "$hit[$i][0]\t$hit[$i][1]\t$hit[$i][2]\t$hit[$i][3]\n";
			
			}
			
			
		print RESULT "hitnum is $hitnum\n";
		
		$leftcondition=$N-30;
		$rightcondition=$N+$repeatlen+30;
		$longleftcondition=$N-70;
		$longrightcondition=$N+$repeatlen+70;

		for($i=0;$i<$hitnum;$i++){
			if($hit[$i][2]<$leftcondition and $hit[$i][3]>$rightcondition){
				#print GRETEMP "all\n";
				print GRE "$hit[$i][0] $hit[$i][2] $hit[$i][3]\n";
				$find=1;
				$loopnum++;
				}

			if($hit[$i][2]<$longleftcondition and $hit[$i][3]>$longrightcondition){
				$islong=1;
				}
			
			}
	
	
	
	 close(TABLE);
	 close(GRE);
	 close(COHIT);
	 
	 	if(($find==1 and $loopnum>2) or ($find==1 and $islong==1 and $loopnum>1)){
	 		      my $prophageseq1=$tempseq1;
          	my $phagelength1=length($prophageseq1);
          	
            system("rm -rf $dir/result/$section/check") if(-e "$dir/result/$section/check");	
            system("rm -rf $dir/result/$section/check/in") if(-e "$dir/result/$section/check/in");	
            system("rm -rf $dir/result/$section/check/out") if(-e "$dir/result/$section/check/out");	
            system("rm -rf $dir/result/$section/check/temp") if(-e "$dir/result/$section/check/temp");	
            system("mkdir $dir/result/$section/check");
            system("mkdir $dir/result/$section/check/in");
            system("mkdir $dir/result/$section/check/out");
            system("mkdir $dir/result/$section/check/temp");
           
            open(CHE,">$dir/result/$section/check/phage.fa");
            print CHE ">$ctg  length=$phagelength1   numreads\n";
            print CHE "$prophageseq1";
            close(CHE);
            system("cp $dir/result/$section/check/phage.fa $dir/result/$section/check/in/phage.fa");
            system("perl $path_programm/prophageIden-20171030.pl $dir/result/$section/check/in $dir/result/$section/check/out $dir/result/$section/check/temp");
						system("python $path_programm/blast+2tab-nn+e+2covery.py $dir/result/$section/check/temp/phage.fa.blastp $dir/result/$section/check/temp/phage.fa.table 1");
						system("perl $path_programm/removefirstline.pl $dir/result/$section/check/temp/phage.fa.table $dir/result/$section/check/temp/phage.table");
						system("perl $path_programm/GetPosProtein.pl $dir/result/$section/check/temp/phage.table $dir/result/$section/check/temp/Possible2");
            
            my $checkpath="$dir/result/$section/check/temp/Possible2";
            if(-s $checkpath){
            
   
            system("perl $path_programm/ace_contig.pl $dir/Assembly/454Contigs.ace $dir/result/$section $ctg $dir/result/$section/out.txt");
            open(ACE,"$dir/result/$section/out.txt");
            $ingenome=<ACE>;
            close(ACE);
            if($ingenome=~/OnBacterialGenome/){
            $phagecount++;
            $phagenamepath=$phagenamepath."_$start1\_$start2\_newactiveprophage2_$phagecount.fa";
            $size=length($tempseq1);
            $allphage="/mount/user/xiexc/allphage/"."$name\_$size\_$phagecount\_a.fa";
            print "allphagepath: $allphage";
          
    				print "It is a real phage!\n";   print RESULT "It is a real phage!\n";
            print "The prophageseq1 is: \n";  print RESULT "The prophageseq1 is: \n";
            print "$phagenamepath\n$prophageseq1\n";   print RESULT "$phagenamepath\n$prophageseq1\n";
            print "The length of prophageseq1 $name is: $phagelength1\n";  print RESULT "The length of prophageseq1 $name is: $phagelength1\n";
 
            open(OUT,">$phagenamepath");
            open(ALL,">$allphage");
            open(AC,"$dir/result/$section/out.txt");
            print OUT "$dir\n";
            print OUT ">$ctg"."_phage2_".$phagecount."\n";  
            print OUT "$prophageseq1";  
            print ALL "$dir\n";
            print ALL ">$ctg"."_phage2_".$phagecount."\n";  
            print ALL "$prophageseq1";  
						
						print OUT "\n\n";
						print ALL "\n\n";
						while(<AC>){
							print OUT $_;
							print ALL $_;
							}          
            
            
            close(AC);
            $vertedfind=1;
            $findlast=1;
            }
						          
            close(ALL);
            close(OUT);
            }
            
						#system("rm -rf check/in");
						#system("rm -rf check/out");
						#system("rm -rf check/temp");
						#system("rm -rf check");
	 		
	 		}
	 	close(RESULT);
	 	}

sub grep1{
	 my ($case,$tempseq1,$checkseq,$stringtemp,$phagenamepath,$ctg,$section,$start1,$start2,$repeatlen,$repeat1,$repeat2)=@_;
	 my $i=0;
	 my $j=0;
	 my $position1;
	 my $position2;
	 my $dist;
	 my $find=0;
	 my $read1temp;
	 my $read2temp;
	 my $read1temprevcom;
	 my $read2temprevcom;
	 my @temp;
	 my @temp1;
	 my @tempid1;
	 my @temp2;
	 my @tempid2;
	 my $id;
	 my @idtemp;
	 my $flag; 
	 my @hit;
	 my @cohit;
	 my $hitnum=0;
	 my $cohitnum=0;
	 my @align1;
	 my @align2;
	 my $alignnum1=0;
	 my $alignnum2=0;
	 my $refseq;
	 my $p1;
	 my $p2;
	 my $t;
	 my $forrev;
	 my $hitnumtemp;
	 my $left;
	 my $right;
	 my $leftcondition;
	 my $rightcondition;
	 my $longleftcondition;
	 my $longrightcondition;
	 my $loopnum=0;
	 my $islong=0;
	 my $allphage;
	 my $size;
	  
	 
	 open(CCC,">$dir/result/$section/checkcontigseq.fa");
	#print "$dir/result/$section/checkcontigseq.fa\n";
	 print CCC "$checkseq\n";
	 close(CCC);
	 print "$ctg\n";
   system("mkdir $dir/result/$section");
	 open(RESULT,">$dir/result/$section/result.txt");
	 open(GRETEMP,">$dir/result/$section/greptemp.txt");
	 print RESULT "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	 my $mysta= substr($tempseq1,0,$N);     
	 my $myend= substr($tempseq1,length($tempseq1)-$N,$N);  
	 my $len=length($tempseq1);
	 $refseq=$myend.$mysta;
	 print "The N is $N\n";
	 print "The length is $len\n";
	 #my $mystarevcom=revcom($mysta);
	 #my $myendrevcom=revcom($myend);
	 print RESULT "If case $case, mysta is: $mysta myend is: $myend\n";
	 #print RESULT "mystarevcom is: $mystarevcom myendrevcom is: $myendrevcom\n";
	 
	 open(REF,">$dir/result/$section/ref.fa");	
	 open(COHIT,">$dir/result/$section/cohit.fa");	
	 print REF ">ref.fa\n";
	 print REF "$refseq";
	 print RESULT "checkcontigseq.fa: $dir/result/$section/checkcontigseq.fa\n";
	 print "ref.fa: $dir/result/$section/ref.fa\n";
	 system("blastn -query $dir/result/$section/ref.fa -db $readspath/read -out $dir/result/$section/ref.table -num_alignments 30000 -num_threads 10 -outfmt \"7 qseqid pident slen length mismatch gapopen qstart qend sstart send sseqid\"");
	 system("blastn -query $dir/result/$section/ref.fa -db $readspath/read -out $dir/result/$section/ref_all.table -num_alignments 30000 -num_threads 10");
	 system("perl $path_programm/handleiden.pl $dir/result/$section/ref.table $dir/result/$section/ref99.table 0");
	 close(REF);
	 
	 system("blastn -query $dir/result/$section/checkcontigseq.fa -db $readspath/read -out $dir/result/$section/contig.table -num_alignments 200000 -num_threads 10 -outfmt \"7 qseqid pident slen length mismatch gapopen qstart qend sstart send sseqid\"");
	 system("perl $path_programm/handleidencon.pl $dir/result/$section/contig.table $dir/result/$section/contig97.table");

	 open(TABLE,"<$dir/result/$section/ref99.table");


	 while(<TABLE>){
	 	chomp;
	 	print RESULT "The line is: $_ ";
	 	@temp=split(/\s+/,$_);
	 	$p1=$temp[6];
	 	$p2=$temp[7];

	 	if($p1>$p2){
	 		$t=$p2;
	 		$p2=$p1;
	 		$p1=$t;
	 		}
	 	@idtemp=split(/&/,$temp[10]);
	 	$id=$idtemp[0];
	 	$forrev=$idtemp[1];
	 	$hit[$hitnum][0]=$id;
	 	$hit[$hitnum][1]=$forrev;
	 	$hit[$hitnum][2]=$p1;
	 	$hit[$hitnum][3]=$p2;
	 	print RESULT "!!!  $hit[$hitnum][0]  $hit[$hitnum][2]  $hit[$hitnum][3]\n";
		$hitnum++;
	 	}
	 	
		print RESULT "hitnum is $hitnum\n";
		print RESULT "...........\n";
		for($i=0;$i<$hitnum;$i++){
			print RESULT "µÚ$iÂÖ\n";
			for($j=$i+1;$j<$hitnum;$j++){

				
			if($hit[$i][0] eq $hit[$j][0]){
				$left=&min($hit[$i][2],$hit[$j][2]);
				$right=&max($hit[$i][3],$hit[$j][3]);
				$hit[$i][2]=$left;
				$hit[$j][2]=$left;
				$hit[$i][3]=$right;
				$hit[$j][3]=$right;
				

				
				}

				}
			}
		for($i=0;$i<$hitnum;$i++){
			print COHIT "$hit[$i][0]\t$hit[$i][1]\t$hit[$i][2]\t$hit[$i][3]\n";
			
			}
			
			
		print RESULT "hitnum is $hitnum\n";
		
		$leftcondition=$N-30;
		$rightcondition=$N+$repeatlen+30;
		$longleftcondition=$N-70;
		$longrightcondition=$N+$repeatlen+70;

		for($i=0;$i<$hitnum;$i++){
			if($hit[$i][2]<$leftcondition and $hit[$i][3]>$rightcondition){
				#print GRETEMP "notall\n";
				print GRETEMP "$hit[$i][0] $hit[$i][1] $hit[$i][2] $hit[$i][3]\n";
				$find=1;
				$loopnum++;
				}

			
			}
	
	
	
	 close(TABLE);
	 close(GRETEMP);
	 close(COHIT);
	 

	 system("perl $path_programm/handlegrep.pl $dir/result/$section/greptemp.txt $dir/result/$section/contig97.table $dir/result/$section/greptemp1.txt $dir/result/$section/grep.txt");	 
	 

	 
	 my $filename1="$dir/result/$section/grep.txt";
	 	if(-s $filename1){
	 		      my $prophageseq1=$tempseq1;
          	my $phagelength1=length($prophageseq1);
          	
            system("rm -rf $dir/result/$section/check") if(-e "$dir/result/$section/check");	
            system("rm -rf $dir/result/$section/check/in") if(-e "$dir/result/$section/check/in");	
            system("rm -rf $dir/result/$section/check/out") if(-e "$dir/result/$section/check/out");	
            system("rm -rf $dir/result/$section/check/temp") if(-e "$dir/result/$section/check/temp");	
            system("mkdir $dir/result/$section/check");
            system("mkdir $dir/result/$section/check/in");
            system("mkdir $dir/result/$section/check/out");
            system("mkdir $dir/result/$section/check/temp");
           
            open(CHE,">$dir/result/$section/check/phage.fa");
            print CHE ">$ctg  length=$phagelength1   numreads\n";
            print CHE "$prophageseq1";
            close(CHE);
            system("cp $dir/result/$section/check/phage.fa $dir/result/$section/check/in/phage.fa");
            system("perl $path_programm/prophageIden-20171030.pl $dir/result/$section/check/in $dir/result/$section/check/out $dir/result/$section/check/temp");
						system("python $path_programm/blast+2tab-nn+e+2covery.py $dir/result/$section/check/temp/phage.fa.blastp $dir/result/$section/check/temp/phage.fa.table 1");
						system("perl $path_programm/removefirstline.pl $dir/result/$section/check/temp/phage.fa.table $dir/result/$section/check/temp/phage.table");
						system("perl $path_programm/GetPosProtein.pl $dir/result/$section/check/temp/phage.table $dir/result/$section/check/temp/Possible2");
            
            my $checkpath="$dir/result/$section/check/temp/Possible2";
            if(-s $checkpath){
            $phagecount++;
            $phagenamepath=$phagenamepath."_$start1\_$start2\_newactiveprophage2_$phagecount.fa";
            $size=$start2-$start1;
            $allphage="/mount/user/xiexc/allphage/"."$name\_$size\_$phagecount\_p.fa";
            print "allphagepath: $allphage";
            
    				print "It is a real phage!\n";   print RESULT "It is a real phage!\n";
            print "The prophageseq1 is: \n";  print RESULT "The prophageseq1 is: \n";
            print "$phagenamepath\n$prophageseq1\n";   print RESULT "$phagenamepath\n$prophageseq1\n";
            print "The length of prophageseq1 $name is: $phagelength1\n";  print RESULT "The length of prophageseq1 $name is: $phagelength1\n";
 
            open(OUT,">$phagenamepath");
            open(ALL,">$allphage");
            print OUT "$dir\n";
            print OUT ">$ctg"."_phage2_".$phagecount."\n";  
            print OUT "repeat1: $repeat1   repeat2: $repeat2\n";
            print OUT "$prophageseq1";  
            print ALL "$dir\n";
            print ALL ">$ctg"."_phage2_".$phagecount."\n";  
            print ALL "repeat1: $repeat1   repeat2: $repeat2\n";
            print ALL "$prophageseq1";  
            $vertedfind=1;
            $findlast=1;
            close(OUT);
            close(ALL);
            }
            
						#system("rm -rf check/in");
						#system("rm -rf check/out");
						#system("rm -rf check/temp");
						#system("rm -rf check");
	 		
	 		}
	 	close(RESULT);
	 	}
	
	sub max{
		my ($a,$b)=@_;
		my $max;
		if($a>=$b){ $max=$a;}
		if($b>$a){ $max=$b;}
		return ($max);
		}
		
	sub min{
		my ($a,$b)=@_;
		my $min;
		if($a<=$b){ $min=$a;}
		if($b<$a){ $min=$b;}
		return ($min);
		}	