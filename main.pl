#usage: perl main.pl $infile1[0] $infile2[1]
use strict;
use warnings;
use threads;
use Cwd;

my $thread;

#Add the paths of supporting software
#glimmer_path
our $path_glimmer = '/mount/user/sunqiang/soft/glimmer3.02/bin';
#blast_path
our $path_blast = '/mount/user/fanhang/software/ncbi-blast-2.2.31+/bin';

#Add paths of databases
#Add a path of database of phage protein
our $phage_prot_db = '/mount/data/backup/biodb/phage-ncbi-201508-prot/phagedatabase-prot.fa';
#Add a path of nucleic acid database of bacteria
our $nt_Bacteria = '/mount/data/database/nt-Bac/nt-Bac';

our $Newbler_runAssemblypath='/opt/454/bin';

my $cwd=getcwd();
our $dir=getcwd;
my $file;
my @dir;
my $dirtemp;
my $flagposition;
my $taglength;
my $name=$ARGV[2];
print "The path is:$dir\n";


print "The CWD PATH= ",$cwd,"\n";
mkdir 'ok';

#Read the location of programs
our ($path_programm) = ($0=~/^(.*)\/main.*$/); 
print $path_programm."\n";
system("rm -rf phage") if(-e "phage");
system("rm -rf contig") if(-e "contig");
system("rm -rf Temp") if(-e "Temp");
system("rm -rf result") if(-e "result");
mkdir 'Temp';
mkdir 'phage';
mkdir 'contig';
system("python $path_programm/Quality.py --type PE --infile1 $ARGV[0] --infile2 $ARGV[1] --outdir Temp");
system("gunzip Temp/trimmomatic_1P.fastq.gz");
system("gunzip Temp/trimmomatic_2P.fastq.gz");
system("head -n 3000000 Temp/trimmomatic_1P.fastq > Temp/1.fastq");
system("head -n 3000000 Temp/trimmomatic_2P.fastq > Temp/2.fastq");

system("$Newbler_runAssemblypath/runAssembly -force -m -o Assembly -cpu 30 Temp/1.fastq Temp/2.fastq");
system("perl $path_programm/aceReader.pl Assembly/454AllContigs.ace Assembly/$name.ace");
system("perl $path_programm/fa2std.pl Assembly/454AllContigs.fna Assembly/454AllContigs.std");

system("python $path_programm/fq2fa.py Temp/trimmomatic_1P.fastq Temp/trimmomatic_1P.fa");
system("python $path_programm/fq2fa.py Temp/trimmomatic_2P.fastq Temp/trimmomatic_2P.fa");
system("head -n 500 Temp/trimmomatic_1P.fa > Temp/check1.fa");
system("head -n 500 Temp/trimmomatic_2P.fa > Temp/check2.fa");

system("perl $path_programm/handle.pl Temp/trimmomatic_1P.fa Temp/handle1.fa 1");
system("perl $path_programm/handle.pl Temp/trimmomatic_2P.fa Temp/handle2.fa 2");	
system("cat Temp/handle1.fa Temp/handle2.fa > Temp/cat2.fa");


open(DATA,"Assembly/454AllContigs.std");
our $ok=0;
our $contig;
mkdir 'contig'; 
our $tag;
our $title;
our $count=0;
our @allcontig;
our $len;
while(<DATA>){
	chomp;
	$title=$_;
	if(/>contig/){
		$tag=$title;
		my @temp=split(/\s+/,$title);
		my @temp1=split(/>/,$temp[0]);		
		$contig=$temp1[1];
		print "111\n";
		if(/length=(\w+)(\s+)/)
		{ 
			print "222\n";
			$len=$1;
			print "The len is $len\n";
		}
		if($len>10000){
			$ok=1;
			$allcontig[$count]=$contig;
			$count++;
			print "The contig $contig length is $len\n";
			system("mkdir contig/$contig");
			open(CONTIG,">contig/$contig/$contig.fa");
			print CONTIG "$tag\n";
			close CONTIG;
			}
		}
	else{
		if($ok==1){
			open(CONTIG,">>contig/$contig/$contig.fa");
			print CONTIG "$title\n";
			system("mkdir contig/$contig/in");
			system("mkdir contig/$contig/out");
			system("cp contig/$contig/$contig.fa contig/$contig/in/$contig.fa");
			close CONTIG;					
			$ok=0;
			}
		}
	}
close DATA;

our $i=0;
			system("mkdir contig/ctgphgNum");
			system("mkdir contig/PossibleSeq");
			system("mkdir contig/Matrix");
			system("mkdir phage");				
					
			my $part2path;
			
our $j=0;
print "allcontig: @allcontig\n";
print "count is: $count\n";
while(){
	last if($i>=$count);	
		while(scalar(threads->list())<20 and $i<$count){
		print "i:  $i\n";
		threads->new(\&part2,$i);
		$i++;
		}
	foreach $thread(threads->list(threads::all)){
		if($thread->is_joinable()){
			$thread->join();
			print scalar(threads->list()),"\t$i\t",localtime(time),"\n";			
			}
	}
}

	foreach $thread(threads->list(threads::all)){
		$thread->join();
		print scalar(threads->list()),"\n";
		}
		
	print localtime(time),"\n";

	system("cat contig/ctgphgNum/* > contig/ctgphgNum.txt");


	system("makeblastdb -in Temp/cat2.fa -dbtype nucl -parse_seqids -out Temp/read");
	system("perl $path_programm/prophage_grep.pl contig/ctgphgNum.txt contig/PossibleSeq contig/Matrix Temp 1000 phage $name");


sub part2(){
	my ($i)=@_;
	$contig=$allcontig[$i];
	print "This is i: $i\n";
	print "This is contig[$i]: $contig\n";
	my $part2path="contig/".$contig."/part2";
	system("perl $path_programm/prophageIden.pl contig/$contig/in contig/$contig/out contig/$contig/temp");
	system("mkdir contig/$contig/part2");	
	system("python $path_programm/blast+2tab-nn+e+2covery.py contig/$contig/temp/$contig.fa.blastp contig/$contig/temp/$contig.fa.table 1");
	system("perl $path_programm/removefirstline.pl contig/$contig/temp/$contig.fa.table $part2path/$contig.table");
	system("perl $path_programm/GetPos.pl contig/$contig/temp/$contig.fa.txt.title $part2path/Possible1");
	system("perl $path_programm/GetPosProtein.pl $part2path/$contig.table $part2path/Possible2");
	system("perl $path_programm/fa2std.pl contig/$contig/in/$contig.fa $part2path/$contig.std");
	system("cat $part2path/Possible1 $part2path/Possible2 > $part2path/catPossible");
	system("perl $path_programm/corange.pl $part2path/catPossible $part2path/coPossible");
	system("perl $path_programm/getprophageseqnew.pl $part2path/$contig.std $part2path/coPossible $part2path/prophage.fa");
	system("perl $path_programm/extractphage.pl $part2path/prophage.fa contig/ctgphgNum contig/PossibleSeq contig/Matrix");
	}
	
sub checkflag{
	my ($a,$b)=@_;
	my $flagposition;
  my $firstline1;
  my $firstline2;
  my $firstlinenum1;
  my $firstlinenum2;
  my @array1;
  my @array2;
  my $i;
  my $taglength;   
  open(CHECK1,$a);
  open(CHECK2,$b);
	$firstline1=<CHECK1>;
	$firstline2=<CHECK2>;
	@array1=split //,$firstline1;
	@array2=split //,$firstline2;
	$firstlinenum1=length($firstline1);
	$firstlinenum2=length($firstline2);
	$taglength=$firstlinenum1;
		print "The array1 is: @array1\n";
		print "The array2 is: @array2\n";
		print "The firstlinenum2 is: $firstlinenum2\n";
	for($i=0;$i<$firstlinenum1;$i++){
		print "The array1 array2 is: $array1[$i]  $array2[$i]\n";
			if(($array1[$i] eq "1")&&($array2[$i] eq "2")){
		print "The final array1 array2 is: $array1[$i]  $array2[$i]\n";
				$flagposition=$i;
				}
		}
	print "The inflagpostion is: $flagposition\n";
	return ($flagposition,$taglength);
	}