
use strict;
use warnings;
use Cwd;

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


opendir(DIRin, "$ARGV[0]") or die "Cannot open:$!";
mkdir $ARGV[1];
opendir(DIRout, "$ARGV[1]") or die "Cannot open:$!";

our $temp=$ARGV[2];
system("mkdir $temp");

#Read the location of programs
#our ($path_programm) = "/data/zhangxianglilan/prophage";
our ($path_programm) = ($0=~/^(.*)prophageIden.*$/); 
if($path_programm eq ""){ 
$path_programm= getcwd; }
else{
	
	}
print $path_programm."\n";


while (our $name = readdir DIRin){
    next if $name =~ /^\./;
	  print "$name\n";
	  system ("$path_glimmer/long-orfs -z 11 -n -t 1.15 -l $ARGV[0]/$name $temp/$name.longorfs");
	  system ("$path_glimmer/extract -t $ARGV[0]/$name $temp/$name.longorfs > $temp/$name.train");
	  system ("$path_glimmer/build-icm -r $temp/$name.icm < $temp/$name.train");
	  system ("$path_glimmer/glimmer3 -z 11 -o 50 -g 110 -t 30 -l $ARGV[0]/$name $temp/$name.icm $temp/$name");
	

	  system ("perl $path_programm/predict2formal.pl $temp/$name.predict $temp/$name.formal");
	  system ("$path_glimmer/multi-extract -t $ARGV[0]/$name $temp/$name.formal > $ARGV[1]/$name.cds");

	  #system ("rm -rf $temp/*");

	  system ("perl $path_programm/fa2std.pl $ARGV[1]/$name.cds $temp/$name.std");
	  system ("perl $path_programm/cds2aa.v2.pl $temp/$name.std $ARGV[1]/$name.aa 11");

	  print "Blasting...\n";
	  system ("$path_blast/blastp -query $ARGV[1]/$name.aa -db $phage_prot_db -num_alignments 5 -evalue 1e-10 -num_threads 16 -out $temp/$name.blastp");

	  system ("perl $path_programm/high-positives.pl $temp/$name.blastp $temp/$name.txt");
}

closedir DIRin;
closedir DIRout;