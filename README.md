# LysoPhD
LysoPhD is a perl/c tool to automatically and accurately predict lysogenic phage sequences within bacterial sequencing data.

Installation

User can download all programs and put them into the same file path.
Before user run program, running environment should be configured.
Specially,the supporting software,glimmer, blast, newbler, FastQC and trimmomatic should be contained in the running computer. The paths of glimmer, blast and newbler should be written in the program "main.pl" and "prophageIden.pl". The path of FastQC and trimmomatic should be written in the program "Quality.py".
The position of path is behind the related annotation in the program.
Addtionally, the protein database of phages and nucleic acid database of bacteria are needed. The paths of this two database should be written in the programs "main.pl" and "prophageIden.pl".
The position of path is behind the related annotation in the program.


Usage:

perl main.pl [sequencing_data_1.fastq.gz] [sequencing_data_2.fastq.gz]

[sequencing_data_1.fastq.gz] = the input file of first sequencing data 
[sequencing_data_2.fastq.gz] = the input file of second sequencing data

Notes:
An example:

perl tools/main.pl 1339_S3_L001_R1_001.fastq.gz 1339_S3_L001_R2_001.fastq.gz

Who to blame: Xiangcheng Xie (xiexiangcheng1(at)163(dot)com)
