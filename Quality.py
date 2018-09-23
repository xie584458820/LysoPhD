#!/usr/bin/python
# -*- coding: utf-8 -*-
#----------------------------------------------------------------------------------------------
#Module Manual
#----------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
#import modules
#----------------------------------------------------------------------------------------------
import sys
import argparse         #Parameters of script
import glob
import os
import time             #Time
import datetime
import os.path
import re               #Regular expression
import threading        #Multithreading
import linecache        
from bs4 import BeautifulSoup
#----------------------------------------------------------------------------------------------
#Argparse
#----------------------------------------------------------------------------------------------
Argparse=argparse.ArgumentParser(description = "The script is used to filter fastq file")
Argparse.add_argument('--type', type = str, default = 'SE', help = 'The type of sequencing sethod, SE or PE, default SE')
Argparse.add_argument('--infile1', type = str, required = True, help = 'The input fastq1 file')
Argparse.add_argument('--infile2', type = str, help = 'The input fastq2 file, if type is PE')
Argparse.add_argument('--outdir', type = str, required = True, help = 'The output dir')
args = Argparse.parse_args()
#----------------------------------------------------------------------------------------------
#define class
#----------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------
#Global 
#----------------------------------------------------------------------------------------------
current_dir = os.getcwd()
target_dir = os.path.abspath(args.outdir)

#----------------------------------------------------------------------------------------------
#define funtion
#----------------------------------------------------------------------------------------------
def adapter(html_dir):
    adapter = html_dir + '/adapter.fa'
    Outfile = open(adapter, 'w')                                    #Open file
    TempHtmls = html_dir + "/*.html"
    Htmls = glob.glob(TempHtmls)
    n = 1
    for Html in Htmls:
        soup = BeautifulSoup(open(Html, "r"), "lxml")               #Read the html file
        strings = soup.find(text="No overrepresented sequences")    #Determine whether there is adpter sequences
        if strings:                                                 #If there is no adapter sequence, then skip.
            continue
        else:
            table = soup.find_all('table')[1]                       
            trs = table.find_all('tr')                              
            for i in range(1, len(trs)):                            
                td = trs[i].find_all('td')[0]                       
                #Write the sequence into the file
                Outfile.write(">" + str(n) + "\n" + td.contents[0] + "\n")
                n = n + 1
def single_end():
    file1 = os.path.abspath(args.infile1)
    Command1 = '/mount/user/zhaoyc/softwares/Trim/FastQC/fastqc -o ' + target_dir + ' -t 80 ' + file1
    Command2 = 'java -jar /mount/user/zhaoyc/softwares/Trim/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 80 -phred33 -trimlog trimmomatic.log' + ' ' + file1 + ' ' + 'trimmomatic.fastq.gz ILLUMINACLIP:adapter.fa:2:30:10 HEADCROP:15 MINLEN:50 AVGQUAL:20'
    Command3 = '/mount/user/zhaoyc/softwares/Trim/FastQC/fastqc -o ' + target_dir + ' -t 80 trimmomatic.fastq.gz'
    os.chdir(target_dir)
    os.system(Command1)
    adapter(target_dir)
    os.system(Command2)
    os.system(Command3)
    os.chdir(current_dir)
def pair_end():
    file1 = os.path.abspath(args.infile1)
    file2 = os.path.abspath(args.infile2)
    tp1 = target_dir + '/trimmomatic_1P.fastq'
    tp2 = target_dir + '/trimmomatic_2P.fastq'
    tu1 = target_dir + '/trimmomatic_1U.fastq'
    tu2 = target_dir + '/trimmomatic_2U.fastq'   
    Command1 = '/mount/user/zhaoyc/softwares/Trim/FastQC/fastqc -o ' + target_dir + ' -t 80 ' + file1 + ' ' + file2
    Command2 = 'java -jar /mount/user/zhaoyc/softwares/Trim/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 20 -phred33' + ' ' + file1 + ' ' + file2 + ' ' + tp1 + ' ' + tu1 + ' ' + tp2 + ' ' + tu2 + ' ' + 'ILLUMINACLIP:' + target_dir + '/adapter.fa:2:30:10 SLIDINGWINDOW:3:20 MINLEN:35'
    Command3 = '/mount/user/zhaoyc/softwares/Trim/FastQC/fastqc -o ' + target_dir + ' -t 80 ' + 'trimmomatic_1P.fastq trimmomatic_2P.fastq'
    #os.chdir(target_dir)
    os.system(Command1)
    adapter(target_dir)
    os.system(Command2)
    os.system(Command3)
    #os.chdir(current_dir)
def main():
    if args.type == 'SE':
        single_end()
    elif args.type == 'PE':
        pair_end()
    else:
        print('Please input correct type')
        sys.exit(1)
#----------------------------------------------------------------------------------------------
#Main function
#----------------------------------------------------------------------------------------------
if "__main__" == __name__:
    main()