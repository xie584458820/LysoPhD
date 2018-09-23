#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <algorithm>
#include <ctime>
#include <iostream>
#include<fstream>
using namespace std;

#define MAXSIZE 10000
#define AF 0x01
#define CF 0x02
#define GF 0x04
#define TF 0x08
#define FULLFLAG 0x0F
char strUnit[] = {'A', 'C', 'G', 'T'};
unsigned char hGeneFlag;
int maxrepeatlen=0,imebegin=0,imeend=0;
char* imeseq;
char maxrepeat[MAXSIZE / 2 + 1];
int **allrepeat=new int*[1000];
int num=0;


int MatchVerted(int MINMATCH, char *strQueue, int nQueueSize,char *contigpath ) {
ofstream outfile2,fout;
strcat(contigpath,"/VertedRepeat.fa");
outfile2.open(contigpath);
char *pQueue, strMatch[MAXSIZE / 2 + 1];
int i, nMatchLen, nBegin, nFlag, nPosition;
printf("\nVerted Repeat:\n");
for(nBegin = MINMATCH; nBegin < nQueueSize - MINMATCH; nBegin++) {
pQueue = strQueue + nBegin;
nMatchLen = 0;
nFlag = 0;
hGeneFlag = 0;
for(i = 0; i <= nQueueSize - nBegin; i++) {
if(strQueue[i] == pQueue[i]) {
if(nFlag == 0) {
nFlag = 1;
nPosition = i;
}
switch(strQueue[i]) {
case 'A':
hGeneFlag |= AF;
break;
case 'C':
hGeneFlag |= CF;
break;
case 'G':
hGeneFlag |= GF;
break;
case 'T':
hGeneFlag |= TF;
break;
}
strMatch[nMatchLen++] = strQueue[i];
if(nMatchLen >= nBegin && hGeneFlag == FULLFLAG) {
strMatch[nMatchLen] = 0;
int temp1=nPosition,temp2=nPosition + nBegin;
if(temp2-temp1>9300){
printf("Repeat:%s, Size: %d, Start Positioins:%d, %d\n", strMatch, nMatchLen, nPosition, nPosition + nBegin);
outfile2<<"Repeat: "<<strMatch<<", Size "<<nMatchLen<<", Start Positioins:"<<nPosition<<", "<<pQueue - strQueue + 1<<endl;
allrepeat[num]=new int[4];
allrepeat[num][0]=nMatchLen;
allrepeat[num][1]=nPosition;
allrepeat[num][2]=nPosition + nBegin;
allrepeat[num][3]=0;
num++;
if(nMatchLen>maxrepeatlen)
	{
		maxrepeatlen=nMatchLen;
		imebegin=nPosition;
		imeend=nPosition + nBegin;
		strcpy(maxrepeat,strMatch);
		}
}
nFlag = 0;
nPosition = 0;
nMatchLen = 0;
}
}
else {
if(nFlag == 1) {
if(nMatchLen >= MINMATCH && hGeneFlag == FULLFLAG) {
strMatch[nMatchLen] = 0;
int temp3=nPosition,temp4=nPosition + nBegin;
if(temp4-temp3>9300){
printf("Repeat:%s, Size: %d, Start Positioins:%d, %d\n", strMatch, nMatchLen, nPosition, nPosition + nBegin);
outfile2<<"Repeat: "<<strMatch<<", Size "<<nMatchLen<<", Start Positioins:"<<nPosition<<", "<<pQueue - strQueue + 1<<endl;
allrepeat[num]=new int[4];
allrepeat[num][0]=nMatchLen;
allrepeat[num][1]=nPosition;
allrepeat[num][2]=nPosition + nBegin;
allrepeat[num][3]=0;
num++;
if(nMatchLen>maxrepeatlen)
	{
		maxrepeatlen=nMatchLen;
		imebegin=nPosition;
		imeend=nPosition + nBegin;
		strcpy(maxrepeat,strMatch);
		}
}
}
nFlag = 0;
nPosition = 0;
nMatchLen = 0;
}
}
}
}
outfile2.close();
}


char* substring(char* ch,int pos,int length){
	char* pch=ch;

	char* subch=(char*)calloc(sizeof(char),length+1);
	int i;
	pch=pch+pos;
	for(i=0;i<length;i++){
		subch[i]=*(pch++);
		}	
	subch[length]='\0';
	return subch;
	}
	
int atoi(char s[]){
	int i,n,sign;
	for(i=0;isspace(s[i]);i++);
	sign=(s[i]=='-')?-1:1;
		if(s[i]=='+'||s[i]=='-')
			i++;
			for(n=0;isdigit(s[i]);i++)
				n=10*n+(s[i]-'0');
				return sign*n;
	}



bool cmp(int *p,int *q){
    if(p[0]==q[0]){
        if(p[1]==q[1]){
            if(p[2]==q[2]){
            		return p[3]<q[3];
            		}
            else return p[2]<q[2];
        }
        else return p[1]<q[1];
    }
    else return p[0]<q[0];
}


void sort(int **a)
{
		    for(int i=0;i<num;i++)
    {
        printf("%d\t%d\t%d\t%d\n",a[i][0],a[i][1],a[i][2],a[i][3]);
    }
		sort(a,a+num,cmp);
    cout<<"After sort"<<endl;
    for(int i=num-1;i>=0;i--){
	for(int j=i-1;j>=0;j--){
		if((abs(a[j][1]-a[i][1])<110)&&(abs(a[j][2]-a[i][2])<110)&&(i!=j)){
			a[i][3]=1;
			a[j][3]=1;
		}
		}
	}
    for(int i=num-1;i>=0;i--)
    {
        printf("%d\t%d\t%d\t%d\n",a[i][0],a[i][1],a[i][2],a[i][3]);
    }
}


int main(int argc,char** argv ) {
int minmatch=atoi(argv[2]);
char* contigpath=argv[4];
//Read the input sequence
FILE *pFile1=fopen(argv[1],"r");
fseek(pFile1,0,SEEK_END);
int len1=ftell(pFile1);
static char *pBuf1=new char[len1+1];
rewind(pFile1);
fread(pBuf1,1,len1,pFile1);
pBuf1[len1]=0;
//printf("The seq1 is:%s\n",pBuf1);
int nQueueSize=strlen(pBuf1);
//Find the repeat
MatchVerted(minmatch, pBuf1, nQueueSize, contigpath);
//Export the the found phage sequence将查找到的噬菌体序列输出
printf("The maxRepeat is :%s, Size: %d, Start Positioins:%d, %d\n", maxrepeat, maxrepeatlen, imebegin, imeend);
//Sort the repeats in the ascending order and output them in the descending order
sort(allrepeat);
imeseq=substring(pBuf1,imebegin,imeend-imebegin);
//printf("The imeseq is: %s\n,length is: %d\n",imeseq,strlen(imeseq));

ofstream outfile,fout;
outfile.open(argv[3]);
int n;
for(int i=num-1;i>=0;i--){
	outfile<<allrepeat[i][0]<<" "<<allrepeat[i][1]<<" "<<allrepeat[i][2]<<" "<<allrepeat[i][3]<<endl;
	n++;
	if(n>6) break;
	}
outfile.close();
fclose(pFile1);
return 0;
}