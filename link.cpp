#include<cstdio>
#include<iostream>
#include<algorithm>
#include<cstdlib>
#include<cmath>
#include<cctype>
#include<string>
#include<cstring>
#include<ctime>
#include<stack>
#include<queue>
#include<vector>
#include<set>
#include<map>
#include<unordered_map>
#include<fstream>
using namespace std;
typedef long long ll;
typedef double dd;
#define mx 20000007
#define maxlen 200007
#define inf 100000007
//#pragma comment(linker, "/STACK:102400000,102400000") 

//hyperparameter
const string datanum="3";  //1-4 

//parameter
const string file_pre="D:/wch/study/PhD/一上/生物信息学/pre/demo";
string short_file1=file_pre+"/data/short_1.fasta";
string short_file2=file_pre+"/data/short_2.fasta";
string contig=file_pre+"/output/dbg_contig.fasta";
string output_file=file_pre+"/output/contig.fasta";
vector<string> sh,ln,ans;
int lsize,ssize,sz;
int a[1007][5],tot;

//read_data
string get_reverse_complemantary(string s)
{
	int l=s.length();
	string ss="";
	char ch=' ';
	for (int i=l-1;i>=0;i--)
	{
		if (s[i]=='A') ch='T';
		else if (s[i]=='T') ch='A';
		else if (s[i]=='C') ch='G';
		else if (s[i]=='G') ch='C';
		else ch='N';
		ss=ss+ch;
	}
	return ss;
}

void read_short(string filepath)
{
	ifstream fin(filepath);
	string s;
	int t=0;
	while (getline(fin,s))
	{
		getline(fin,s);
		sh.push_back(s);
		s=get_reverse_complemantary(s);
		sh.push_back(s);
		t+=2;
	}
	cout<<"short "<<t<<endl;
	ssize=sh.size();
}

void read_contig(string filepath)
{
	ifstream fin(filepath);
	string s;
	int t=0;
	while (getline(fin,s))
	{
		getline(fin,s);
		ans.push_back(s);
		cout<<s.size()<<endl;
		t+=1;
	}
	cout<<"contig "<<t<<endl;
	sz=ans.size();
}

void read_data()
{
	read_short(short_file1);
	read_short(short_file2);
	read_contig(contig);
}

void link()
{
	int K=10;
	int p,ad;
	string a,b,c="";
	for (int ii=0;ii<sz;ii++)
	{
		if (ans[ii]=="") continue;
		ad=1;
		while (ad)
		{	
			ad=0;
			for (int jj=0;jj<ssize;jj++)
			{	
				a=ans[ii];
				b=sh[jj];
				int la=a.size(),lb=b.size();
				c=""; 
				for (int k=K;k<100;k++)
				{
					p=0;
					for (int i=la-k,j=0;j<k;i++,j++)
						if (a[i]==b[j]) p++; 
					if (p>=k) 
					{
						c=a.substr(0,la-k)+b;
						cout<<"link short: "<<ii<<' '<<a.size()<<"->"<<c.size()<<' '<<k<<endl;
						ad=1;
						break; 
					}			
				}
				for (int k=K;k<100;k++)
				{
					p=0;
					for (int i=0,j=lb-k;i<k;i++,j++)
						if (a[i]==b[j]) p++;
					if (p>=k) 
					{
						c=b.substr(0,lb-k)+a;
						cout<<"link short: "<<ii<<' '<<a.size()<<"->"<<c.size()<<' '<<k<<endl;
						ad=1;
						break;
					}
				}
				if (c!="") ans[ii]=c;
			}
		}
	}
}

bool cmpstring(string a,string b)
{
	return a.size()>b.size();
}

//output
void output()
{
	//sort(ln.begin(),ln.end(),cmpstring);
	string filepath=output_file;
	ofstream fout(filepath);
	for (int i=0;i<sz;i++)
	{
		if (ans[i]!="")
		{
			fout<<">contigs_"<<i<<endl;
			fout<<ans[i]<<endl;
			cout<<ans[i].size()<<endl;
		}
	}
}

int main()
{
	cout<<"Start!"<<endl<<endl;
	clock_t t1,t2;
	t1=clock();
	
	read_data();
	cout<<"Read data complete!"<<endl;
	t2=clock();
	cout<<"time: "<<(double)(t2-t1)/CLOCKS_PER_SEC<<"s"<<endl<<endl;
	
	link();
	cout<<"Link complete!"<<endl;

	output();
	cout<<"Finish!"<<endl;
	t2=clock();
	cout<<"time: "<<(double)(t2-t1)/CLOCKS_PER_SEC<<"s"<<endl<<endl;

	system("pause"); 
    return 0;
}

