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
const int limit=80; 

//parameter
const string file_pre="D:/wch/study/PhD/一上/生物信息学/pre/demo";
string short_file1=file_pre+"/data/short_1.fasta";
string short_file2=file_pre+"/data/short_2.fasta";
string long_file=file_pre+"/data/long.fasta";
string long_fix=file_pre+"/output/long_fix.fasta";
vector<string> sh,ln,ans;
int a[2007][5],tot,lsize,ssize;

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
}

void read_long(string filepath)
{
	ifstream fin(filepath);
	string s;
	int t=0;
	while (getline(fin,s))
	{
		getline(fin,s);
		ln.push_back(s);
		t+=1;
	}
	cout<<"long "<<t<<endl;
}

void load_fix(string filepath)
{
	ifstream fin(filepath);
	string s;
	while (getline(fin,s))
	{
		getline(fin,s);
		ans.push_back(s);
	}
}

void read_data()
{
	read_short(short_file1);
	read_short(short_file2);
	read_long(long_file);
	load_fix(long_fix);	
}

void correct()
{
	ofstream fout(long_fix);
	lsize=ln.size(),ssize=sh.size();
	int tmp,id,ls,lt;
	string s,t;
	int sz=ans.size();
	for (int ii=0;ii<sz;ii++)
	{
		fout<<">contigs_"<<ii<<endl;
		fout<<ans[ii]<<endl;
	}
	cout<<sz<<" strings have been fixed before"<<endl;
	cout<<endl;
	for (int ii=sz;ii<lsize;ii++)
	{
		memset(a,0,sizeof(a));
		s=ln[ii];
		ls=s.size();
		for (int i=0;i<ls;i++)
		{
			if (s[i]=='A') a[i][0]+=1;
			if (s[i]=='T') a[i][1]+=1;
			if (s[i]=='C') a[i][2]+=1;
			if (s[i]=='G') a[i][3]+=1;
		}
		for (int jj=0;jj<ssize;jj++)
		{
			t=sh[jj];
			lt=t.size();
			for (int i=0;i<=ls-lt+1;i++)
			{
				tmp=lt;
				for (int j=0;j<lt;j++)
				{
					if (s[i+j]!=t[j]) tmp--;
					if (tmp<limit) break;
				}
				if (tmp>=limit)
				{
					for (int j=0;j<lt;j++)
					{
						if (t[j]=='A') a[i+j][0]+=1;
						if (t[j]=='T') a[i+j][1]+=1;
						if (t[j]=='C') a[i+j][2]+=1;
						if (t[j]=='G') a[i+j][3]+=1;
					}
				}
			}
		}
		s="";
		for (int i=0;i<ls;i++)
		{
			tmp=0,id=0;
			for (int j=0;j<4;j++)
			{
				if (a[i][j]>tmp) tmp=a[i][j],id=j;
			}
			if (id==0) s=s+'A';
			if (id==1) s=s+'T';
			if (id==2) s=s+'C';
			if (id==3) s=s+'G';		
		}
		fout<<">contigs_"<<ii<<endl;
		fout<<s<<endl;
		tmp=0;
		for (int i=0;i<ls;i++)
			if (ln[ii][i]!=s[i]) tmp++;
		cout<<"long string "<<ii<<" fixed"<<endl;
		cout<<tmp<<"/"<<ls<<" letters were corrected"<<endl;
		tot=tot+tmp; 
		cout<<endl;
	}	
	cout<<1.0*tot/lsize<<endl;
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

	correct();
	cout<<"Correct data complete!"<<endl;
	t2=clock();
	cout<<"time: "<<(double)(t2-t1)/CLOCKS_PER_SEC<<"s"<<endl<<endl;
		
	system("pause"); 
    return 0;
}

