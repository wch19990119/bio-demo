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
#define inf 100000007
//#pragma comment(linker, "/STACK:102400000,102400000") 

//hyperparameter
const int K=25;
const int contig_num=10;   //输出的路径数量 

//parameter
const string file_pre="D:/wch/study/PhD/一上/生物信息学/pre/demo";
string short_file1=file_pre+"/data/short_1.fasta";
string short_file2=file_pre+"/data/short_2.fasta";
string long_file=file_pre+"/output/long_fix.fasta";
string contig=file_pre+"/output/dbg_contig.fasta";
vector<string> str,label,ans;
vector<int> path,cur_node,cur_edge,max_edge;
int n,m,cnt,tot,maxlength;
unordered_map<string,int> mp;
int valid[mx];

int last[mx];
struct edge
{
	int nxt;  
	int b;  
	int c;  //出现次数 
	int p;  //剩余可访问次数 
}e[mx];
struct node
{
	int ind;  //入度 
	int outd;  //出度
	int num;  //节点编号 
	int c;  //出现次数 
	int l;  //最远距离 
	int t;  //走向 
	int p;  //是否存在 
	int v;  //是否被访问过 
}nd[mx];


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

void read_file(string filepath)
{
	ifstream fin(filepath);
	string s;
	while (getline(fin,s))
	{
		getline(fin,s);
		str.push_back(s);
		s=get_reverse_complemantary(s);
		str.push_back(s);
	}
}

void read_data()
{
	read_file(short_file1);
	read_file(short_file2);
	read_file(long_file);	
}

//build graph
void add(int x,int y)
{
	for (int i=last[x];i;i=e[i].nxt)
	{
		if (e[i].b==y)
		{
			e[i].c++;
			return;
		}
	}
	nd[x].outd++;
	nd[y].ind++;
	m++;
	e[m].nxt=last[x];
	last[x]=m;
	e[m].b=y;
	e[m].c=1;
	e[m].p=1;
}

void build_graph()
{
	n=m=0;
	label.clear();
	label.push_back("0");
	int sz=str.size(),l=0,lst,now;
	string s,kmer;
	for (int i=0;i<sz;i++)
	{
		s=str[i];
		l=s.length();
		for (int j=0;j<=l-K;j++)
		{
			kmer=s.substr(j,K);
			if (mp.find(kmer)==mp.end())
			{
				n++;
				nd[n].num=n;
				nd[n].c=0;
				nd[n].p=1;
				mp[kmer]=n;
				label.push_back(kmer);
			}
			now=mp[kmer];
			nd[now].c++;
			if (j>0) add(lst,now);
			lst=now;
		}
	} 
}

//get contig
bool cmp(node a,node b)
{
	return a.ind<b.ind;
}

void deal_path()
{
	int sz=path.size();
	string s="";
	s=s+label[path[0]];
	for (int i=1;i<sz;i++)
	{
		s=s+label[path[i]].back();		
	}
	ans.push_back(s);
	cout<<s.size()<<endl;	
}

void dfs(int x)
{
	nd[x].v=1;
	int tmp=0,tt=0,y=0;
	for (int j=last[x];j;j=e[j].nxt)
	{
		y=e[j].b;
		if (nd[y].p)
		{
			if (!nd[y].v) dfs(y);
			if (tmp<nd[y].l)
			{
				tmp=nd[y].l;
				tt=y;
			}
			else if (tmp==nd[y].l&&nd[y].c>nd[tt].c) tt=y;
		}
	}
	nd[x].l=tmp+1;
	nd[x].t=tt;
}

void get_contig()
{		
	for (int ii=1;ii<=contig_num;ii++)
	//while (1)
	{
		for (int i=1;i<=n;i++) 
			nd[i].v=0,nd[i].l=0,nd[i].t=0;
		for (int i=1;i<=n;i++)
		{
			if (nd[i].p&&!nd[i].v) dfs(i);
		}
		int id=0;
		for (int i=1;i<=n;i++)
		{
			if (nd[i].p&&nd[i].l>nd[id].l) id=i;
		}
		if (id==0) break;
		path.clear();
		while (id)
		{
			if (nd[id].p==0) break;
			path.push_back(id);
			nd[id].p=0;
			id=nd[id].t;
		}
		deal_path();
	}
}

//contig_simplification
bool cmpstring(string a,string b)
{
	return a.size()>b.size();
}

//output
void output()
{
	string filepath=contig;
	ofstream fout(filepath);
	int sz=ans.size();
	for (int i=0;i<sz;i++)
	{
		{
			fout<<">contigs_"<<i<<endl;
			fout<<ans[i]<<endl;
			tot++;
			//if (i<20) cout<<ans[i].size()<<endl;
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
		
	build_graph();
	cout<<"Build graph complete!"<<endl;
	cout<<"This graph has "<<n<<" nodes and "<<m<<" edges"<<endl;
	t2=clock();
	cout<<"time: "<<(double)(t2-t1)/CLOCKS_PER_SEC<<"s"<<endl<<endl;
	
	get_contig();
	cout<<"Get config complete!"<<endl;
	cout<<ans.size()<<" strings have been found!"<<endl;
	t2=clock();
	cout<<"time: "<<(double)(t2-t1)/CLOCKS_PER_SEC<<"s"<<endl<<endl;
	
	output();
	cout<<tot<<" strings total!"<<endl;
	cout<<"Finish!"<<endl;
	t2=clock();
	cout<<"time: "<<(double)(t2-t1)/CLOCKS_PER_SEC<<"s"<<endl<<endl;
		
	system("pause"); 
    return 0;
}

