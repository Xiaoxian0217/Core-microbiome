/*
	Reference:
	1.Ulrik Brandes--A Faster Algorithm for Betweenness Centrality
	2.Hirokazu Toju--Core microbiomes for sustainable agroecosystems	
	Version1:Designed for Calculation of  Core Microbiome
	Date:2018/09/24
*/
#include <iostream>
#include <fstream>
#include <cstring>
#include<string>
#include <sstream>
#include <stack>
#include <queue>
#include <vector>
#include<stdlib.h>
#include<algorithm>
using namespace std;

//使用邻接矩阵完成图的遍历
struct Graph_array {
    int vexnum;  //图的顶点数
    int edge;    //图的边数
    int ** arc;  //邻接矩阵
	int ** ns;  //共有节点数
	int ** cij;  //关联系数
    int kind;    //0,为有向图，1，为无向图
	double * weight; //表示每个vertex的权重
	int * nlink; //链接点数目

    string * infromation; //表示每个vertex的信息
};

//使用邻接矩阵表示的图
void createGraph_by_array(int **edge, Graph_array & g) 
{
    int i = 0;
    g.arc = new int*[g.vexnum];//为邻接矩阵开辟空间
	g.ns = new int*[g.vexnum];//为ns矩阵开辟空间
	g.cij = new int*[g.vexnum];//为cij矩阵开辟空间
    for (i = 0; i < g.vexnum; i++)
    {
        g.arc[i] = new int[g.vexnum];
		g.ns[i] = new int[g.vexnum];
		g.cij[i] = new int[g.vexnum];
        for (int j = 0; j < g.vexnum; j++)
		{
			g.arc[i][j] = 0;
			g.ns[i][j] = 0;
			g.cij[i][j] = 0;
		}
    }
    for (i = 0; i < g.edge; i++)
    {
        //对矩阵进行赋值
		//以对角线为零元素的对称矩阵
        g.arc[edge[i][0]][edge[i][1]] = 1;
		g.arc[edge[i][1]][edge[i][0]] = 1;
    }
	int it;
	//统计Ri	
	for (int iv = 0; iv< g.vexnum; iv++)
    {
		it=0;	
		for(int j = 0; j < g.vexnum; j++)
		{							
			if(g.arc[iv][j] > 0)
			{
				it=it+1;
			}			
		}
		g.nlink[iv] = it;
    }
	//计算共有节点数S
	int ic,ir;
	 for (i = 0; i < g.edge; i++)
    {   
		it=0;
		ir=edge[i][0];
		ic=edge[i][1];
		for(int im=0; im < g.vexnum; im++)
		{
			if(g.arc[ir][im]*g.arc[ic][im] > 0)
			{
				it=it+1;
			}			
		}
		g.ns[ir][ic] = it;
		g.ns[ic][ir] = it;
    }
	 //计算交互作用系数Cij
	 for (i = 0; i < g.edge; i++)
    {   
		ir=edge[i][0];
		ic=edge[i][1];
		g.cij[ir][ic]=(g.nlink[ir]-g.ns[ir][ic])*(g.nlink[ic]-g.ns[ir][ic]);		
    }
}

//打印邻接矩阵
void print_array(Graph_array g) {

    int i = 0;
	ofstream os("Array.dat");
    for (i = 0; i <g.vexnum; i++)
	{
        for (int j = 0; j < g.vexnum; j++)
		{
            os << g.arc[i][j] << " ";
        }
        os << endl;
    }
}

//Ulrik Brandes Algorithm
void FBCA(int **edge , Graph_array g )
{
	double* CB=new double[g.vexnum] ;
	double* Cij=new double[g.edge] ;
	int n=g.vexnum;
    for(int i = 0 ;i < n; i++)
	{
		CB[i] = 0;
	}
	for(int i = 0 ;i < g.edge; i++)
	{
		Cij[i] = 0;
	}
	for(int s = 0; s < n; s++)
	{
	  //用于存储最短路径到达节点在最短路径上的前一个节点
	  //一个链表数组
	  vector< vector<int> > p(n); 
	  stack<int> S;
	  queue<int> Q;
    
	  double* a=new double[n];
   	  for(int h = 0 ; h < n; h++)
	  {
		 //从s出发，到达顶点t的最短路径数目
		a[h] = 0.0;
	  }
	  a[s] = 1.0;
	  //从s出发，到达t的路径的长度
	  int* b=new int[n];
	  for(int e = 0 ; e < n; e++)
	  {
		b[e] = -1;
	  }
	  b[s] = 0;

	  Q.push(s);

	  while(!Q.empty())
	  {
		  int v = Q.front();
		  Q.pop();
		  S.push(v);

		  for(int w = 0; w < n ;w++)
			  if(g.arc[v][w]!=0)
			  {
				  if(b[w] < 0)
				  {
					  Q.push(w);
					  b[w] = b[v] +1;
				  }

				  //shortest path to w via v
				   
				  if(b[w] == b[v] +1)
				  {
					  a[w] = a[w] +a[v];
					  p[w].push_back(v);
				  }
			  }
  
	  }//while循环
	 //S中存的是从s出发到各定点的最短路径。
	  double* sum=new double[n];
	  int v;
	  for(v = 0; v < n; v++)
	  {
		  sum[v]=0;
	  }
	  while(!S.empty())
	  {
		  int w = S.top();
		  S.pop();
		  
		  for(vector<int>::iterator ix = p[w].begin(); ix != p[w].end();++ix)
		 {
			 sum[*ix] = sum[*ix] + (g.weight[s]*g.weight[w]*a[*ix]/a[w])*(1.0+sum[w]);
		 }
		  if(w != s)
			  CB[w] = CB[w] + sum[w]/2;
	  }
	}
	ofstream osc("Betweeness.csv");
   string header="Vertex，Betweeness Centralty";
   osc<<"Vertex"<<"," <<"Betweeness Centralty"<<endl;
   for(int iv = 0 ; iv < n; iv ++)
   {
	   osc <<g.infromation[iv] <<"," << CB[iv] <<"\n";
   }
   osc <<endl;
   osc.close();
	vector< double > bec;
	vector< double > nbec;
	for(int iv = 0 ; iv < n; iv ++)
   {
	   bec.push_back( CB[iv]);
   }
	std::vector< double>::iterator max=std::max_element(std::begin(bec),std::end(bec));
	std::vector< double>::iterator min=std::min_element(std::begin(bec),std::end(bec));
	double temp,minus;
	minus=1.0/(*max-*min);
	for(int iv = 0 ; iv < n; iv ++)
   {
	   temp=(CB[iv]-*min)*minus;
	   nbec.push_back(temp);
   }
	//计算Fi
	int ir,ic;
	for (int i = 0; i < g.edge; i++)
    {   
		ir=edge[i][0];
		ic=edge[i][1];
		Cij[i]=nbec[ir]*nbec[ic]*g.cij[ir][ic];	
    }
	osc.open("core_reinforcement.csv");
   osc<<"Source"<<"," <<"Target"<<","<<"Core_reinforcement"<<endl;
   for(int i = 0 ; i < g.edge; i ++)
   {
	   ir=edge[i][0];
	   ic=edge[i][1];
	   osc <<g.infromation[ir] <<","<<g.infromation[ic] <<"," << Cij[i] <<"\n";
   }
   osc <<endl;
   osc.close();
}

//main function
int main()
{
	Graph_array g;
	string input_file;
	char str[100];
	int num=0;
	vector< string > Vertex;
	vector< string > Edge_p1;
	vector< string > Edge_p2;
	vector< double >weight;
	cout<<">>1.Please Input Vertex Dataname:< eg:Vertex.dat >"<<endl;
	cin>>str;
	input_file=str;
	ifstream is(input_file.c_str());
	is.clear();//将流刷新
	is.seekg(0,ios::beg);
	string tmp;
	getline(is,tmp);
	while(!is.eof())
	{
		string tmp,tmp1,tmp2;
		getline(is,tmp);
		istringstream stream(tmp);
		stream>>tmp1>>tmp2;
		Vertex.push_back(tmp1);
		weight.push_back(atof(tmp2.data()));
	}
	is.close();
	//kind=1带表无向图；kind=0代表有向图
	g.kind=1;
	g.vexnum=Vertex.size();
	cout<<">>2.Please Input Edge Dataname:< eg:Edge.dat >"<<endl;
	cin>>str;
	input_file=str;
	is.open(input_file.c_str());
	is.clear();//将流刷新
	is.seekg(0,ios::beg);
	getline(is,tmp);
	while(!is.eof())
	{
		string tmp,tmp1,tmp2;
		getline(is,tmp);
		istringstream stream(tmp);
		stream>>tmp1>>tmp2;
		Edge_p1.push_back(tmp1);
		Edge_p2.push_back(tmp2);
	}
	is.close();
	//计算图的边的个数
	g.edge = Edge_p1.size();
	//每个节点的信息
	g.infromation = new string[g.vexnum];
	//每个节点的权重
	g.weight = new double[g.vexnum]; 
	g.nlink = new int[g.vexnum];
	for (int i = 0; i < g.vexnum; i++ )
	{
		g.nlink[i]=0;
		g.infromation[i]=Vertex[i];
		g.weight[i] =weight[i];
	}	
	int ** edge_information;
	edge_information = new int * [g.edge];
	//搜索确定每条边两个顶点的编号
	for (int i = 0; i<g.edge; i++)
	{
		edge_information[i]= new int[2];
		for ( int iv=0;iv < g.vexnum ; iv++)
		{
			if(Vertex[iv] == Edge_p1[i])
			{
				edge_information[i][0]=iv;
				break;
			}			
		}
		for ( int iv=0;iv < g.vexnum ; iv++)
		{
			if(Vertex[iv] == Edge_p2[i])
			{
				edge_information[i][1]=iv;
				break;
			}			
		}
	}
	//构造邻接矩阵
	createGraph_by_array(edge_information, g);
	print_array(g);
	//计算Betweenness Centrality	
	FBCA(edge_information , g);
	cout<<">>3.Calculation Completed"<<endl;

	system("pause");

	return 0;	
}

