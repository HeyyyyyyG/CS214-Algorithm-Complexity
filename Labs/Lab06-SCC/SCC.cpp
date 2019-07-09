#include <vector>
#include <iostream>
#include <fstream>
using namespace std;
//Please put this source code in the same directory with scc.in
//And do NOT change the file name.
/*
This function computes the number of Strongly Connected Components in a graph
Args:
    n: The number of nodes in the graph. The nodes are indexed as 0~n-1
    edge: The edges in the graph. For each element (a,b) in edge, it means
          there is a directed edge from a to b
          Notice that there may exists multiple edge and self-loop
Return:
    The number of strongly connected components in the graph.
*/
void previsit(int v,int &clock,int *pre)
{
    pre[v]=clock;
    clock++;
}
void postvisit(int v,int &clock,int *post)
{
    post[v]=clock;
    clock++;
}
void explore(int v,vector<pair<int,int> >& edge,int *&visited,int nodeNum,int &clock,int*& pre,int*& post)
{
    visited[v]=true;//cout<<v;

    previsit(v,clock,pre);
    for (int i=0;i<edge.size();i++)
    {
        if(edge.at(i).first==v)
        {
            int u=edge.at(i).second;
            if(visited[u]==0)
            explore(u,edge,visited,nodeNum,clock,pre,post);
        }

    }
    postvisit(v,clock,post);
}

void dfs(int nodeNum,vector<pair<int,int> >& edge,int *&pre,int *&post)
{
    int *visited=new int[nodeNum];
    int clock=1;

    for(int i=0;i<nodeNum;i++)  //initialize
    {
        visited[i]=0;
        pre[i]=0;
        post[i]=0;
    }

    for (int i=0;i<nodeNum;i++)
        if (visited[i]==0)
            explore(i,edge,visited,nodeNum,clock,pre,post);

    //for(int i=0;i<nodeNum;i++)
      //  cout<<pre[i]<<","<<post[i]<<endl;
}
void explore2(int v,vector<pair<int,int> >& edge,int *&visited,int nodeNum)
{
    visited[v]=true;//cout<<v;


    for (int i=0;i<edge.size();i++)
    {
        if(edge.at(i).first==v)
        {
            int u=edge.at(i).second;
            if(visited[u]==0)
            explore2(u,edge,visited,nodeNum);
        }

    }

}
void deleteNode(int n,vector<pair<int,int> >& edge,int *&pre,int *&post,int &nodeNumInStrong)
{

    int postmax=0;int postmaxnode=0;
    for(int i=0;i<n&&post[i]>0;i++) //find the max  post node
    {
        if(post[i]>postmax)
        {
            postmax=post[i];postmaxnode=i;
        }

    }

    //dfs from postmaxnode
    vector<int> nodeInStrong;
    int *visited=new int[n];

    for(int i=0;i<n;i++)  //initialize
    {
        visited[i]=0;
    }


        explore2(postmaxnode,edge,visited,n);



    for(int j=0;j<n;j++)
    {
        if(visited[j])
        {
            nodeNumInStrong++;
            //cout<<j<<" ";
            for (int i=0;i<edge.size();i++)  //delete all edges related to the strong connected components,i.e. delete postmax1node
            {
                if(edge.at(i).first==j||edge.at(i).second==j)
                    edge.erase(edge.begin()+i);
            }
        }
    }



}
int SCC2(int n, vector<pair<int,int> >& edge,int &strong) {     //n:nodeNum
    vector<pair<int,int> > edgeR;
    int *pre=new int[n];
    int *post=new int[n];

    for(int i=0;i<edge.size();i++) //reverse the previous graph
        edgeR.push_back(pair<int,int>(edge.at(i).second,edge.at(i).first));

    dfs(n,edgeR,pre,post);

    int nodeNumInStrong=0;
    deleteNode(n,edge,pre,post,nodeNumInStrong); //delete the nodes in the first strongly connected components
    strong++;
//cout<<endl;
    delete pre;
    delete post;

    if(n-nodeNumInStrong==0)
        return strong;
    else return SCC2(n-nodeNumInStrong,edge,strong);
}
int SCC(int n,vector<pair<int,int> >& edge)
{
    int strong=0; //the number of strongly connected components
    return SCC2(n,edge,strong);
}
//Please do NOT modify anything in main(). Thanks!
int main()
{
    int m,n;
    vector<pair<int,int> > edge;
    ifstream fin;
    ofstream fout;
    fin.open("scc.in.txt");
    fin>>n>>m;
    int tmp1,tmp2;
    for(int i=0;i<m;i++)
    {
        fin>>tmp1>>tmp2;
        edge.push_back(pair<int,int>(tmp1,tmp2));

    }
    fin.close();
    int ans=SCC(n,edge);
    fout.open("scc.out");
    fout<<ans<<'\n';
    fout.close();

    return 0;
}
