#define PLOGP(p) (p > 0.0 ? (p*log(p)) : 0.0)
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include "MersenneTwister.h"
#include "GreedyBase.h" 
#include "Greedy.h" 
#include "Node.h" 
#define PI 3.14159265
using namespace std;

unsigned stou(char *s);

template <class T>
inline std::string to_string (const T& t){
  std::stringstream ss;
  ss << t;
  return ss.str();
}

class treeNode{
 public:
//  bool stop;
  int level;
  double codeLength;
  set<int> members;
//  vector<int> cluster; // Two-level partition
  vector<int> rev_renumber; // Vectors to reorganize node numbers
  map<int,int> renumber;
  multimap<double,treeNode,greater<double> > nextLevel;
};

class printTreeNode{
public:
  int rank;
  double size;
  multimap<double,int,greater<double> > members;
};

class treeStats{
public:
  double twoLevelCodeLength;
  int Nmodules;
  int largeModuleLimit;
  int NlargeModules;
  double aveDepth;
  double aveSize;
};

//void delTree(treeNode &map){
//  map.level = 0;
//  map.codeLength = 0.0;
// // map.stop = true;
//  set<int>().swap(map.members);
// // vector<int>().swap(map.cluster);
//  vector<int>().swap(map.rev_renumber);
//  vector<int>().swap(map.renumber);
//  for(multimap<double,treeNode,greater<double> >::iterator it = map.nextLevel.begin(); it != map.nextLevel.end(); it++){
//    delTree(it->second);
//  }
//  multimap<double,treeNode,greater<double> >().swap(map.nextLevel);
//}

//void copyTree(treeNode &map,treeNode &cpy_map){
//  
//  cpy_map.level = map.level;
//  cpy_map.codeLength = map.codeLength;
// // cpy_map.stop = map.stop;
//  cpy_map.members = set<int>(map.members.begin(),map.members.end());
// // copy(map.cluster.begin(),map.cluster.end(),back_inserter(cpy_map.cluster));
//  copy(map.rev_renumber.begin(),map.rev_renumber.end(),back_inserter(cpy_map.rev_renumber));
//  copy(map.renumber.begin(),map.renumber.end(),back_inserter(cpy_map.renumber));
//  
//  for(multimap<double,treeNode,greater<double> >::iterator it = map.nextLevel.begin(); it != map.nextLevel.end(); it++){
//    treeNode nextLevelMap;
//    copyTree(it->second,nextLevelMap);
//    cpy_map.nextLevel.insert(make_pair(it->first,nextLevelMap));
//  }
//  
//}

void genSubNet(Node **orig_node,int Nnode, Node **sub_node,int sub_Nnode, treeNode &map,double totalDegree){
  
  vector<int>(sub_Nnode).swap(map.rev_renumber);
  //vector<int>(Nnode,-1).swap(map.renumber);
  std::map<int,int>().swap(map.renumber); 

  vector<double> degree = vector<double>(sub_Nnode,0.0);
  double exit = 0.0;
  double flow = 0.0;
  
  // Construct sub network
  set<int>::iterator it_mem = map.members.begin();
  for(int i=0;i<sub_Nnode;i++){
    int orig_nr = (*it_mem);
    map.renumber[orig_nr] = i;
    map.rev_renumber[i] = orig_nr;
    it_mem++;
  }
  it_mem = map.members.begin();
  for(int i=0;i<sub_Nnode;i++){
    int orig_nr = (*it_mem);
    int orig_Nlinks = orig_node[orig_nr]->links.size();
    sub_node[i] = new Node(i);
    for(int j=0;j<orig_Nlinks;j++){
      int orig_link = orig_node[orig_nr]->links[j].first;
      int orig_link_newnr = -1;
      std::map<int,int>::iterator it = map.renumber.find(orig_link);
      if(it != map.renumber.end())
        orig_link_newnr = it->second;
      //int orig_link_newnr = map.renumber[orig_link];
      double orig_weight = orig_node[orig_nr]->links[j].second;
      degree[i] += orig_weight;
      flow += orig_weight;
      if(orig_link_newnr < 0){
        sub_node[i]->outDegree += orig_weight;
        exit += orig_weight;
        flow += orig_weight;  
      }
      else if(orig_link < orig_nr){ // Should be <= if self-links are included
        if(map.members.find(orig_link) != map.members.end()){
          sub_node[i]->links.push_back(make_pair(orig_link_newnr,orig_weight));
          sub_node[orig_link_newnr]->links.push_back(make_pair(i,orig_weight));
        }
        
      }
    }
    it_mem++;
  }
  
  double codeLength = 0.0;
  for(int i=0;i<sub_Nnode;i++)
    codeLength -= PLOGP(degree[i]/flow);
  codeLength -= PLOGP(exit/flow);
  codeLength *= flow/totalDegree/log(2.0);
  map.codeLength = codeLength;
}

void printTree(string s,treeNode &map,vector<string> &nodeNames,vector<double> &degree,double totalDegree,ofstream *outfile,int depth,treeStats &stats){
  
//  if(map.nextLevel.size() > 0){
//    if(s.empty())
//      cout << "Code structure:" << endl << map.codeLength << " " << map.nextLevel.size() << " " << map.members.size() << endl;
//    else
//      cout << s << " " << map.codeLength << " " << map.nextLevel.size() << " " << map.members.size() << endl;
//  }
//  else{
//    cout << s << " " << map.codeLength << " " << map.nextLevel.size() << " " << map.members.size() << endl;
//  }
  if (map.nextLevel.size() > 0) {
  
    int i=1;
    for(multimap<double,treeNode,greater<double> >::iterator it = map.nextLevel.begin(); it != map.nextLevel.end(); it++){
    
      stats.Nmodules++;
      if(1.0*it->second.members.size() > 1.0*stats.largeModuleLimit)
        stats.NlargeModules++;
      
      string cpy_s(s + to_string(i) + ":");
      printTree(cpy_s,it->second,nodeNames,degree,totalDegree,outfile,depth+1,stats);
      i++;
    }
  }
  else{
    
    stats.aveDepth += 1.0*map.members.size()*depth;
		stats.aveSize += 1.0*map.members.size()*map.members.size();
    
    multimap<double,int,greater<double> > sortedMem;
    for(set<int>::iterator mem = map.members.begin(); mem != map.members.end(); mem++){
      sortedMem.insert(make_pair(degree[(*mem)],(*mem)));
    }
    int i = 1;
    for(multimap<double,int,greater<double> >::iterator mem = sortedMem.begin(); mem != sortedMem.end(); mem++){
      string cpy_s(s + to_string(i) + " " + to_string(1.0*mem->first/totalDegree) + " \"" + nodeNames[mem->second] + "\"");
      (*outfile) << cpy_s << endl;
      i++;
    } 
  }
}

void cpyNode(Node *newNode,Node *oldNode){
  
  newNode->index = oldNode->index;
  newNode->exit = oldNode->exit;
  newNode->outDegree = oldNode->outDegree;
  newNode->degree = oldNode->degree;

  int Nmembers = oldNode->members.size();
  newNode->members = vector<int>(Nmembers);
  for(int i=0;i<Nmembers;i++)
    newNode->members[i] = oldNode->members[i];
  
  int Nlinks = oldNode->links.size();
  newNode->links = vector<pair<int,double> >(Nlinks);
  for(int i=0;i<Nlinks;i++){
    newNode->links[i].first = oldNode->links[i].first;
    newNode->links[i].second = oldNode->links[i].second;
  }

}

void setCodeLength(Node **orig_node,Node **sub_node,int mod,treeNode &map, double totalDegree){
  
  double exit =  sub_node[mod]->exit;
  double flow = sub_node[mod]->degree + sub_node[mod]->exit;
  double codeLength = 0.0;
  
  for(set<int>::iterator mem = map.members.begin(); mem != map.members.end(); mem++){
    double p = orig_node[(*mem)]->degree/flow;
  	codeLength -= PLOGP(p);	
  }
  double p = exit/flow;
  codeLength -= PLOGP(p);
    
  // Weight by usage and convert to base 2
  
  codeLength *= flow/totalDegree/log(2.0);
  
  map.codeLength = codeLength;
  
}

void addNodesToMap(treeNode &map,vector<double> &size){
  
  if (map.nextLevel.size() > 0) {
    for(multimap<double,treeNode,greater<double> >::iterator it = map.nextLevel.begin(); it != map.nextLevel.end(); it++){
      addNodesToMap(it->second,size);
    }
  }
  else{
    for(set<int>::iterator mem = map.members.begin(); mem != map.members.end(); mem++){
      treeNode tmp;
      tmp.members.insert((*mem));
      map.nextLevel.insert(make_pair(size[(*mem)],tmp));
    }
  }
}

void collapseTree(multimap<double,printTreeNode,greater<double> > &collapsedmap,treeNode &map,vector<double> &size,int level){
  
  if ( ((level != 0) && (map.level <= level)) || (level == 0 && map.nextLevel.size() > 0)) {
    for(multimap<double,treeNode,greater<double> >::iterator it = map.nextLevel.begin(); it != map.nextLevel.end(); it++){
      // Check if Barabasi is in the module
      //if(map.members.find(30) != map.members.end()) 
      collapseTree(collapsedmap,it->second,size,level);
    }
  }
  else{
    printTreeNode tmp;
    double s = 0.0;
    for(set<int>::iterator mem = map.members.begin(); mem != map.members.end(); mem++){
      s += size[(*mem)];
      tmp.members.insert(make_pair(size[(*mem)],(*mem)));
    }
    tmp.size = s;
    
    collapsedmap.insert(make_pair(s,tmp));
  }
}




