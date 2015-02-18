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

class Network{
  
public:
  
  Network(string netname);
  string name;
  int Nnode;
  int Nlinks;
  double totNodeWeights;
  vector<string> nodeNames;
  vector<double> nodeWeights;
  map<int,map<int,double> > Links;
  
};

Network::Network(string netname){
  name = netname;
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
//class treeNode{
// public:
//  set<int> members;
//  multimap<double,treeNode> nextLevel;
//};

template <class T>
inline std::string to_string (const T& t){
  std::stringstream ss;
  ss << t;
  return ss.str();
}

void genSubNet(Node **orig_node,int Nnode, Node **sub_node,int sub_Nnode, treeNode &map){
    
  vector<int>(sub_Nnode).swap(map.rev_renumber);
  //  vector<int>(Nnode,-1).swap(map.renumber);
  std::map<int,int>().swap(map.renumber); 
  
  double outFlow = 0.0;
  double size = 0.0;

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
    int orig_NoutLinks = orig_node[orig_nr]->outLinks.size();
    int orig_NinLinks = orig_node[orig_nr]->inLinks.size();
    sub_node[i] = new Node(i,orig_node[orig_nr]->teleportWeight);
    sub_node[i]->size = orig_node[orig_nr]->size;
    sub_node[i]->selfLink = orig_node[orig_nr]->selfLink; // Take care of self-link
    size += sub_node[i]->size;
    
    for(int j=0;j<orig_NoutLinks;j++){
      int orig_link = orig_node[orig_nr]->outLinks[j].first;
      //int orig_link_newnr = map.renumber[orig_link];
      int orig_link_newnr = -1;
      std::map<int,int>::iterator it = map.renumber.find(orig_link);
      if(it != map.renumber.end())
        orig_link_newnr = it->second;
      double orig_weight = orig_node[orig_nr]->outLinks[j].second;
      if(orig_link_newnr < 0){
        sub_node[i]->outFlow += orig_weight;
        outFlow += orig_weight;
      }
      else if(orig_link < orig_nr){
        if(map.members.find(orig_link) != map.members.end()){
          sub_node[i]->outLinks.push_back(make_pair(orig_link_newnr,orig_weight));
          sub_node[orig_link_newnr]->inLinks.push_back(make_pair(i,orig_weight));
        }
      }
    }
    
    for(int j=0;j<orig_NinLinks;j++){
      int orig_link = orig_node[orig_nr]->inLinks[j].first;
      //int orig_link_newnr = map.renumber[orig_link];
      int orig_link_newnr = -1;
      std::map<int,int>::iterator it = map.renumber.find(orig_link);
      if(it != map.renumber.end())
        orig_link_newnr = it->second;
      double orig_weight = orig_node[orig_nr]->inLinks[j].second;
      if(orig_link < orig_nr){
        if(map.members.find(orig_link) != map.members.end()){
          sub_node[i]->inLinks.push_back(make_pair(orig_link_newnr,orig_weight));
          sub_node[orig_link_newnr]->outLinks.push_back(make_pair(i,orig_weight));
        }
      }
    }
    
    it_mem++;
  }
      
  double totFlow = size + outFlow;
  
  double codeLength = 0.0;
  for(int i=0;i<sub_Nnode;i++)
    codeLength -= PLOGP(sub_node[i]->size/totFlow);
  codeLength -= PLOGP(outFlow/totFlow);
  codeLength *= totFlow;
    
//  if(sub_Nnode == 1)
//    cout << "Level with one node: " << sub_node[0]->size << " " << outFlow << " " << codeLength/log(2.0) << endl;
  
  map.codeLength = codeLength;
}

void setCodeLength(Node **orig_node,Node **sub_node,int mod,treeNode &map){
  
  //double enter =  sub_node[mod]->enter;
  double exit =  sub_node[mod]->exit;
  double flow = sub_node[mod]->size + sub_node[mod]->exit;
  double codeLength = 0.0;
    
  for(set<int>::iterator mem = map.members.begin(); mem != map.members.end(); mem++){
    double p = orig_node[(*mem)]->size/flow;
    codeLength -= PLOGP(p);
  }
  double p = exit/flow;
  codeLength -= PLOGP(p);
    
  // Weight by usage
  codeLength *= flow;
  
  map.codeLength = codeLength;
  
}

void cpyNode(Node *newNode,Node *oldNode){
  
  newNode->index = oldNode->index;
  newNode->enter = oldNode->enter;
  newNode->exit = oldNode->exit;
  newNode->size = oldNode->size;
  newNode->outFlow = oldNode->outFlow;
  newNode->teleportWeight = oldNode->teleportWeight;
  
  int Nmembers = oldNode->members.size();
  newNode->members = vector<int>(Nmembers);
  for(int i=0;i<Nmembers;i++)
    newNode->members[i] = oldNode->members[i];
  
  newNode->selfLink = oldNode->selfLink;
  
  int NoutLinks = oldNode->outLinks.size();
  newNode->outLinks = vector<pair<int,double> >(NoutLinks);
  for(int i=0;i<NoutLinks;i++){
    newNode->outLinks[i].first = oldNode->outLinks[i].first;
    newNode->outLinks[i].second = oldNode->outLinks[i].second;
  }
  
  int NinLinks = oldNode->inLinks.size();
  newNode->inLinks = vector<pair<int,double> >(NinLinks);
  for(int i=0;i<NinLinks;i++){
    newNode->inLinks[i].first = oldNode->inLinks[i].first;
    newNode->inLinks[i].second = oldNode->inLinks[i].second;
  }
  
}

void printTree(string s,treeNode &map,vector<string> &nodeNames,vector<double> &size,ofstream *outfile,int depth,treeStats &stats){
  
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
      printTree(cpy_s,it->second,nodeNames,size,outfile,depth+1,stats);
      i++;
    }
  }
  else{
    
    stats.aveDepth += 1.0*map.members.size()*depth;
		stats.aveSize += 1.0*map.members.size()*map.members.size();

    multimap<double,int,greater<double> > sortedMem;
    for(set<int>::iterator mem = map.members.begin(); mem != map.members.end(); mem++){
      sortedMem.insert(make_pair(size[(*mem)],(*mem)));
    }
    int i = 1;
    for(multimap<double,int,greater<double> >::iterator mem = sortedMem.begin(); mem != sortedMem.end(); mem++){
      string cpy_s(s + to_string(i) + " " + to_string(1.0*mem->first) + " \"" + nodeNames[mem->second] + "\"");
      (*outfile) << cpy_s << endl;
      i++;
    } 
  }
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
  
  //if (map.nextLevel.size() > 0) {
  if ( ((level != 0) && (map.level <= level)) || (level == 0 && map.nextLevel.size() > 0)) {
    for(multimap<double,treeNode,greater<double> >::iterator it = map.nextLevel.begin(); it != map.nextLevel.end(); it++){
      // Check if xxxx is in the module
      //  if(map.members.find(5299) != map.members.end()) 
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

void loadPajekNet(Network &network){
  
  string line;
  string buf;
  
  /* Read network in Pajek format with nodes ordered 1, 2, 3, ..., N,            */
  /* each directed link occurring only once, and link weights > 0.               */
  /* For more information, see http://vlado.fmf.uni-lj.si/pub/networks/pajek/.   */
  /* Node weights are optional and sets the relative proportion to which         */
  /* each node receives teleporting random walkers. Default value is 1.          */
  /* Example network with three nodes and four directed and weighted links:      */
  /* *Vertices 3                                                                 */
  /* 1 "Name of first node" 1.0                                                  */
  /* 2 "Name of second node" 2.0                                                 */
  /* 3 "Name of third node" 1.0                                                  */
  /* *Arcs 4                                                                     */
  /* 1 2 1.0                                                                     */
  /* 1 3 1.7                                                                     */
  /* 2 3 2.0                                                                     */
  /* 3 2 1.2                                                                     */
  
  cout << "Reading network " << network.name << "..." << flush;
  ifstream net(network.name.c_str());
  network.Nnode = 0;
  istringstream ss;
  while(network.Nnode == 0){ 
    if(getline(net,line) == NULL){
      cout << "the network file is not in Pajek format...exiting" << endl;
      exit(-1);
    }
    else{
      ss.clear();
      ss.str(line);
      ss >> buf;
      if(buf == "*Vertices" || buf == "*vertices" || buf == "*VERTICES"){
        ss >> buf;
        network.Nnode = atoi(buf.c_str());
      }
      else{
        cout << "the network file is not in Pajek format...exiting" << endl;
        exit(-1);
      }
    }
  }
  
  network.nodeNames = vector<string>(network.Nnode);
  network.nodeWeights = vector<double>(network.Nnode,1.0);
  network.totNodeWeights = 0.0;
  
  // Read node names, assuming order 1, 2, 3, ...
  for(int i=0;i<network.Nnode;i++){
    getline(net,line);
    int nameStart = line.find_first_of("\"");
    int nameEnd = line.find_last_of("\"");
    if(nameStart < nameEnd){
      network.nodeNames[i] =  string(line.begin() + nameStart + 1,line.begin() + nameEnd);
      line = string(line.begin() + nameEnd + 1, line.end());
      ss.clear();
      ss.str(line);
    }
    else{
      ss.clear();
      ss.str(line);
      ss >> buf; 
      ss >> network.nodeNames[i];
    }
    
    buf = "1";
    ss >> buf;
    double nodeWeight = atof(buf.c_str());
    if(nodeWeight <= 0.0)
      nodeWeight = 1.0;
    network.nodeWeights[i] = nodeWeight;
    network.totNodeWeights += nodeWeight; 
  }
  
  // Read the number of links in the network
  getline(net,line);
  ss.clear();
  ss.str(line);
  ss >> buf;
  
  if(buf != "*Edges" && buf != "*edges" && buf != "*Arcs" && buf != "*arcs"){
    cout << endl << "Number of nodes not matching, exiting" << endl;
    exit(-1);
  }

  network.Nlinks = 0;
  int NdoubleLinks = 0;
	//map<int,map<int,double> > Links;
  
  // Read links in format "from to weight", for example "1 3 0.7"
  while(getline(net,line) != NULL){
    ss.clear();
    ss.str(line);
    ss >> buf;
    int linkEnd1 = atoi(buf.c_str());
    ss >> buf;
    int linkEnd2 = atoi(buf.c_str());
    buf.clear();
    ss >> buf;
    double linkWeight;
    if(buf.empty()) // If no information 
      linkWeight = 1.0;
    else
      linkWeight = atof(buf.c_str());
    linkEnd1--; // Nodes start at 1, but C++ arrays at 0.
    linkEnd2--;
     
    // Aggregate link weights if they are definied more than once
    map<int,map<int,double> >::iterator fromLink_it = network.Links.find(linkEnd1);
    if(fromLink_it == network.Links.end()){ // new link
      map<int,double> toLink;
      toLink.insert(make_pair(linkEnd2,linkWeight));
      network.Links.insert(make_pair(linkEnd1,toLink));
      network.Nlinks++;
    }
    else{
      map<int,double>::iterator toLink_it = fromLink_it->second.find(linkEnd2);
      if(toLink_it == fromLink_it->second.end()){ // new link
        fromLink_it->second.insert(make_pair(linkEnd2,linkWeight));
        network.Nlinks++;
      }
      else{
        toLink_it->second += linkWeight;
        NdoubleLinks++;
      }
    }
    
  }

  net.close();
  
  cout << "done! (found " << network.Nnode << " nodes and " << network.Nlinks << " links";
  if(NdoubleLinks > 0)
    cout << ", aggregated " << NdoubleLinks << " link(s) defined more than once";
  
}

void loadLinkList(Network &network){
  
  string line;
  string buf;
  istringstream ss;
  
  /* Read network in the format "FromNodeId    ToNodeId    [LinkWeight] ",    */
  /* not assuming a complete list of nodes 1..maxnode                         */
  
  cout << "Reading network " << network.name << "..." << flush;
  ifstream net(network.name.c_str());
  
  network.Nlinks = 0;
  int NdoubleLinks = 0;
	//map<int,map<int,double> > Links;
    
  // Read links in format "from to weight", for example "1 3 0.7"
  while(getline(net,line) != NULL){
    
    if(line[0] != '#'){
      
      ss.clear();
      ss.str(line);
      ss >> buf;
      
      int linkEnd1 = atoi(buf.c_str());
      ss >> buf;
      int linkEnd2 = atoi(buf.c_str());
      buf.clear();
      ss >> buf;
      double linkWeight;
      if(buf.empty()) // If no information 
        linkWeight = 1.0;
      else
        linkWeight = atof(buf.c_str());
      //      linkEnd1--; // Nodes start at 1, but C++ arrays at 0.
      //      linkEnd2--;
    
      //    set<int> blacklist;
      //    blacklist.insert(4789);
      //    blacklist.insert(5813);
      //    blacklist.insert(5125);
      
      //    if(blacklist.find(linkEnd1) == blacklist.end() && blacklist.find(linkEnd2) == blacklist.end()){ 
      
      // Aggregate link weights if they are definied more than once
      map<int,map<int,double> >::iterator fromLink_it = network.Links.find(linkEnd1);
      if(fromLink_it == network.Links.end()){ // new link
        map<int,double> toLink;
        toLink.insert(make_pair(linkEnd2,linkWeight));
        network.Links.insert(make_pair(linkEnd1,toLink));
        network.Nlinks++;
      }
      else{
        map<int,double>::iterator toLink_it = fromLink_it->second.find(linkEnd2);
        if(toLink_it == fromLink_it->second.end()){ // new link
          fromLink_it->second.insert(make_pair(linkEnd2,linkWeight));
          network.Nlinks++;
        }
        else{
          toLink_it->second += linkWeight;
          NdoubleLinks++;
        }
      }
      
    }
    //  }
    
  }
    
  
  net.close();
    
  network.Nnode = network.Links.size();
  network.nodeNames = vector<string>(network.Nnode);
  network.nodeWeights = vector<double>(network.Nnode,1.0);
  network.totNodeWeights = 1.0*network.Nnode;
  
  int nodeCounter = 0;
  map<int,int> renumber;
  bool renum = false;
  for(map<int,map<int,double> >::iterator it = network.Links.begin(); it != network.Links.end(); it++){
    network.nodeNames[nodeCounter] = to_string(it->first); 
    renumber.insert(make_pair(it->first,nodeCounter));
    if(nodeCounter != it->first)
      renum = true;
    nodeCounter++;
  }
  
  if(renum){
    
    map<int,map<int,double> > newLinks;
    map<int,map<int,double> >::iterator newLinks_it = newLinks.begin();
    for(map<int,map<int,double> >::iterator it = network.Links.begin(); it != network.Links.end(); it++){
      int fromLink = renumber.find(it->first)->second;
      map<int,double> toLink;
      map<int,double>::iterator toLink_it = toLink.begin();
      for(map<int,double>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++){
        toLink_it = toLink.insert(toLink_it,make_pair(renumber.find(it2->first)->second,it2->second));
      }
      newLinks_it = newLinks.insert(newLinks_it,make_pair(fromLink,toLink));  
    }
   
    network.Links.swap(newLinks);
    
  }
  
  cout << "done! (found " << network.Nnode << " nodes and " << network.Nlinks << " links";
  if(NdoubleLinks > 0)
    cout << ", aggregated " << NdoubleLinks << " link(s) defined more than once";
  
}




