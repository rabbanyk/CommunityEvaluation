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
  map<pair<int,int>,double> Links;
  
};

Network::Network(string netname){
  name = netname;
}

class treeNode{
 public:
  double exit;
  multimap<double,pair<int,string>,greater<double> > members;
  multimap<double,treeNode,greater<double> > nextLevel;
};

template <class T>
inline std::string to_string (const T& t){
  std::stringstream ss;
  ss << t;
  return ss.str();
}

void cpyNode(Node *newNode,Node *oldNode){
  
  newNode->index = oldNode->index;
  newNode->exit = oldNode->exit;
  newNode->size = oldNode->size;
  newNode->teleportWeight = oldNode->teleportWeight;
  newNode->danglingSize = oldNode->danglingSize;
  
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

  double newLinkWeight;
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
     
    newLinkWeight = network.Links[make_pair(linkEnd1,linkEnd2)] += linkWeight;
    if(newLinkWeight > linkWeight)
      NdoubleLinks++;

  //   // Aggregate link weights if they are definied more than once
  //   map<int,map<int,double> >::iterator fromLink_it = network.Links.find(linkEnd1);
  //   if(fromLink_it == network.Links.end()){ // new link
  //     map<int,double> toLink;
  //     toLink.insert(make_pair(linkEnd2,linkWeight));
  //     network.Links.insert(make_pair(linkEnd1,toLink));
  //     network.Nlinks++;
  //   }
  //   else{
  //     map<int,double>::iterator toLink_it = fromLink_it->second.find(linkEnd2);
  //     if(toLink_it == fromLink_it->second.end()){ // new link
  //       fromLink_it->second.insert(make_pair(linkEnd2,linkWeight));
  //       network.Nlinks++;
  //     }
  //     else{
  //       toLink_it->second += linkWeight;
  //       NdoubleLinks++;
  //     }
  //   }
  //   
  }

  net.close();

  network.Nlinks = network.Links.size();
  
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
  
  double newLinkWeight;
	set<int> Nodes;
  
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
      
			// To keep track of all nodes for renaming 0...N-1
			Nodes.insert(linkEnd1);
			Nodes.insert(linkEnd2);
      
      newLinkWeight = network.Links[make_pair(linkEnd1,linkEnd2)] += linkWeight;
			if(newLinkWeight > linkWeight)
        NdoubleLinks++;
    
    }
    
  }
  
  net.close();
  
  network.Nnode = Nodes.size();
  network.Nlinks = network.Links.size();
  network.nodeNames = vector<string>(network.Nnode);
  
  int nodeCounter = 0;
  map<int,int> renumber;
  bool renum = false;
  for(set<int>::iterator it = Nodes.begin(); it != Nodes.end(); it++){
    network.nodeNames[nodeCounter] = to_string((*it)); 
    renumber.insert(make_pair((*it),nodeCounter));
    if(nodeCounter != (*it))
      renum = true;
    nodeCounter++;	
  }	
  
  if(renum){
        
    map<pair<int,int>,double> newLinks;
    for(map<pair<int,int>,double>::iterator it = network.Links.begin(); it != network.Links.end(); it++)
      newLinks.insert(make_pair(make_pair(renumber.find(it->first.first)->second,renumber.find(it->first.second)->second),it->second));
    
    network.Links.swap(newLinks);
    
  }
  
  cout << "done! (found " << network.Nnode << " nodes and " << network.Nlinks << " links";
  if(NdoubleLinks > 0)
    cout << ", aggregated " << NdoubleLinks << " link(s) defined more than once";
    
}





  




