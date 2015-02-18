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

class treeNode{
 public:
  set<int> members;
  multimap<int,treeNode> nextLevel;
};

template <class T>
inline std::string to_string (const T& t){
  std::stringstream ss;
  ss << t;
  return ss.str();
}


void cpyNode(Node *newNode,Node *oldNode){
  
  newNode->inlinks = oldNode->inlinks;

  newNode->index = oldNode->index;

  int Nmembers = oldNode->members.size();
  newNode->members = vector<int>(Nmembers);
  for(int i=0;i<Nmembers;i++)
    newNode->members[i] = oldNode->members[i];
  
  int Nlinks = oldNode->links.size();
  newNode->links = vector<pair<int,int> >(Nlinks);
  for(int i=0;i<Nlinks;i++){
    //newNode->links[i] = oldNode->links[i];
    newNode->links[i].first = oldNode->links[i].first;
    newNode->links[i].second = oldNode->links[i].second;
  }

}



