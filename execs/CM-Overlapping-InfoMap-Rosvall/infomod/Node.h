#ifndef NODE_H
#define NODE_H

#include <cmath>
#include <iostream>
#include <vector>
#include <set>
#include <list>
#include <map>
#include <stack>
#include <algorithm>
#include <iomanip>
#include <fstream>
#include <cstring>
class Node;
using namespace std;

class Node{
 public:
  ~Node();
  Node();
  
  Node(int modulenr);

  vector<int> members; // If module, lists member nodes in module
  vector<pair<int,int> > links; // List of identities and link weight of connected nodes/modules  

  int inlinks; // Number of links to its own module

  /* int exit; // total weight of links to other nodes / modules */
  /*   int degree; // total degree of node / module */
  int index; // the node / module identity
  
 protected:
  
};

#endif
