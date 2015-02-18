#ifndef GREEDY_H
#define GREEDY_H

#include "MersenneTwister.h"
#include "GreedyBase.h"
#include "Node.h"
#include <cmath>
#include <iostream>
#include <vector>
#include <queue>
#include <deque>
#include <set>
#include <stack>
#include <map>
#include <algorithm>
using namespace std;


class Greedy : public GreedyBase{
 public:
  Greedy(MTRand *RR,int nnode,int deg,Node **node);
  virtual ~Greedy();
  virtual void initiate(void);
  virtual void calibrate(void);
  virtual void prepare(bool sort);
  virtual void level(Node ***, bool sort);
  virtual void move(bool &moved);
  virtual void determMove(int *moveTo);  
  virtual void genLogTable(int maxsize);
  int Nempty;
  vector<int> mod_empty;
  
  vector<map<int,int> > mod_links;
  vector<int> mod_inlinks; // Equals degree - exit
  
  vector<int> mod_exit;
  vector<int> mod_degree;
  vector<int> mod_members;
  
 protected:
  double logChoose(int n,int k);
  int theta(int pen);
  vector<pair<int,int> >::iterator link;
  vector<int> modWnode;
};

#endif
