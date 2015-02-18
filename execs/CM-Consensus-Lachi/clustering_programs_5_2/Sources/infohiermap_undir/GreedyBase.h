#ifndef GREEDYBASE_H
#define GREEDYBASE_H
#include "MersenneTwister.h"
#include <cstdio>
#include <vector>
using namespace std;
// forward declaration
class Node;
class GreedyBase{
 public:
  GreedyBase(){};
  virtual ~GreedyBase(){};
  virtual void initiate(void){};
  virtual void tune(void){};
  virtual void calibrate(void){};
  virtual void prepare(bool sort){};
  virtual void level(Node ***, bool sort){};
  virtual void move(bool &moved){};
  virtual void determMove(vector<int> &moveTo){};
  int Nmod;
  int Nnode;
 
  double degree;
  double log2;

  double outDegree;
  double out_log_out;
  double exit;
  double exitDegree;
  double exit_log_exit;
  double degree_log_degree;
  double nodeDegree_log_nodeDegree;
    
  double indexLength;
  double moduleLength;
  double codeLength;
  double twoLevelCodeLength;
 
  Node **node;
 
 protected:

  MTRand *R;
 
};

#endif
