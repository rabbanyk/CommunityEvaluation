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
  virtual void level(Node ***,bool sort){};
  virtual void move(bool &moved){};
  virtual void determMove(vector<int> &moveTo){};
  virtual void eigenvector(void){};
  virtual void eigenfactor(void){};
  virtual void collapseNodes(void){};
  int Nmod;
  int Nnode;
  
  
  double outFlow;
  double out_log_out;
  double enter;
  double enterFlow;
  double enter_log_enter;
  double exit_log_exit;
  double size_log_size;
  double nodeSize_log_nodeSize;
  
  double indexLength;
  double moduleLength;
  double codeLength;
 
  Node **node;
  double alpha,beta;

  bool initRun;
  
 protected:

  MTRand *R;
  
};

#endif
