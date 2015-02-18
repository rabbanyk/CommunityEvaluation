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
  virtual void calibrate(void){};
  virtual void prepare(bool sort){};
  virtual void level(Node ***, bool sort){};
  virtual void move(bool &moved){};
  virtual void determMove(int *moveTo){};
  virtual void genLogTable(int maxsize){};  
  int Nmod;
  int Nnode;
  int Nmem;
  int Nlinks;
  
  vector <double> logFac;  
    
  double modelLength;
  double networkLength;
  double codeLength;  
  int penalty;
  double pF;
  double score;
  
  Node **node;
  
 protected:
  
  MTRand *R;
  
};

#endif
