#include "Node.h"

Node::~Node(){

  members.clear();
  links.clear();

}

Node::Node(){
}


Node::Node(int nodenr){
  
  index = nodenr;
  inlinks = 0;
  //   exit = 0;
  //   degree = 0;
  members.push_back(nodenr);
  
}
