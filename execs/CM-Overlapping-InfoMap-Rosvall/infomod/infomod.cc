#include "infomod.h"

using namespace std;
using std::cout;
using std::cin;
using std::endl;

unsigned stou(char *s){
  return strtoul(s,(char **)NULL,10);
}

void partition(MTRand *R,Node ***node, GreedyBase *greedy, bool silent);
void repeated_partition(MTRand *R, Node ***node, GreedyBase *greedy, bool silent,int Ntrials);
void printTree(string s,multimap<int,treeNode>::reverse_iterator it_tM,string *nodeNames,int *degree,int totalDegree,ofstream *outfile);

// Call: trade <seed> <Ntries>
int main(int argc,char *argv[]){
  
  if( argc !=4 ){
    cout << "Call: ./infomap <seed> <network.net> <# attempts>" << endl;
    exit(-1);
  }

  int Ntrials = atoi(argv[3]);  // Set number of partition attempts
  string infile = string(argv[2]);
  string networkName(infile.begin(),infile.begin() + infile.find_last_of("."));
  string line;
  string buf;

  MTRand *R = new MTRand(stou(argv[1]));

  /* Read network in Pajek format with nodes ordered 1, 2, 3, ..., N,            */
  /* each undirected link occurring only once, and integer link weights > 0.     */
  /* For more information, see http://vlado.fmf.uni-lj.si/pub/networks/pajek/.   */
  /* Example network with three nodes and                                        */
  /* three undirected and integer weighted links:                                */
  /* *Vertices 3                                                                 */
  /* 1 "Name of first node"                                                      */
  /* 2 "Name of second node"                                                     */
  /* 3 "Name of third node"                                                      */
  /* *Arcs 3                                                                     */
  /* 1 2 1                                                                       */
  /* 1 3 3                                                                       */
  /* 2 3 2                                                                       */

  cout << "Reading network " << argv[2] << "..." << flush;
  ifstream net(argv[2]);
  int Nnode = 0;
  istringstream ss;
  while(Nnode == 0){ 
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
	Nnode = atoi(buf.c_str());
      }
      else{
	cout << "the network file is not in Pajek format...exiting" << endl;
	exit(-1);
      }
    }
  }
  
  string *nodeNames = new string[Nnode];
  
 // Read node names, assuming order 1, 2, 3, ...
  for(int i=0;i<Nnode;i++){
    getline(net,line);
    int nameStart = line.find_first_of("\"");
    int nameEnd = line.find_last_of("\"");
    if(nameStart < nameEnd){
      nodeNames[i] =  string(line.begin() + nameStart + 1,line.begin() + nameEnd);
    }
    else{
      ss.clear();
      ss.str(line);
      ss >> buf; 
      ss >> nodeNames[i];
    }
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
    
  int Nlinks = 0;
  int NdoubleLinks = 0;
  map<int,set<int> > Links; 

  // Read links in format "from to", for example "1 3" (all integers) and each undirected link only ones.
  while(getline(net,line) != NULL){
    ss.clear();
    ss.str(line);
    ss >> buf;
    int linkEnd1 = atoi(buf.c_str());
    ss >> buf;
    int linkEnd2 = atoi(buf.c_str());
    
    linkEnd1--; // Nodes start at 1, but C++ arrays at 0.
    linkEnd2--;
    
    if(linkEnd2 < linkEnd1){
      int tmp = linkEnd1;
      linkEnd1 = linkEnd2;
      linkEnd2 = tmp;
    }

    // Only include multiply definied links once
    map<int,set<int> >::iterator fromLink_it = Links.find(linkEnd1);
    if(fromLink_it == Links.end()){ // new link
      set<int> toLink;
      toLink.insert(linkEnd2);
      Links.insert(make_pair(linkEnd1,toLink));
      Nlinks++;
    }
    else{
      set<int>::iterator toLink_it = fromLink_it->second.find(linkEnd2);
      if(toLink_it == fromLink_it->second.end()){ // new link
        fromLink_it->second.insert(linkEnd2);
        Nlinks++;
      }
      else{
        NdoubleLinks++;
      }
    }    
  }
  
  net.close();
  
  cout << "done! (found " << Nnode << " nodes and " << Nlinks << " links";
  if(NdoubleLinks > 0)
    cout << ", ignoring " << NdoubleLinks << " link(s) defined more than once";

  /////////// Partition network /////////////////////
  Node **node = new Node*[Nnode];
  for(int i=0;i<Nnode;i++){
    node[i] = new Node(i);
  }
  
  int NselfLinks = 0;
 
  for(map<int,set<int> >::iterator fromLink_it = Links.begin(); fromLink_it != Links.end(); fromLink_it++){
    for(set<int>::iterator toLink_it = fromLink_it->second.begin(); toLink_it != fromLink_it->second.end(); toLink_it++){
 
    	int from = fromLink_it->first;
    	int to = *toLink_it;
    	if(from == to){
      	NselfLinks++;
    	}
    	else{
      	node[from]->links.push_back(make_pair(to,1));
      	node[to]->links.push_back(make_pair(from,1));
    	}
    }
  }
  cout << ", ignoring " <<  NselfLinks << " self link(s)." << endl;
  Nlinks -= NselfLinks;

  //Swap maps to free memory
  for(map<int,set<int> >::iterator it = Links.begin(); it != Links.end(); it++)
    set<int>().swap(it->second);
  map<int,set<int> >().swap(Links);  

  // Initiation
  GreedyBase* greedy;
  greedy = new Greedy(R,Nnode,Nlinks,node);
  greedy->initiate();
  
  cout << "Now partition the network:" << endl;
  repeated_partition(R,&node,greedy,false,Ntrials);
     
  int Nmod = greedy->Nnode;
  cout << "Done! Code length " << greedy->codeLength << " in " << Nmod << " modules." << endl; 
    
  // Print partitions in Pajek's .clu format
  vector<int> clusterVec = vector<int>(Nnode);
  vector<set<int> > sortedMembers = vector<set<int> >(Nmod);
  for(int i=0;i<Nmod;i++){
    int Nmem = node[i]->members.size();
    for(int j=0;j<Nmem;j++){
      clusterVec[node[i]->members[j]] = i;
      sortedMembers[i].insert(node[i]->members[j]);
    }
  }
  ostringstream oss;
  oss << networkName << ".clu";
  ofstream outfile;
  outfile.open(oss.str().c_str());
	
  outfile << "*Vertices " << Nnode << "\x0D\x0A";
  for(int i=0;i<Nnode;i++)
    outfile << clusterVec[i]+1 << "\x0D\x0A";
  outfile.close();


  oss.str("");
  oss << networkName << ".mod";
  outfile.open(oss.str().c_str());
  outfile << "# " << Nnode << " nodes in " << Nmod << " modules.\x0D\x0A";
  for(int i=0;i<Nmod;i++)
    for(set<int>::iterator it = sortedMembers[i].begin();it != sortedMembers[i].end();it++)
      outfile << i+1 << " \"" << nodeNames[(*it)] << "\"\x0D\x0A";
  outfile.close();

 
  delete [] nodeNames;
  for(int i=0;i<greedy->Nnode;i++){
    delete node[i];
  }
  delete [] node;
  delete greedy;
  delete R;
     
}

void partition(MTRand *R,Node ***node, GreedyBase *greedy, bool silent){
  
  int Nnode = greedy->Nnode;
  int Nmem = greedy->Nmem;
  int Nlinks = greedy->Nlinks;

  Node **cpy_node = new Node*[Nnode];
  for(int i=0;i<Nnode;i++){
    cpy_node[i] = new Node();
    cpyNode(cpy_node[i],(*node)[i]);
  }
  
  int iteration = 0;
  double outer_oldCodeLength;
  do{
    outer_oldCodeLength = greedy->codeLength;
    
    if((iteration > 0) && (iteration % 2 == 0) && (greedy->Nnode > 1)){  // Partition the partition
      
      if(!silent)
	cout << "Iteration " << iteration+1 << ", moving " << flush;
      
      Node **rpt_node = new Node*[Nnode];
      for(int i=0;i<Nnode;i++){
	rpt_node[i] = new Node();
	cpyNode(rpt_node[i],cpy_node[i]);
      }
      int *subMoveTo = new int[Nnode];
      int *moveTo = new int[Nnode];
      int subModIndex = 0;
      
      for(int i=0;i<greedy->Nnode;i++){
	
	int sub_Nnode = (*node)[i]->members.size();
	
	if(sub_Nnode > 1){
	  
	  Node **sub_node = new Node*[sub_Nnode]; 
	  set<int> sub_mem;
	  for(int j=0;j<sub_Nnode;j++)
	    sub_mem.insert((*node)[i]->members[j]);
	  set<int>::iterator it_mem = sub_mem.begin();
	  int *sub_renumber = new int[Nnode];
	  int *sub_rev_renumber = new int[sub_Nnode];
	  int sub_Nlinks = 0;
	  for(int j=0;j<sub_Nnode;j++){

	    int orig_nr = (*it_mem);
	    int orig_Nlinks = cpy_node[orig_nr]->links.size();
	    sub_renumber[orig_nr] = j;
	    sub_rev_renumber[j] = orig_nr;
	    sub_node[j] = new Node(j);
	    for(int k=0;k<orig_Nlinks;k++){
	      int orig_link = cpy_node[orig_nr]->links[k].first;
	      int orig_link_newnr = sub_renumber[orig_link];
	      int orig_weight = cpy_node[orig_nr]->links[k].second;
	      if(orig_link < orig_nr){
		if(sub_mem.find(orig_link) != sub_mem.end()){
		  sub_node[j]->links.push_back(make_pair(orig_link_newnr,orig_weight));
		  sub_node[orig_link_newnr]->links.push_back(make_pair(j,orig_weight));
		  sub_Nlinks += orig_weight;
		}
	      }
	    }
	    it_mem++;
	  }
	  

	  GreedyBase* sub_greedy;
	  sub_greedy = new Greedy(R,sub_Nnode,sub_Nlinks,sub_node);
	  sub_greedy->initiate();
	  partition(R,&sub_node,sub_greedy,true);

	  for(int j=0;j<sub_greedy->Nnode;j++){
	    int Nmembers = sub_node[j]->members.size();
	    for(int k=0;k<Nmembers;k++){
	      subMoveTo[sub_rev_renumber[sub_node[j]->members[k]]] = subModIndex;
	    }
	    moveTo[subModIndex] = i;
	    subModIndex++;
	    delete sub_node[j];
	  }
	  
	  delete [] sub_node;
	  delete sub_greedy;
	  delete [] sub_renumber;
	  delete [] sub_rev_renumber;
	  
	}
	else{
	  
	  subMoveTo[(*node)[i]->members[0]] = subModIndex;
	  moveTo[subModIndex] = i;
	  
	  subModIndex++;
	  
	}
      }
      
      for(int i=0;i<greedy->Nnode;i++)
	delete (*node)[i];
      delete [] (*node);
      
      greedy->Nnode = Nnode;
      greedy->Nmod = Nnode;
      greedy->Nmem = Nmem;
      greedy->Nlinks = Nlinks;
      greedy->node = rpt_node;
      greedy->initiate();
      greedy->determMove(subMoveTo);
      greedy->level(node,false); 
      greedy->determMove(moveTo);
      delete [] subMoveTo;
      (*node) = rpt_node;
      delete [] moveTo;
      
      outer_oldCodeLength = greedy->codeLength;
      
      if(!silent)
	cout << greedy->Nnode << " modules, looping ";
      
    }
    else if(iteration > 0){
 
      if(!silent)
	cout << "Iteration " << iteration+1 << ", moving " << Nnode << " nodes, looping ";

      Node **rpt_node = new Node*[Nnode];
      for(int i=0;i<Nnode;i++){
	rpt_node[i] = new Node();
	cpyNode(rpt_node[i],cpy_node[i]);
      }
      
     

      int *moveTo = new int[Nnode];
      for(int i=0;i<greedy->Nnode;i++){
	int Nmembers = (*node)[i]->members.size();
	for(int j=0;j<Nmembers;j++){
	  moveTo[(*node)[i]->members[j]] = i;
	}
      }

           
      
      for(int i=0;i<greedy->Nnode;i++)
	delete (*node)[i];
      delete [] (*node);

      greedy->Nnode = Nnode;
      greedy->Nmod = Nnode;
      greedy->Nmem = Nmem;
      greedy->Nlinks = Nlinks;
      greedy->node = rpt_node;
      greedy->initiate();
      greedy->determMove(moveTo);
      delete [] moveTo;
      (*node) = rpt_node;
    
   

    }
    else{
      
      if(!silent)
	cout << "Iteration " << iteration+1 << ", moving " << Nnode << " nodes, looping ";
      
    }

    double oldCodeLength;
    do{

      oldCodeLength = greedy->codeLength;
      bool moved = true;
      int Nloops = 0;
      while(moved){

	moved = false;
	double inner_oldCodeLength = greedy->codeLength;
	greedy->move(moved);
	Nloops++;
	if(inner_oldCodeLength-greedy->codeLength < 1.0e-10)
	  moved = false;
      }
      
      greedy->level(node,true);

      if(!silent)
	cout << Nloops << " ";

    } while(oldCodeLength - greedy->codeLength >  1.0e-10);
 
    iteration++;
    if(!silent)
      cout << "times between mergings to code length " <<  greedy->codeLength << " in " << greedy->Nmod << " modules." << " Link penalty is " << greedy->penalty << " with penalty factor " << greedy->pF << "." <<  endl;
  
  } while(outer_oldCodeLength - greedy->codeLength > 1.0e-10);

  
  for(int i=0;i<Nnode;i++)
    delete cpy_node[i];
  delete [] cpy_node;
 
}

void repeated_partition(MTRand *R, Node ***node, GreedyBase *greedy, bool silent,int Ntrials){
  
  double shortestCodeLength = 1.0e10;
  int Nnode = greedy->Nnode;
  int *cluster = new int[Nnode];

  for(int trial = 0; trial<Ntrials;trial++){
    
    if(!silent && greedy->pF < 0.5)
      cout << "Attempt " << trial+1 << "/" << Ntrials << endl;
    
    Node **cpy_node = new Node*[Nnode];
    for(int i=0;i<Nnode;i++){
      cpy_node[i] = new Node();
      cpyNode(cpy_node[i],(*node)[i]);
    }
    
    greedy->Nnode = Nnode;
    greedy->Nmod = Nnode;
    greedy->node = cpy_node;
    greedy->initiate();
    
    partition(R,&cpy_node,greedy,silent);
 
    if(!silent && greedy->penalty == 0){
      if(greedy->codeLength < shortestCodeLength){
	
	shortestCodeLength = greedy->codeLength;
	
	// Store best partition
	for(int i=0;i<greedy->Nnode;i++){
	  for(vector<int>::iterator mem = cpy_node[i]->members.begin(); mem != cpy_node[i]->members.end(); mem++){
	    cluster[(*mem)] = i;
	  }
	}
      }
      greedy->pF = 0.0;
    }
    else{
            
      if(greedy->pF < 5){
	trial--;
	greedy->pF += 0.5;
	cout << "One or more modules have more links across boundary than within, increasing the penalty factor and trying again..." << endl;
      }
      else{
	cout << "One or more modules still have more links across boundary than within, giving up..." << endl;
	// Store partition
	for(int i=0;i<greedy->Nnode;i++){
	  for(vector<int>::iterator mem = cpy_node[i]->members.begin(); mem != cpy_node[i]->members.end(); mem++){
	    cluster[(*mem)] = i;
	  }
	}
      }
    }
    
    for(int i=0;i<greedy->Nnode;i++){
      delete cpy_node[i];
    }
    delete [] cpy_node;
    
  }
  
  // Commit best partition
  greedy->Nnode = Nnode;
  greedy->Nmod = Nnode;
  greedy->node = (*node);
  greedy->initiate();
  greedy->determMove(cluster);
  greedy->level(node,true);
  delete [] cluster;
  
}

void printTree(string s,multimap<int,treeNode>::reverse_iterator it_tM,string *nodeNames,int *degree,int totalDegree,ofstream *outfile){
 
  multimap<int,treeNode>::reverse_iterator it;
  if(it_tM->second.nextLevel.size() > 0){
    int i=1;
    for(it = it_tM->second.nextLevel.rbegin(); it != it_tM->second.nextLevel.rend(); it++){
      string cpy_s(s + to_string(i) + ":");
      printTree(cpy_s,it,nodeNames,degree,totalDegree,outfile);
      i++;
    }
  }
  else{
    multimap<int,int> sortedMem;
    for(set<int>::iterator mem = it_tM->second.members.begin(); mem != it_tM->second.members.end(); mem++){
      sortedMem.insert(make_pair(degree[(*mem)],(*mem)));
    }
    int i = 1;
    for(multimap<int,int>::reverse_iterator mem = sortedMem.rbegin(); mem != sortedMem.rend(); mem++){
      string cpy_s(s + to_string(i) + " " + to_string(1.0*mem->first/totalDegree) + " \"" + nodeNames[mem->second] + "\"");
      (*outfile) << cpy_s << endl;
      i++;
    } 
  }  
}
