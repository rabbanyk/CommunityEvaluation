#include "infomap.h"

using namespace std;
using std::cout;
using std::cin;
using std::endl;

unsigned stou(char *s){
  return strtoul(s,(char **)NULL,10);
}
double fast_hierarchical_partition(MTRand *R, Node **node, treeNode &map, double totalDegree, int Nnode,double &twoLevelCodeLength, bool deep);
double hierarchical_partition(MTRand *R, Node **node, treeNode &map, double totalDegree, int Nnode, double recursive);
double repeated_hierarchical_partition(string networkName,vector<double> &degree, vector<string> &nodeNames, MTRand *R, Node **orig_node, treeNode &map, double totalDegree, int Nnode,int Ntrials, double recursive, treeStats &stats);
void partition(MTRand *R,Node ***node, GreedyBase *greedy, bool silent);
void repeated_partition(MTRand *R, Node ***node, GreedyBase *greedy, bool silent,int Ntrials);

// Call: trade <seed> <Ntries>
int main(int argc,char *argv[]){
  
  if(argc < 4 || argc > 5){
    cout << "Call: ./infomap <seed> <network.net> <# attempts> <recursive[0-1]>" << endl;
    exit(-1);
  }
  
  int Ntrials = atoi(argv[3]);  // Set number of partition attempts
  double recursive = 0.0;
  if(argc == 5)
    recursive = atoi(argv[4]);
  string infile = string(argv[2]);
  string networkName(infile.begin(),infile.begin() + infile.find("."));
  string line;
  string buf;
  
  MTRand *R = new MTRand(stou(argv[1]));
  
  /* Read network in Pajek format with nodes ordered 1, 2, 3, ..., N,            */
  /* each undirected link occurring only once, and link weights > 0.             */
  /* (if a link is defined more than once, weights are aggregated)               */   
  /* For more information, see http://vlado.fmf.uni-lj.si/pub/networks/pajek/.   */
  /* Example network with three nodes and                                        */
  /* three undirected weighted links:                                            */
  /* *Vertices 3                                                                 */
  /* 1 "Name of first node"                                                      */
  /* 2 "Name of second node"                                                     */
  /* 3 "Name of third node"                                                      */
  /* *Edges 3                                                                    */
  /* 1 2 1.0                                                                     */
  /* 1 3 3.3                                                                     */
  /* 2 3 2.2                                                                     */
  
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
  
  vector<string> nodeNames(Nnode);
  
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
  map<int,map<int,double> > Links; 
  // Read links in format "from to weight", for example "1 3 2" (all integers) and each undirected link only ones (weight is optional).
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
    
    if(linkEnd2 < linkEnd1){
      int tmp = linkEnd1;
      linkEnd1 = linkEnd2;
      linkEnd2 = tmp;
    }
    
    // Aggregate link weights if they are definied more than once
    map<int,map<int,double> >::iterator fromLink_it = Links.find(linkEnd1);
    if(fromLink_it == Links.end()){ // new link
      map<int,double> toLink;
      toLink.insert(make_pair(linkEnd2,linkWeight));
      Links.insert(make_pair(linkEnd1,toLink));
      Nlinks++;
    }
    else{
      map<int,double>::iterator toLink_it = fromLink_it->second.find(linkEnd2);
      if(toLink_it == fromLink_it->second.end()){ // new link
        fromLink_it->second.insert(make_pair(linkEnd2,linkWeight));
        Nlinks++;
      }
      else{
        toLink_it->second += linkWeight;
        NdoubleLinks++;
      }
    }
  }
  
  net.close();
  
  cout << "done! (found " << Nnode << " nodes and " << Nlinks << " links";
  if(NdoubleLinks > 0)
    cout << ", aggregated " << NdoubleLinks << " link(s) defined more than once";
  
  /////////// Partition network /////////////////////
  double totalDegree = 0.0;
  vector<double> degree(Nnode,0.0);
  Node **node = new Node*[Nnode];
  for(int i=0;i<Nnode;i++){
    node[i] = new Node(i);
  }
  
  int NselfLinks = 0;
  for(map<int,map<int,double> >::iterator fromLink_it = Links.begin(); fromLink_it != Links.end(); fromLink_it++){
    for(map<int,double>::iterator toLink_it = fromLink_it->second.begin(); toLink_it != fromLink_it->second.end(); toLink_it++){
      
      int from = fromLink_it->first;
      int to = toLink_it->first;
      double weight = toLink_it->second;
      if(weight > 0.0){
        if(from == to){
          NselfLinks++;
        }
        else{
          node[from]->links.push_back(make_pair(to,weight));
          node[to]->links.push_back(make_pair(from,weight));
          node[from]->degree += weight;
          node[to]->degree += weight;
          totalDegree += 2*weight;
          degree[from] += weight;
          degree[to] += weight;
        }
      }
    }
  }
  if(NselfLinks > 0)
    cout << ", ignoring " <<  NselfLinks << " self link(s)." << endl;
  else
    cout << ")" << endl;
  
  //Swap maps to free memory
  for(map<int,map<int,double> >::iterator it = Links.begin(); it != Links.end(); it++)
    map<int,double>().swap(it->second);
  map<int,map<int,double> >().swap(Links);
  
  // Calculate uncompressed code length
  double uncompressedCodeLength = 0.0;
  for(int i=0;i<Nnode;i++){
    double p = degree[i]/totalDegree;
    if(p > 0.0)
      uncompressedCodeLength -= p*log(p)/log(2.0);
  }
  
  // Partition network hierarchically
  treeNode map;
  treeStats stats;
  double codeLength = repeated_hierarchical_partition(networkName,degree,nodeNames,R,node,map,totalDegree,Nnode,Ntrials,recursive,stats);
  
  cout << endl << "Best codelength = " << codeLength << " bits." << endl;
  cout << "Compression: " << 100.0*(1.0-codeLength/uncompressedCodeLength) << " percent." << endl;
  cout << "Number of modules: " << stats.Nmodules << endl;
  cout << "Number of large modules (> 1 percent of total number of nodes): " << stats.NlargeModules << endl;
  cout << "Average depth: " << stats.aveDepth << endl;
  cout << "Average size: " << stats.aveSize << endl;
  cout << "Gain over two-level code: " << 100.0*(stats.twoLevelCodeLength-codeLength)/codeLength << " percent." << endl;
  
  for(int i=0;i<Nnode;i++)
    delete node[i];
  delete [] node;
  delete R;
  
}

double fast_hierarchical_partition(MTRand *R, Node **orig_node, treeNode &map, double totalDegree, int Nnode, double &twoLevelCodeLength, bool deep){
  
  //MEMBERS FASTER WITH VECTOR?
  
  // Construct sub network
  int sub_Nnode = map.members.size();
  Node **sub_node = new Node*[sub_Nnode];
  genSubNet(orig_node,Nnode,sub_node,sub_Nnode,map,totalDegree);
  
  if(sub_Nnode == 1){
    // Clean up
    for(int i=0;i<sub_Nnode;i++)
      delete sub_node[i];
    delete [] sub_node;
    
    return map.codeLength;
  }
  
  // Store best map
  double best_codeLength = map.codeLength;
  treeNode best_map = map;
  
  // Initiate solver
  GreedyBase* sub_greedy;
  sub_greedy = new Greedy(R,sub_Nnode,totalDegree,sub_node);
  sub_greedy->initiate();
  
  // If a subtree exists, use this information
  if(!map.nextLevel.empty()){ 
    vector<int> cluster = vector<int>(sub_Nnode);
    int sub_Nmod = 0;
    for(multimap<double,treeNode,greater<double> >::iterator subsub_it = map.nextLevel.begin(); subsub_it != map.nextLevel.end(); subsub_it++){
      set<int> subsub_members = subsub_it->second.members;
      for(set<int>::iterator mem = subsub_members.begin(); mem != subsub_members.end(); mem++)
        cluster[map.renumber[(*mem)]] = sub_Nmod;
      sub_Nmod++;
    }
    sub_greedy->determMove(cluster);
  }
  
  if(map.level == 1)
    cout << "Partition first level into..." << flush;
  
  // Partition partition --> index codebook + module codebook
  partition(R,&sub_node,sub_greedy,true);
  
  if(map.level == 1){
    cout << sub_greedy->Nnode << " module(s) with code length " << sub_greedy->codeLength << " bits." << endl;
    if(sub_greedy->codeLength < twoLevelCodeLength)
      twoLevelCodeLength = sub_greedy->codeLength;
  }
  
  double subIndexLength = sub_greedy->indexLength;
  double subCodeLength = sub_greedy->codeLength;
  
  // Continue only if the network splits
  if(sub_greedy->Nnode == 1){
    
    // Clean up
    for(int i=0;i<sub_greedy->Nnode;i++)
      delete sub_node[i];
    delete [] sub_node;
    delete sub_greedy;
    
    return map.codeLength;
    
  }
  else{
    
    // Store temporary result
    vector<vector<int> > members = vector<vector<int> >(sub_greedy->Nnode);
    for (int i=0; i<sub_greedy->Nnode; i++){
      int Nmembers = sub_node[i]->members.size();
      members[i] = vector<int>(Nmembers);
      for(int j=0; j<Nmembers; j++){
        members[i][j] = map.rev_renumber[sub_node[i]->members[j]];
      }
    }
    
    // Clear deeper hierarchy of previous map
    multimap<double,treeNode,greater<double> >().swap(map.nextLevel);
    
    // Extend map with the new level
    int Nm = members.size();
    for(int i=0;i<Nm;i++){
      treeNode sub_map;
      double degree = 0.0;
      int Nn = members[i].size();
      for(int j=0;j<Nn;j++){
        int member = members[i][j];
        degree += orig_node[member]->degree;
        sub_map.members.insert(member);
      }
      sub_map.level = map.level + 1;
      setCodeLength(orig_node,sub_node,i,sub_map,totalDegree);
      
      map.nextLevel.insert(make_pair(degree/totalDegree,sub_map));
    }
    
    map.codeLength = subIndexLength;
    best_codeLength = subCodeLength;
    best_map = map;
    
    
    if(!deep){
      
      // Clean up
      for(int i=0;i<sub_greedy->Nnode;i++)
        delete sub_node[i];
      delete [] sub_node;
      delete sub_greedy;
      
      return subCodeLength;
      
    }
    
    double codeLength = map.codeLength;
    
    // Add index codebooks as long as the code gets shorter
    int Nmod = 2*sub_greedy->Nnode;
    while(sub_greedy->Nnode > 1 && sub_greedy->Nnode != Nmod){ // Continue as long as the network can be partitioned and the result is non-trivial
      
      // Add index codebook <--> move up in hierarchy
      Nmod = sub_greedy->Nnode; 
      for(int i=0; i<sub_greedy->Nnode; i++)
        vector<int>(1,i).swap(sub_node[i]->members);
      sub_greedy->initiate();
    
      if(map.level == 1)
        cout << "Trying to add index codebook..." << flush;
        
      partition(R,&sub_node,sub_greedy,true);
          
      map.codeLength = sub_greedy->moduleLength;
      subIndexLength = sub_greedy->moduleLength;  // Because the module has been collapsed
      
      // If trivial result
      if(sub_greedy->Nnode > 1 && sub_greedy->Nnode != Nmod){
        
        // Store temporary result
        vector<vector<int> > newMembers = vector<vector<int> >(sub_greedy->Nnode);
        for(int i=0; i<sub_greedy->Nnode; i++){
          int Nmembers = sub_node[i]->members.size();
          for(int j=0; j<Nmembers; j++){
            copy(members[sub_node[i]->members[j]].begin(),members[sub_node[i]->members[j]].end(),back_inserter(newMembers[i]));
          }
        }
        vector<vector<int> >(newMembers).swap(members);
        
        // Delete subtrees before generating new ones
        multimap<double,treeNode,greater<double> >().swap(map.nextLevel);
        
        // Store result as one-level subtrees
        int Nm = members.size();
        for(int i=0;i<Nm;i++){
          treeNode sub_map;
          double degree = 0.0;
          int Nn = members[i].size();
          for(int j=0;j<Nn;j++){
            int member = members[i][j];
            degree += orig_node[member]->degree;
            sub_map.members.insert(member);
          }
          sub_map.level = map.level + 1;
          map.nextLevel.insert(make_pair(degree/totalDegree,sub_map));
        }
        
        if(map.level == 1)
          cout << "succeeded. " << sub_greedy->Nnode << " modules with estimated code length... " << flush;
        
        // Create hierarchical tree under current level recursively 
        codeLength = sub_greedy->indexLength;
      
        for(multimap<double,treeNode,greater<double> >::iterator it = map.nextLevel.begin(); it != map.nextLevel.end(); it++)
          codeLength += fast_hierarchical_partition(R,orig_node,it->second,totalDegree,Nnode,twoLevelCodeLength,false);
      
        if(map.level == 1)
          cout << codeLength << " bits." << endl;
        
        if(codeLength < best_codeLength - 1.0e-10) { // Improvement
          
          map.codeLength = sub_greedy->indexLength;
          subIndexLength = sub_greedy->indexLength;
          best_codeLength = codeLength;
          best_map = map;
          
        }
        else{ // Longer code, restore best result and stop
        
          if(map.level == 1)
            cout << "No improvement. Now trying to create deeper structure based on best result." << endl;
          
          break;
        }
      
      }
      else {
        if(map.level == 1)
          cout << "failed. Now trying to create deeper structure based on best result." << endl;
      }
      
    }   
    
    // Clean up
    for(int i=0;i<sub_greedy->Nnode;i++)
      delete sub_node[i];
    delete [] sub_node;
    delete sub_greedy;   
    
    // Restore best map
    map = best_map;      
    
  }
  
  // Create hierarchical tree under current level recursively 
  double codeLength = map.codeLength;
  for(multimap<double,treeNode,greater<double> >::iterator it = map.nextLevel.begin(); it != map.nextLevel.end(); it++)
    codeLength += fast_hierarchical_partition(R,orig_node,it->second,totalDegree,Nnode,twoLevelCodeLength,true);
  
  // Update best map if improvements
  if(codeLength < best_codeLength - 1.0e-10){
    best_codeLength = codeLength;
  }
  else {
    map = best_map;
  }
  
  return best_codeLength;
  
}


double hierarchical_partition(MTRand *R, Node **orig_node, treeNode &map, double totalDegree, int Nnode, double recursive){
  
  //MEMBERS FASTER WITH VECTOR?
  
  // Store best map
  treeNode return_map;
  double return_codeLength = 1000.0;
  bool improvement = true;
  
  
  do{
    
    // Construct sub network
    int sub_Nnode = map.members.size();
    Node **sub_node = new Node*[sub_Nnode];
    genSubNet(orig_node,Nnode,sub_node,sub_Nnode,map,totalDegree);
    
    if(sub_Nnode == 1){
      // Clean up
      for(int i=0;i<sub_Nnode;i++)
        delete sub_node[i];
      delete [] sub_node;
      
      return map.codeLength;
    }
    
    // Store best map
    double best_codeLength = map.codeLength;
    treeNode best_map = map;
    
    // Initiate solver
    GreedyBase* sub_greedy;
    sub_greedy = new Greedy(R,sub_Nnode,totalDegree,sub_node);
    sub_greedy->initiate();
    
    // If a subtree exists, use this information
    if(!map.nextLevel.empty()){ 
      vector<int> cluster = vector<int>(sub_Nnode);
      int sub_Nmod = 0;
      for(multimap<double,treeNode,greater<double> >::iterator subsub_it = map.nextLevel.begin(); subsub_it != map.nextLevel.end(); subsub_it++){
        set<int> subsub_members = subsub_it->second.members;
        for(set<int>::iterator mem = subsub_members.begin(); mem != subsub_members.end(); mem++)
          cluster[map.renumber[(*mem)]] = sub_Nmod;
        sub_Nmod++;
      }
      sub_greedy->determMove(cluster);
    }
    
    // Partition partition --> index codebook + module codebook
    partition(R,&sub_node,sub_greedy,true);
    
    double subIndexLength = sub_greedy->indexLength;
    double subCodeLength = sub_greedy->codeLength;
    
    // Continue only if the network splits
    if(sub_greedy->Nnode == 1){
      
      // Clean up
      for(int i=0;i<sub_greedy->Nnode;i++)
        delete sub_node[i];
      delete [] sub_node;
      delete sub_greedy;
      
      return map.codeLength;
      
    }
    else{
      
      // Store temporary result
      vector<vector<int> > members = vector<vector<int> >(sub_greedy->Nnode);
      for (int i=0; i<sub_greedy->Nnode; i++){
        int Nmembers = sub_node[i]->members.size();
        members[i] = vector<int>(Nmembers);
        for(int j=0; j<Nmembers; j++){
          members[i][j] = map.rev_renumber[sub_node[i]->members[j]];
        }
      }
      
      // Clear deeper hierarchy of previous map
      multimap<double,treeNode,greater<double> >().swap(map.nextLevel);
      
      // Extend map with the new level
      int Nm = members.size();
      for(int i=0;i<Nm;i++){
        treeNode sub_map;
        double degree = 0.0;
        int Nn = members[i].size();
        for(int j=0;j<Nn;j++){
          int member = members[i][j];
          degree += orig_node[member]->degree;
          sub_map.members.insert(member);
        }
        sub_map.level = map.level + 1;
        setCodeLength(orig_node,sub_node,i,sub_map,totalDegree);
        
        map.nextLevel.insert(make_pair(degree/totalDegree,sub_map));
      }
      
      // Shouldn't be necessary
      if(subCodeLength < best_codeLength - 1.0e-10) { // Improvement
        best_codeLength = subCodeLength;
        best_map = map;
      }
      
      double codeLength = map.codeLength;
      
      if(R->rand() < recursive){ // Skipping saves time
        
        // Create hierarchical tree under current level recursively 
        map.codeLength = subIndexLength;
        codeLength = subIndexLength;
        for(multimap<double,treeNode,greater<double> >::iterator it = map.nextLevel.begin(); it != map.nextLevel.end(); it++){
          codeLength += hierarchical_partition(R,orig_node,it->second,totalDegree,Nnode,recursive);
        }
        
        // Update best map if improvements
        if(codeLength < best_codeLength){
          best_codeLength = codeLength;
          best_map = map;
        }
        
      }
      
      // Add index codebooks as long as the code gets shorter
      int Nmod = 2*sub_greedy->Nnode;
      while(sub_greedy->Nnode > 1 && sub_greedy->Nnode != Nmod){ // Continue as long as the network can be partitioned and the result is non-trivial
        
        // Add index codebook <--> move up in hierarchy
        Nmod = sub_greedy->Nnode; 
        for(int i=0; i<sub_greedy->Nnode; i++)
          vector<int>(1,i).swap(sub_node[i]->members);
        sub_greedy->initiate();
        partition(R,&sub_node,sub_greedy,true);
        
        map.codeLength = sub_greedy->moduleLength;
        subIndexLength = sub_greedy->moduleLength;  // Because the module has been collapsed
        
        // If trivial result
        if(sub_greedy->Nnode > 1 && sub_greedy->Nnode != Nmod){
          
          // Store temporary result
          vector<vector<int> > newMembers = vector<vector<int> >(sub_greedy->Nnode);
          for(int i=0; i<sub_greedy->Nnode; i++){
            int Nmembers = sub_node[i]->members.size();
            for(int j=0; j<Nmembers; j++){
              copy(members[sub_node[i]->members[j]].begin(),members[sub_node[i]->members[j]].end(),back_inserter(newMembers[i]));
            }
          }
          vector<vector<int> >(newMembers).swap(members);
          
          // Delete subtrees before generating new ones
          multimap<double,treeNode,greater<double> >().swap(map.nextLevel);
          
          // Store result as one-level subtrees
          int Nm = members.size();
          for(int i=0;i<Nm;i++){
            treeNode sub_map;
            double degree = 0.0;
            int Nn = members[i].size();
            for(int j=0;j<Nn;j++){
              int member = members[i][j];
              degree += orig_node[member]->degree;
              sub_map.members.insert(member);
            }
            sub_map.level = map.level + 1;
            map.nextLevel.insert(make_pair(degree/totalDegree,sub_map));
          }
          
          // Create hierarchical tree under current level recursively 
          codeLength = sub_greedy->indexLength;
          
          for(multimap<double,treeNode,greater<double> >::iterator it = map.nextLevel.begin(); it != map.nextLevel.end(); it++)
            codeLength += hierarchical_partition(R,orig_node,it->second,totalDegree,Nnode,recursive);
          if(codeLength < best_codeLength - 1.0e-10) { // Improvement
            
            best_codeLength = codeLength;
            best_map = map;
            
          }
          else{ // Longer code, restore best result and stop
            
            break;
            
          }
          
        }
        
      }   
      
      // Clean up
      for(int i=0;i<sub_greedy->Nnode;i++)
        delete sub_node[i];
      delete [] sub_node;
      delete sub_greedy;   
      
      // Restore best map
      map = best_map;      
      
    }
    
    if(map.level == 1 && return_codeLength > 100.0 && recursive > 0.0)
      cout << "Tuning" << flush;
    
    if(best_codeLength < return_codeLength - 1.0e-10){
      // Save best map
      return_map = best_map;
      return_codeLength = best_codeLength;
      
    }
    else {
      improvement = false;
    }
    
    if(map.level == 1 && recursive > 0.0)
      cout << "." << flush;
    
  } while(improvement && R->rand() < recursive);
  
  if(map.level == 1 && recursive > 0.0)
    cout << endl;
  
  // Restore best map
  map = return_map;
  
  return return_codeLength;
  
}

double repeated_hierarchical_partition(string networkName,vector<double> &degree, vector<string> &nodeNames, MTRand *R, Node **orig_node, treeNode &best_map, double totalDegree, int Nnode,int Ntrials,double recursive, treeStats &stats){
  
  double shortestCodeLength = 1000.0;
  stats.twoLevelCodeLength = 1000.0;
  
  for(int trial = 0; trial<Ntrials;trial++){
    
    cout << "Attempt " << trial+1 << "/" << Ntrials << ":" << endl;
    treeNode map;
    map.level = 1;
    //map.stop = false;
    for(int i=0;i<Nnode;i++){
      map.members.insert(i);
    }
    
    //double codeLength = hierarchical_partition(R,orig_node,map,totalDegree,Nnode,recursive);
    double codeLength = fast_hierarchical_partition(R,orig_node,map,totalDegree,Nnode,stats.twoLevelCodeLength,true);   
    
    cout << "Code length = " << codeLength << " bits." << endl;
    
    if(codeLength < shortestCodeLength){
      
      shortestCodeLength = codeLength;
      best_map = map;
      
      //Print hierarchical partition
      ostringstream oss;
      oss << networkName << ".tree";
      cout << endl << "New best result. Writing hierarchy to " << networkName << ".tree ... " << flush; 
      ofstream outfile;
      outfile.open(oss.str().c_str());
      outfile << "# Codelength = " << codeLength << " bits." << endl;
      string s;
      int depth=1;
      stats.aveDepth = 0.0;
      stats.aveSize = 0.0;
      stats.Nmodules = 0;
      stats.NlargeModules = 0;
      stats.largeModuleLimit = static_cast<int>(0.01*Nnode);
      printTree(s,map,nodeNames,degree,totalDegree,&outfile,depth,stats);
      outfile.close();
      stats.aveDepth /= 1.0*Nnode;
			stats.aveSize /= 1.0*Nnode;
      cout << "done!" << endl;
      cout << "Average depth: " << stats.aveDepth << endl;
      cout << "Average size: " << stats.aveSize << endl;
      cout << "Number of modules: " << stats.Nmodules << endl;
      cout << "Number of large modules (> 1 percent of total number of nodes): " << stats.NlargeModules << endl;
      cout << "Gain over two-level code: " << 100.0*(stats.twoLevelCodeLength-codeLength)/codeLength << " percent." << endl;

      vector<double> size(Nnode);
      for(int i=0;i<Nnode;i++)
        size[i] = degree[i]/totalDegree;
      addNodesToMap(map,size);
      
      for(int level=0;level<=2;level++){ 
        // Print map in .map format for the Map Generator at www.mapequation.org
        // Collapse to two levels
        multimap<double,printTreeNode,greater<double> > collapsedmap;
        collapseTree(collapsedmap,map,size,level);
        int Nmod = collapsedmap.size();
        vector<int> cluster(Nnode,-1);
        int cluNr = 0;
        for(multimap<double,printTreeNode,greater<double> >::iterator it = collapsedmap.begin(); it != collapsedmap.end(); it++){
          it->second.rank = cluNr;
          for (multimap<double,int,greater<double> >::iterator mem = it->second.members.begin(); mem != it->second.members.end(); mem++) {
            cluster[mem->second] = cluNr;
          }
          cluNr++;
        }
        // Generate modular network
        int Nlinks = 0;
        multimap<int,multimap<int,double> > unsortedLinks;
        for(int i=0;i<Nnode;i++){
          int NLinks = orig_node[i]->links.size();
          for(int j=0;j<NLinks;j++){
            int from = cluster[i];
            int to = cluster[orig_node[i]->links[j].first];
            double linkFlow = orig_node[i]->links[j].second/totalDegree;
            if(from < to && from >= 0){
              Nlinks++;
              multimap<int,multimap<int,double> >::iterator fromLink_it = unsortedLinks.find(from);
              if(fromLink_it == unsortedLinks.end()){ // new link
                multimap<int,double> toLink;
                toLink.insert(make_pair(to,linkFlow));
                unsortedLinks.insert(make_pair(from,toLink));
              }
              else{
                multimap<int,double>::iterator toLink_it = fromLink_it->second.find(to);
                if(toLink_it == fromLink_it->second.end()){ // new link
                  fromLink_it->second.insert(make_pair(to,linkFlow));
                }
                else{
                  toLink_it->second += linkFlow;
                }
              }
            }
          }
        }
        // Order links by size
        vector<double> exit(Nmod,0.0);
        multimap<double,pair<int,int>,greater<double> > sortedLinks;
        for(multimap<int,multimap<int,double> >::iterator it = unsortedLinks.begin(); it != unsortedLinks.end(); it++){
          for(multimap<int,double>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++){
            int from = it->first;
            int to = it2->first;
            double linkFlow = it2->second;
            sortedLinks.insert(make_pair(linkFlow,make_pair(from+1,to+1)));
            exit[from] += linkFlow;
            exit[to] += linkFlow;
          }
        }
                
        // Print map in .map format for the Map Generator at www.mapequation.org
        oss.str("");
        oss << networkName << "_level" << level << ".map";
        outfile.open(oss.str().c_str());
        outfile << "# modules: " << Nmod << endl;
        outfile << "# modulelinks: " << sortedLinks.size() << endl;
        outfile << "# nodes: " << Nnode << endl;
        outfile << "# links: " << Nlinks << endl;
        outfile << "# codelength: " << codeLength << endl;
        outfile << "*Undirected" << endl;
        outfile << "*Modules " << Nmod << endl;
        for(multimap<double,printTreeNode,greater<double> >::iterator it = collapsedmap.begin(); it != collapsedmap.end(); it++){
          outfile << it->second.rank+1 << " \"" << nodeNames[it->second.members.begin()->second] << "\" " << it->second.size << " " << exit[it->second.rank] << endl;
        }
        int Nmem = 0;
        for(multimap<double,printTreeNode,greater<double> >::iterator it = collapsedmap.begin(); it != collapsedmap.end(); it++){
          Nmem += it->second.members.size();
        }
        outfile << "*Nodes " << Nmem << endl;
        for(multimap<double,printTreeNode,greater<double> >::iterator it = collapsedmap.begin(); it != collapsedmap.end(); it++){
          int k=1;
          for (multimap<double,int,greater<double> >::iterator mem = it->second.members.begin(); mem != it->second.members.end(); mem++) {
            outfile << it->second.rank+1 << ":" << k << " \"" << nodeNames[mem->second] << "\" " << mem->first << endl;
            k++;
          }
        }
        outfile << "*Links " << sortedLinks.size() << endl;
        for(multimap<double,pair<int,int>,greater<double> >::iterator it = sortedLinks.begin();it != sortedLinks.end();it++)   
          outfile << it->second.first << " " << it->second.second << " " << 1.0*it->first << endl;
        outfile.close();
      }
    }
    
  }
  
  return shortestCodeLength;
  
}

void partition(MTRand *R,Node ***node, GreedyBase *greedy, bool silent){
  
  int Nnode = greedy->Nnode;
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
      vector<int> subMoveTo = vector<int>(Nnode);
      vector<int> moveTo = vector<int>(Nnode);
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
          double totalDegree = 0.0;
          for(int j=0;j<sub_Nnode;j++){
            
            //    fprintf(stderr,"%d %d\n",j,(*it_mem));
            int orig_nr = (*it_mem);
            int orig_Nlinks = cpy_node[orig_nr]->links.size(); // ERROR HERE
            sub_renumber[orig_nr] = j;
            sub_rev_renumber[j] = orig_nr;
            sub_node[j] = new Node(j);
            for(int k=0;k<orig_Nlinks;k++){
              int orig_link = cpy_node[orig_nr]->links[k].first;
              int orig_link_newnr = sub_renumber[orig_link];
              double orig_weight = cpy_node[orig_nr]->links[k].second;
              if(orig_link < orig_nr){
                if(sub_mem.find(orig_link) != sub_mem.end()){
                  sub_node[j]->links.push_back(make_pair(orig_link_newnr,orig_weight));
                  sub_node[orig_link_newnr]->links.push_back(make_pair(j,orig_weight));
                  totalDegree += 2.0*orig_weight;
                }
              }
            }
            it_mem++;
          }
          
          GreedyBase* sub_greedy;
          sub_greedy = new Greedy(R,sub_Nnode,totalDegree,sub_node);
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
      greedy->node = rpt_node;
      greedy->initiate();
      greedy->determMove(subMoveTo);
      greedy->level(node,false); 
      greedy->determMove(moveTo);
      (*node) = rpt_node;
      
      outer_oldCodeLength = greedy->codeLength;
      
      if(!silent)
        cout << greedy->Nnode << " modules, looping " << flush;
      
    }
    else if(iteration > 0){
      
      if(!silent)
        cout << "Iteration " << iteration+1 << ", moving " << Nnode << " nodes, looping " << flush;
      
      
      Node **rpt_node = new Node*[Nnode];
      for(int i=0;i<Nnode;i++){
        rpt_node[i] = new Node();
        cpyNode(rpt_node[i],cpy_node[i]);
      }
      
      vector<int> moveTo = vector<int>(Nnode);
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
      greedy->node = rpt_node;
      greedy->initiate();
      greedy->determMove(moveTo);
      
      (*node) = rpt_node;
    }
    else{
      
      if(!silent)
        cout << "Iteration " << iteration+1 << ", moving " << Nnode << " nodes, looping " << flush;
      
    }
    
    double oldCodeLength;
    do{
      oldCodeLength = greedy->codeLength;
      bool moved = true;
      int Nloops = 0;
      int count = 0;
      while(moved){
        moved = false;
        double inner_oldCodeLength = greedy->codeLength;
        greedy->move(moved);
        Nloops++;
        count++;
        if(inner_oldCodeLength-greedy->codeLength < 1.0e-10)
          moved = false;
        
        if(count == 10){	  
          greedy->tune();
          count = 0;
        }
        // 	if(!silent){
        // 	  cerr << Nloops;
        // 	  int loopsize = to_string(Nloops).length();
        // 	  for(int i=0;i<loopsize;i++)
        // 	    cerr << "\b";
        // 	}
      }
      
      greedy->level(node,true);
      
      if(!silent)
        cout << Nloops << " " << flush;
      
    } while(oldCodeLength - greedy->codeLength >  1.0e-10);
    
    iteration++;
    if(!silent)
      cout << "times between mergings to code length " <<  greedy->codeLength << " in " << greedy->Nmod << " modules." << endl;
    
  } while(outer_oldCodeLength - greedy->codeLength > 1.0e-10);
  
  for(int i=0;i<Nnode;i++)
    delete cpy_node[i];
  delete [] cpy_node;
  
}

void repeated_partition(MTRand *R, Node ***node, GreedyBase *greedy, bool silent,int Ntrials){
  
  double shortestCodeLength = 1000.0;
  int Nnode = greedy->Nnode;
  vector<int> cluster = vector<int>(Nnode);
  
  for(int trial = 0; trial<Ntrials;trial++){
    
    if(!silent)
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
    
    if(greedy->codeLength < shortestCodeLength){
      
      shortestCodeLength = greedy->codeLength;
      
      // Store best partition
      for(int i=0;i<greedy->Nnode;i++){
        for(vector<int>::iterator mem = cpy_node[i]->members.begin(); mem != cpy_node[i]->members.end(); mem++){
          cluster[(*mem)] = i;
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
  
}

// Sloppy version - no accordion test

//double hierarchical_partition(MTRand *R, Node **orig_node, treeNode &map, double totalDegree, int Nnode, double recursive){
//  
//  // Construct sub network
//  int sub_Nnode = map.members.size();
//  Node **sub_node = new Node*[sub_Nnode];
//  genSubNet(orig_node,Nnode,sub_node,sub_Nnode,map,totalDegree);
//  
//  if(sub_Nnode == 1){
//    // Clean up
//    for(int i=0;i<sub_Nnode;i++)
//      delete sub_node[i];
//    delete [] sub_node;
//    
//    return map.codeLength;
//  }
//  
//  // Store best map
//  double best_codeLength = map.codeLength;
//  treeNode best_map = map;
//  
//  // Initiate solver
//  GreedyBase* sub_greedy;
//  sub_greedy = new Greedy(R,sub_Nnode,totalDegree,sub_node);
//  sub_greedy->initiate();
//  
//  // Partition partition --> index codebook + module codebook
//  partition(R,&sub_node,sub_greedy,true);
//  
//  double subIndexLength = sub_greedy->indexLength;
//  double subCodeLength = sub_greedy->codeLength;
//  
//  // Continue only if the network splits
//  if(sub_greedy->Nnode == 1){
//    
//    // Clean up
//    for(int i=0;i<sub_greedy->Nnode;i++)
//      delete sub_node[i];
//    delete [] sub_node;
//    delete sub_greedy;
//    
//    return map.codeLength;
//    
//  }
//  else{
//    
//    // Store temporary result
//    vector<vector<int> > members = vector<vector<int> >(sub_greedy->Nnode);
//    for (int i=0; i<sub_greedy->Nnode; i++){
//      int Nmembers = sub_node[i]->members.size();
//      members[i] = vector<int>(Nmembers);
//      for(int j=0; j<Nmembers; j++){
//        members[i][j] = map.rev_renumber[sub_node[i]->members[j]];
//      }
//    }
//    
//    // Clear deeper hierarchy of previous map
//    multimap<double,treeNode,greater<double> >().swap(map.nextLevel);
//    
//    // Extend map with the new level
//    int Nm = members.size();
//    for(int i=0;i<Nm;i++){
//      treeNode sub_map;
//      double degree = 0.0;
//      int Nn = members[i].size();
//      for(int j=0;j<Nn;j++){
//        int member = members[i][j];
//        degree += orig_node[member]->degree;
//        sub_map.members.insert(member);
//      }
//      sub_map.level = map.level + 1;
//      setCodeLength(orig_node,sub_node,i,sub_map,totalDegree);
//      
//      map.nextLevel.insert(make_pair(degree/totalDegree,sub_map));
//    }
//    
//    // Shouldn't be necessary
//    if(subCodeLength < best_codeLength - 1.0e-10) { // Improvement
//      best_codeLength = subCodeLength;
//      best_map = map;
//    }
//    
//    //Add index codebooks as long as the code gets shorter
//    int Nmod = 2*sub_greedy->Nnode;
//    int count = 0;
//    while(sub_greedy->Nnode > 1 && sub_greedy->Nnode != Nmod){ // Continue as long as the network can be partitioned and the result is non-trivial
//      
//      // Add index codebook <--> move up in hierarchy
//      Nmod = sub_greedy->Nnode; 
//      for(int i=0; i<sub_greedy->Nnode; i++)
//        vector<int>(1,i).swap(sub_node[i]->members);
//      sub_greedy->initiate();
//      partition(R,&sub_node,sub_greedy,true);
//      
//      map.codeLength = sub_greedy->moduleLength;
//      subIndexLength = sub_greedy->moduleLength;  // Because the module has collapsed
//      
//      // If trivial result
//      if(sub_greedy->Nnode > 1 && sub_greedy->Nnode != Nmod){
//        
//        // Store temporary result
//        vector<vector<int> > newMembers = vector<vector<int> >(sub_greedy->Nnode);
//        for(int i=0; i<sub_greedy->Nnode; i++){
//          int Nmembers = sub_node[i]->members.size();
//          for(int j=0; j<Nmembers; j++){
//            copy(members[sub_node[i]->members[j]].begin(),members[sub_node[i]->members[j]].end(),back_inserter(newMembers[i]));
//          }
//        }
//        vector<vector<int> >(newMembers).swap(members);
//        
//        // Delete subtrees before generating new ones
//        multimap<double,treeNode,greater<double> >().swap(map.nextLevel);
//        
//        // Store result as one-level subtrees
//        int Nm = members.size();
//        for(int i=0;i<Nm;i++){
//          treeNode sub_map;
//          double degree = 0.0;
//          int Nn = members[i].size();
//          for(int j=0;j<Nn;j++){
//            int member = members[i][j];
//            degree += orig_node[member]->degree;
//            sub_map.members.insert(member);
//          }
//          sub_map.level = map.level + 1;
//          map.nextLevel.insert(make_pair(degree/totalDegree,sub_map));
//        }
//        
//      }
//      
//    }   
//    
//    // Create hierarchical tree under current level recursively 
//    map.codeLength = subIndexLength;
//    double codeLength = subIndexLength;
//    for(multimap<double,treeNode,greater<double> >::iterator it = map.nextLevel.begin(); it != map.nextLevel.end(); it++){
//      codeLength += hierarchical_partition(R,orig_node,it->second,totalDegree,Nnode,recursive);
//    }
//    
//    // Update best map if improvements
//    if(codeLength < best_codeLength){
//      best_codeLength = codeLength;
//      best_map = map;
//    }
//    
//    // Clean up
//    for(int i=0;i<sub_greedy->Nnode;i++)
//      delete sub_node[i];
//    delete [] sub_node;
//    delete sub_greedy;   
//    
//    // Restore best map
//    map = best_map;      
//    
//    return best_codeLength;
//    
//  }
//  
//}

// no sloppy recursive < 0 step

//double hierarchical_partition(MTRand *R, Node **orig_node, treeNode &map, double totalDegree, int Nnode, double recursive){
//  
//  //MEMBERS FASTER WITH VECTOR?
//  
//  // Store best map
//  treeNode return_map;
//  double return_codeLength = 1000.0;
//  bool improvement = true;
//  
//  do{
//    
//    // Construct sub network
//    int sub_Nnode = map.members.size();
//    Node **sub_node = new Node*[sub_Nnode];
//    genSubNet(orig_node,Nnode,sub_node,sub_Nnode,map,totalDegree);
//    
//    if(sub_Nnode == 1){
//      // Clean up
//      for(int i=0;i<sub_Nnode;i++)
//        delete sub_node[i];
//      delete [] sub_node;
//      
//      return map.codeLength;
//    }
//    
//    // Store best map
//    double best_codeLength = map.codeLength;
//    treeNode best_map = map;
//    
//    // Initiate solver
//    GreedyBase* sub_greedy;
//    sub_greedy = new Greedy(R,sub_Nnode,totalDegree,sub_node);
//    sub_greedy->initiate();
//    
//    // If a subtree exists, use this information
//    if(!map.nextLevel.empty()){ 
//      vector<int> cluster = vector<int>(sub_Nnode);
//      int sub_Nmod = 0;
//      for(multimap<double,treeNode,greater<double> >::iterator subsub_it = map.nextLevel.begin(); subsub_it != map.nextLevel.end(); subsub_it++){
//        set<int> subsub_members = subsub_it->second.members;
//        for(set<int>::iterator mem = subsub_members.begin(); mem != subsub_members.end(); mem++)
//          cluster[map.renumber[(*mem)]] = sub_Nmod;
//        sub_Nmod++;
//      }
//      sub_greedy->determMove(cluster);
//    }
//    
//    // Partition partition --> index codebook + module codebook
//    partition(R,&sub_node,sub_greedy,true);
//    
//    double subIndexLength = sub_greedy->indexLength;
//    double subCodeLength = sub_greedy->codeLength;
//    
//    // Continue only if the network splits
//    if(sub_greedy->Nnode == 1){
//      
//      // Clean up
//      for(int i=0;i<sub_greedy->Nnode;i++)
//        delete sub_node[i];
//      delete [] sub_node;
//      delete sub_greedy;
//      
//      return map.codeLength;
//      
//    }
//    else{
//      
//      // Store temporary result
//      vector<vector<int> > members = vector<vector<int> >(sub_greedy->Nnode);
//      for (int i=0; i<sub_greedy->Nnode; i++){
//        int Nmembers = sub_node[i]->members.size();
//        members[i] = vector<int>(Nmembers);
//        for(int j=0; j<Nmembers; j++){
//          members[i][j] = map.rev_renumber[sub_node[i]->members[j]];
//        }
//      }
//      
//      // Clear deeper hierarchy of previous map
//      multimap<double,treeNode,greater<double> >().swap(map.nextLevel);
//      
//      // Extend map with the new level
//      int Nm = members.size();
//      for(int i=0;i<Nm;i++){
//        treeNode sub_map;
//        double degree = 0.0;
//        int Nn = members[i].size();
//        for(int j=0;j<Nn;j++){
//          int member = members[i][j];
//          degree += orig_node[member]->degree;
//          sub_map.members.insert(member);
//        }
//        sub_map.level = map.level + 1;
//        setCodeLength(orig_node,sub_node,i,sub_map,totalDegree);
//        
//        map.nextLevel.insert(make_pair(degree/totalDegree,sub_map));
//      }
//      
//      // Shouldn't be necessary
//      if(subCodeLength < best_codeLength - 1.0e-10) { // Improvement
//        best_codeLength = subCodeLength;
//        best_map = map;
//      }
//      
//      double codeLength = map.codeLength;
//      
//      if(R->rand() < recursive){
//      
//        // Create hierarchical tree under current level recursively 
//        map.codeLength = subIndexLength;
//        codeLength = subIndexLength;
//        for(multimap<double,treeNode,greater<double> >::iterator it = map.nextLevel.begin(); it != map.nextLevel.end(); it++){
//          codeLength += hierarchical_partition(R,orig_node,it->second,totalDegree,Nnode,recursive);
//        }
//        
//        // Update best map if improvements
//        if(codeLength < best_codeLength){
//          best_codeLength = codeLength;
//          best_map = map;
//        }
//        
//      }
//        
//      // Add index codebooks as long as the code gets shorter
//      int Nmod = 2*sub_greedy->Nnode;
//      int count = 0;
//      while(sub_greedy->Nnode > 1 && sub_greedy->Nnode != Nmod){ // Continue as long as the network can be partitioned and the result is non-trivial
//        
//        // Add index codebook <--> move up in hierarchy
//        Nmod = sub_greedy->Nnode; 
//        for(int i=0; i<sub_greedy->Nnode; i++)
//          vector<int>(1,i).swap(sub_node[i]->members);
//        sub_greedy->initiate();
//        partition(R,&sub_node,sub_greedy,true);
//        
//        map.codeLength = sub_greedy->moduleLength;
//        subIndexLength = sub_greedy->moduleLength;  // Because the module has been collapsed
//        
//        // If trivial result
//        if(sub_greedy->Nnode > 1 && sub_greedy->Nnode != Nmod){
//          
//          // Store temporary result
//          vector<vector<int> > newMembers = vector<vector<int> >(sub_greedy->Nnode);
//          for(int i=0; i<sub_greedy->Nnode; i++){
//            int Nmembers = sub_node[i]->members.size();
//            for(int j=0; j<Nmembers; j++){
//              copy(members[sub_node[i]->members[j]].begin(),members[sub_node[i]->members[j]].end(),back_inserter(newMembers[i]));
//            }
//          }
//          vector<vector<int> >(newMembers).swap(members);
//          
//          // Delete subtrees before generating new ones
//          multimap<double,treeNode,greater<double> >().swap(map.nextLevel);
//          
//          // Store result as one-level subtrees
//          int Nm = members.size();
//          for(int i=0;i<Nm;i++){
//            treeNode sub_map;
//            double degree = 0.0;
//            int Nn = members[i].size();
//            for(int j=0;j<Nn;j++){
//              int member = members[i][j];
//              degree += orig_node[member]->degree;
//              sub_map.members.insert(member);
//            }
//            sub_map.level = map.level + 1;
//            map.nextLevel.insert(make_pair(degree/totalDegree,sub_map));
//          }
//          
//          // Create hierarchical tree under current level recursively 
//          codeLength = sub_greedy->indexLength;
//          
//          for(multimap<double,treeNode,greater<double> >::iterator it = map.nextLevel.begin(); it != map.nextLevel.end(); it++)
//            codeLength += hierarchical_partition(R,orig_node,it->second,totalDegree,Nnode,recursive);
//          if(codeLength < best_codeLength - 1.0e-10) { // Improvement
//            
//            best_codeLength = codeLength;
//            best_map = map;
//            
//            if(map.level == 1)
//              cout << "Successfully added index codebook, codelength = " << codeLength << " bits." << endl;
//            count++;
//            
//          }
//          else{ // Longer code, restore best result and stop
//            
//            break;
//            
//          }
//          
//        }
//        
//      }   
//      
//      // Clean up
//      for(int i=0;i<sub_greedy->Nnode;i++)
//        delete sub_node[i];
//      delete [] sub_node;
//      delete sub_greedy;   
//      
//      // Restore best map
//      map = best_map;      
//      
//    }
//    
//    if(map.level == 1 && return_codeLength > 100.0 && recursive > 0.0)
//      cout << "Tuning" << flush;
//    
//    if(best_codeLength < return_codeLength - 1.0e-10){
//      // Save best map
//      return_map = best_map;
//      return_codeLength = best_codeLength;
//      
//    }
//    else {
//      improvement = false;
//    }
//    
//    if(map.level == 1 && recursive > 0.0)
//      cout << "." << flush;
//    
//  } while(improvement && R->rand() < recursive);
//  
//  if(map.level == 1 && recursive > 0.0)
//    cout << endl;
//  
//  // Restore best map
//  map = return_map;
//  
//  return return_codeLength;
//  
//}

