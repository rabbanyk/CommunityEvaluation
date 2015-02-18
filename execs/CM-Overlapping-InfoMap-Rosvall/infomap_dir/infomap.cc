#include "infomap.h"

using namespace std;
using std::cout;
using std::cin;
using std::endl;

unsigned stou(char *s){
  return strtoul(s,(char **)NULL,10);
}

void printTree(string s,multimap<double,treeNode,greater<double> >::iterator it_tM,ofstream *outfile,bool flip);
void repeated_partition(MTRand *R, Node ***node, GreedyBase *greedy, bool silent,int Ntrials);
void partition(MTRand *R, Node ***node, GreedyBase *greedy, bool silent);

// Call: trade <seed> <Ntries>
int main(int argc,char *argv[]){
  
  if( argc < 4){
    cout << "Call: ./infomap <seed> <network.net> <# attempts> [selflinks]" << endl;
    exit(-1);
  }
  
  int Ntrials = atoi(argv[3]);  // Set number of partition attempts
  string line;
  string buf;
  
  MTRand *R = new MTRand(stou(argv[1]));

  string infile = string(argv[2]);
  string networkFile = string(argv[2]);
  string networkName(networkFile.begin(),networkFile.begin() + networkFile.find_last_of("."));
  string networkType(infile.begin() + infile.find_last_of("."),infile.end());

	bool includeSelfLinks = false;
	if(argc == 5)
		if(to_string(argv[4]) == "selflinks")
			includeSelfLinks = true;

  Network network(networkFile);
  
  if(networkType == ".net"){
    loadPajekNet(network);    
  }
  else{
    loadLinkList(network); 
  }

  int Nnode = network.Nnode;
  
  /////////// Partition network /////////////////////
  Node **node = new Node*[Nnode];
  for(int i=0;i<Nnode;i++){
    node[i] = new Node(i,network.nodeWeights[i]/network.totNodeWeights);
  }
  
  int NselfLinks = 0;
  for(map<pair<int,int>,double>::iterator it = network.Links.begin(); it != network.Links.end(); it++){
    
    int from = it->first.first;
    int to = it->first.second;
    double weight = it->second;
    if(weight > 0.0){
      if(from == to){
        NselfLinks++;
				if(includeSelfLinks)
					node[from]->selfLink += weight;
      }
      else{
        node[from]->outLinks.push_back(make_pair(to,weight));
        node[to]->inLinks.push_back(make_pair(from,weight));
      }
    }
  }
  
	if(includeSelfLinks)
  	cout << ", including " <<  NselfLinks << " self link(s)." << endl;	
	else
		cout << ", ignoring " <<  NselfLinks << " self link(s)." << endl;

  //Swap vector to free memory
  map<pair<int,int>,double>().swap(network.Links);
    
  // Initiation
  GreedyBase* greedy;
  greedy = new Greedy(R,Nnode,node,Nnode);
  greedy->initiate();
  
  vector<double> size(Nnode);
  for(int i=0;i<Nnode;i++)
    size[i] = node[i]->size;

  cout << "Now partition the network:" << endl;
  repeated_partition(R,&node,greedy,false,Ntrials);
  int Nmod = greedy->Nnode;
  cout << "Done! Code length " << greedy->codeLength/log(2.0) << " in " << Nmod << " modules." << endl; 
      
  // Order links by size
  vector<double> exit(Nmod,0.0);
  multimap<double,pair<int,int>,greater<double> > sortedLinks;
  for(int i=0;i<Nmod;i++){
    int NoutLinks = node[i]->outLinks.size();
    for(int j=0;j<NoutLinks;j++){
      double linkFlow = node[i]->outLinks[j].second/greedy->beta;
      sortedLinks.insert(make_pair(linkFlow,make_pair(i+1,node[i]->outLinks[j].first+1)));
      exit[i] += linkFlow;
    }
  }
  
  // Order modules by size
  multimap<double,treeNode,greater<double> > treeMap;
  multimap<double,treeNode,greater<double> >::iterator it_tM;
  for(int i=0;i<greedy->Nnode;i++){
    int Nmembers = node[i]->members.size();
    treeNode tmp_tN;
    it_tM = treeMap.insert(make_pair(node[i]->size,tmp_tN));
    it_tM->second.exit = exit[i];
    for(int j=0;j<Nmembers;j++){
      it_tM->second.members.insert(make_pair(size[node[i]->members[j]],make_pair(node[i]->members[j],network.nodeNames[node[i]->members[j]])));
    }
  }
  
  //Print partition in format "module:rank size name"
  ofstream outfile;
  ostringstream oss;
  oss << networkName << ".tree";
  outfile.open(oss.str().c_str());
  outfile << "# Code length " << greedy->codeLength/log(2.0) << " in " << Nmod << " modules." << endl; 
  int k = 1;
  for(multimap<double,treeNode,greater<double> >::iterator it = treeMap.begin(); it != treeMap.end(); it++){
    string s;
    s.append(to_string(k));
    s.append(":");
    printTree(s,it,&outfile,false);
    k++;
  }
  outfile.close();
  
  // Create cluster vector 
  vector<int> clusterVec = vector<int>(Nnode);
  int clusterNr = 0;  
  for(multimap<double,treeNode,greater<double> >::iterator mod = treeMap.begin(); mod != treeMap.end(); mod++){
	for(multimap<double,pair<int,string>,greater<double> >::iterator mem = mod->second.members.begin(); mem != mod->second.members.end(); mem++){
	  clusterVec[mem->second.first] = clusterNr;
	}
	clusterNr++;
  }
  // Print partition in Pajek's .clu format
  oss.str("");
  oss << networkName << ".clu";
  outfile.open(oss.str().c_str());
  outfile << "*Vertices " << Nnode << "\x0D\x0A";
  for(int i=0;i<Nnode;i++)
    outfile << clusterVec[i]+1 << "\x0D\x0A";
  outfile.close();
  
  // Print map in Pajek's .net format (links sorted in descending order)
  oss.str("");
  oss << networkName << "_map.net";
  outfile.open(oss.str().c_str());
  outfile << "*Vertices " << Nmod << "\x0D\x0A";
  for(int i=0;i<Nmod;i++)
    outfile << i+1 << " \"" << i+1 << "\"" << "\x0D\x0A";
  outfile << "*Arcs " << sortedLinks.size() << "\x0D\x0A";
  for(multimap<double,pair<int,int>,greater<double> >::iterator it = sortedLinks.begin();it != sortedLinks.end();it++)   
    outfile << "  " << it->second.first << " " << it->second.second << " " << it->first << "\x0D\x0A";
  outfile.close();
  
  // Print size of modules in Pajek's .vec format
  oss.str("");
  oss << networkName << "_map.vec";
  outfile.open(oss.str().c_str());
  outfile << "*Vertices " << Nmod << "\x0D\x0A";
  for(int i=0;i<Nmod;i++)
    outfile << node[i]->size << "\x0D\x0A";
  outfile.close();
  
  // Print map in .map format for the Map Generator at www.mapequation.org
  oss.str("");
  oss << networkName << ".map";
  outfile.open(oss.str().c_str());
  outfile << "# modules: " << Nmod << endl;
  outfile << "# modulelinks: " << sortedLinks.size() << endl;
  outfile << "# nodes: " << Nnode << endl;
  outfile << "# links: " << network.Nlinks << endl;
  outfile << "# codelength: " << greedy->codeLength/log(2.0) << endl;
  outfile << "*Directed" << endl;
  outfile << "*Modules " << Nmod << endl;
  k = 0;
  for(multimap<double,treeNode,greater<double> >::iterator it = treeMap.begin(); it != treeMap.end(); it++){
    outfile << k+1 << " \"" << it->second.members.begin()->second.second << "\" " << it->first << " " << it->second.exit << endl;
    k++;
  }
  outfile << "*Nodes " << Nnode << endl;
  k = 1;
  for(multimap<double,treeNode,greater<double> >::iterator it = treeMap.begin(); it != treeMap.end(); it++){
    string s;
    s.append(to_string(k));
    s.append(":");
    printTree(s,it,&outfile,true);
    k++;
  }
  outfile << "*Links " << sortedLinks.size() << endl;
  for(multimap<double,pair<int,int>,greater<double> >::iterator it = sortedLinks.begin();it != sortedLinks.end();it++)   
    outfile << it->second.first << " " << it->second.second << " " << 1.0*it->first << endl;
  outfile.close();
  
  //   // print size of vertices (imported as a vector in Pajek)
  //   strcpy(netname,"");
  //   sprintf(netname,"%s-%dmodules.vec",outfile1,greedy->Nmodules);
  //   ofw1 = fopen(netname,"w");
  //   fprintf(ofw1,"*Vertices %d\015\012",greedy->Nmodules);
  //   for(int i=0;i<greedy->Nmodules;i++)
  //       fprintf(ofw1,"%e \015\012",module[i]->prob - module[i]->exit);
  //   fclose(ofw1);
  
  
  for(int i=0;i<greedy->Nnode;i++){
    delete node[i];
  }
  delete [] node;
  
  delete greedy;
  delete R;
}

void partition(MTRand *R, Node ***node, GreedyBase *greedy, bool silent){
  
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
        cout << "Iteration " << iteration+1 << ", moving ";
      
      Node **rpt_node = new Node*[Nnode];
      for(int i=0;i<Nnode;i++){
        rpt_node[i] = new Node();
        cpyNode(rpt_node[i],cpy_node[i]);
      }
      vector<int> subMoveTo(Nnode);
      vector<int> moveTo(Nnode);
      int subModIndex = 0;
      
      for(int i=0;i<greedy->Nnode;i++){
        
        int sub_Nnode = (*node)[i]->members.size();
        
        if(sub_Nnode > 1){
          Node **sub_node = new Node*[sub_Nnode]; 
          set<int> sub_mem;
          for(int j=0;j<sub_Nnode;j++)
            sub_mem.insert((*node)[i]->members[j]);
          set<int>::iterator it_mem = sub_mem.begin();
          vector<int> sub_renumber = vector<int>(Nnode);
          vector<int> sub_rev_renumber = vector<int>(sub_Nnode);
          for(int j=0;j<sub_Nnode;j++){
            int orig_nr = (*it_mem);
            int orig_NoutLinks = cpy_node[orig_nr]->outLinks.size();
            int orig_NinLinks = cpy_node[orig_nr]->inLinks.size();
            sub_renumber[orig_nr] = j;
            sub_rev_renumber[j] = orig_nr;
            sub_node[j] = new Node(j,cpy_node[orig_nr]->teleportWeight/(*node)[i]->teleportWeight);
            sub_node[j]->selfLink =  cpy_node[orig_nr]->selfLink; // Take care of self-link
            for(int k=0;k<orig_NoutLinks;k++){
              int orig_link = cpy_node[orig_nr]->outLinks[k].first;
              int orig_link_newnr = sub_renumber[orig_link];
              double orig_weight = cpy_node[orig_nr]->outLinks[k].second;
              if(orig_link < orig_nr){
                if(sub_mem.find(orig_link) != sub_mem.end()){
                  sub_node[j]->outLinks.push_back(make_pair(orig_link_newnr,orig_weight));
                  sub_node[orig_link_newnr]->inLinks.push_back(make_pair(j,orig_weight));
                }
              }
            }
            for(int k=0;k<orig_NinLinks;k++){
              int orig_link = cpy_node[orig_nr]->inLinks[k].first;
              int orig_link_newnr = sub_renumber[orig_link];
              double orig_weight = cpy_node[orig_nr]->inLinks[k].second;
              if(orig_link < orig_nr){
                if(sub_mem.find(orig_link) != sub_mem.end()){
                  sub_node[j]->inLinks.push_back(make_pair(orig_link_newnr,orig_weight));
                  sub_node[orig_link_newnr]->outLinks.push_back(make_pair(j,orig_weight));
                }
              }
            }
            it_mem++;
          }
          
          GreedyBase* sub_greedy;
          sub_greedy = new Greedy(R,sub_Nnode,sub_node,sub_Nnode);
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
      greedy->Ndanglings = 0;
      greedy->node = rpt_node;
      greedy->calibrate();
      greedy->determMove(subMoveTo);
      greedy->level(node,false);
      greedy->determMove(moveTo);
      (*node) = rpt_node;
      
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
      
      vector<int>moveTo(Nnode);
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
      greedy->Ndanglings = 0;
      greedy->node = rpt_node;
      greedy->calibrate();
      greedy->determMove(moveTo);
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
      int count = 0;
      while(moved){
        moved = false;
        double inner_oldCodeLength = greedy->codeLength;
        greedy->move(moved);
        Nloops++;
        count++;
        if(fabs(greedy->codeLength - inner_oldCodeLength) < 1.0e-10)
          moved = false;
        
        if(count == 10){	  
          greedy->tune();
          count = 0;
        }
      }
      
      greedy->level(node,true);
      
      if(!silent)
        cout << Nloops << " ";
      
    } while(oldCodeLength - greedy->codeLength >  1.0e-10);
    
    iteration++;
    if(!silent)
      cout << "times between mergings to code length " <<  greedy->codeLength/log(2.0) << " in " << greedy->Nmod << " modules." << endl;
    
  } while(outer_oldCodeLength - greedy->codeLength > 1.0e-10);
  
  for(int i=0;i<Nnode;i++)
    delete cpy_node[i];
  delete [] cpy_node;
  
}

void repeated_partition(MTRand *R, Node ***node, GreedyBase *greedy, bool silent,int Ntrials){
  
  double shortestCodeLength = 1000.0;
  int Nnode = greedy->Nnode;
  vector<int> cluster(Nnode);
  
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
    greedy->Ndanglings = 0;
    greedy->node = cpy_node;
    greedy->calibrate();
    
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
  greedy->Ndanglings = 0;
  greedy->node = (*node);
  greedy->calibrate();
  greedy->determMove(cluster);
  greedy->level(node,true);
  
}


void printTree(string s,multimap<double,treeNode,greater<double> >::iterator it_tM,ofstream *outfile,bool flip){
  
  multimap<double,treeNode,greater<double> >::iterator it;
  if(it_tM->second.nextLevel.size() > 0){
    int i=1;
    for(it = it_tM->second.nextLevel.begin(); it != it_tM->second.nextLevel.end(); it++){
      string cpy_s(s + to_string(i) + ":");
      printTree(cpy_s,it,outfile,flip);
      i++;
    }
  }
  else{
    int i = 1;
    for(multimap<double,pair<int,string>,greater<double> >::iterator mem = it_tM->second.members.begin(); mem != it_tM->second.members.end(); mem++){
      if(flip){
        string cpy_s(s + to_string(i) + " \"" + mem->second.second + "\" " + to_string(mem->first));
        (*outfile) << cpy_s << endl;
      }
      else{
        string cpy_s(s + to_string(i) + " " + to_string(mem->first) + " \"" + mem->second.second + "\"");
        (*outfile) << cpy_s << endl;
      }
      i++;
    }  
  }
}

