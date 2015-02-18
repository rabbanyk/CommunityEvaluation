#include "Greedy.h"

Greedy::~Greedy(){  
  
  logFac = vector<double>(0);
  modWnode.clear();
  
}

Greedy::Greedy(MTRand *RR,int nnode,int nlinks,Node **ah){

  R = RR;
  Nnode = nnode;  
  Nlinks = nlinks;
  node = ah;
  Nmem = Nnode;
  Nmod = Nnode;
  pF = 0.0; // Penalty factor to obtain a solution with more links within than between modules (positive if not fulfilled).
}

void Greedy::move(bool &moved){
  
  if(Nnode > 1){
    
    //   Generate random enumeration of nodes
    vector<int> randomOrder(Nnode);
    for(int i=0;i<Nnode;i++)
      randomOrder[i] = i;
    for(int i=0;i<Nnode-1;i++){
      int randPos = i + R->randInt(Nnode-i-1);
      int tmp = randomOrder[i];
      randomOrder[i] = randomOrder[randPos];
      randomOrder[randPos] = tmp;
    }
    
  

    for(int k=0;k<Nnode;k++){
            
      // Pick nodes in random order
      int flip = randomOrder[k]; 
      
    // Create map with module links
      map<int,int> wNtoM;
      map<int,int>::iterator it_M;
      int nlinks = node[flip]->links.size();
      
      for(int j=0;j<nlinks;j++){
	int nb_M = node[node[flip]->links[j].first]->index;
	int nb_w = node[flip]->links[j].second;
	
	it_M = wNtoM.find(nb_M);
	if (it_M != wNtoM.end())
	  it_M->second += nb_w;
	else
	  wNtoM.insert(make_pair(nb_M,nb_w));
	
      }
      
      // Calculate exit weight to own module
      int fromM = node[flip]->index; // 
      int fromM_weight = 0;
      it_M = wNtoM.find(fromM);
      if (it_M != wNtoM.end()){
	fromM_weight = it_M->second;
      }
            
      int bestM = fromM;
      int best_weight = 0;
      int best_penalty = 0;
      double best_deltaScore = 0.0;
      double best_deltaModelLength = 0.0;
      double best_networkLength = 0.0;
      
      // Calculate change in description length when node is removed 
      double remove_networkLength = 0.0;
      int remove_penalty = 0;
      int remove_Nmem = mod_members[fromM] - node[flip]->members.size();
      int remove_inlinks = mod_inlinks[fromM] - node[flip]->inlinks - fromM_weight;
      // Change associated to the same module
      remove_networkLength += logChoose(remove_Nmem*(remove_Nmem-1)/2,remove_inlinks) - logChoose(mod_members[fromM]*(mod_members[fromM]-1)/2,mod_inlinks[fromM]);
      
      // Change associated to neighboring modules
      map<int,int>::iterator it_nodelink =  wNtoM.begin();
      if(it_nodelink != wNtoM.end() && it_nodelink->first == fromM)
	it_nodelink++;
  
      for(map<int,int>::iterator it_modulelink = mod_links[fromM].begin(); it_modulelink != mod_links[fromM].end(); it_modulelink++){
	int neighbor = it_modulelink->first;
	int remove_linkw = it_modulelink->second;
	
	if(it_nodelink != wNtoM.end()){
	  if(it_nodelink->first == neighbor){
	    remove_linkw -= it_nodelink->second; // New weight between neighbors
	    remove_penalty += theta(remove_linkw-mod_inlinks[neighbor]) - theta(it_modulelink->second-mod_inlinks[neighbor]); // Change at modules's neighbor
	    it_nodelink++;
	    if(it_nodelink != wNtoM.end() && it_nodelink->first == fromM)
	      it_nodelink++;
	    
	  }
	}
	remove_networkLength += logChoose(remove_Nmem*mod_members[neighbor],remove_linkw) - logChoose(mod_members[fromM]*mod_members[neighbor],it_modulelink->second);
	remove_penalty += theta(remove_linkw-remove_inlinks) - theta(it_modulelink->second-mod_inlinks[fromM]); // Change at module
      }
        
      // Find the move that minimizes the description length
      for(map<int,int>::iterator it_M = wNtoM.begin(); it_M != wNtoM.end(); it_M++){
	
	int toM = it_M->first;
	int wtoM = it_M->second;
	
	double delta_modelLength = 0.0;
	if(mod_members[fromM] == static_cast<int>(node[flip]->members.size())) // If move results in fewer modules
	  delta_modelLength = Nmem*log(1.0*(Nmod-1)/Nmod)/log(2.0)-Nmod*log(1.0*Nlinks)/log(2.0);
	else if(mod_members[toM] == 0) // If move results in more modules
	  delta_modelLength = Nmem*log(1.0*(Nmod+1)/Nmod)/log(2.0)+(Nmod+1)*log(1.0*Nlinks)/log(2.0);

	if(toM != fromM){
	  
	  // Calculate change in description length when node is removed 
	  double add_networkLength = 0.0;
	  int add_penalty = 0;
	  int add_Nmem = mod_members[toM] + node[flip]->members.size();
	  int add_inlinks = mod_inlinks[toM] + node[flip]->inlinks + wtoM;
	  // Change associated to the same module
	  add_networkLength += logChoose(add_Nmem*(add_Nmem-1)/2,add_inlinks) - logChoose(mod_members[toM]*(mod_members[toM]-1)/2,mod_inlinks[toM]);
	  
	  map<int,int>::iterator it_nodelink =  wNtoM.begin();
	  map<int,int>::iterator it_modulelink = mod_links[toM].begin();	
	  while(it_nodelink != wNtoM.end() || it_modulelink != mod_links[toM].end()){
	    
	    bool newLink = false;
	    bool oldLink = false;
	    
	    // Check if new module and node have link 
	    if(it_nodelink == wNtoM.end()){
	      oldLink = true;
	    }
	    else if(it_modulelink == mod_links[toM].end()){
	      newLink = true;
	    }
	    else{
	      if(it_nodelink->first == it_modulelink->first){
		oldLink = true;
		newLink = true;
	      }
	      else if(it_nodelink->first > it_modulelink->first){   
		oldLink = true;
	      }
	      else{
		newLink = true;
	      }
	    }
	    
	    if(oldLink & newLink){
	      int neighbor = it_modulelink->first;
	      int add_linkw = it_modulelink->second + it_nodelink->second;
	      if(neighbor != fromM){
		add_networkLength += logChoose(add_Nmem*mod_members[neighbor],add_linkw) - logChoose(mod_members[toM]*mod_members[neighbor],it_modulelink->second);
		add_penalty += theta(add_linkw-add_inlinks) - theta(it_modulelink->second-mod_inlinks[toM]); // Change at module
		add_penalty += theta(add_linkw-mod_inlinks[neighbor]) - theta(it_modulelink->second-mod_inlinks[neighbor]); // Change at modules's neighbor
	      }
	      else{ // Change associated to connection between fromM and toM
		add_linkw -= wtoM;
		add_networkLength += logChoose(add_Nmem*remove_Nmem,add_linkw) - logChoose(mod_members[toM]*remove_Nmem,add_linkw-fromM_weight);
		add_penalty += theta(add_linkw-remove_inlinks) - theta(add_linkw-fromM_weight-remove_inlinks); // Change at module
		add_penalty += theta(add_linkw-add_inlinks) - theta(add_linkw-fromM_weight-mod_inlinks[toM]); // Change at modules's neighbor
	      }
	      
	      it_nodelink++;
	      it_modulelink++;
	      
	    }
	    else if(oldLink){
	      int neighbor = it_modulelink->first;
	      int add_linkw = it_modulelink->second;
	      if(neighbor != fromM){
		add_networkLength += logChoose(add_Nmem*mod_members[neighbor],add_linkw) - logChoose(mod_members[toM]*mod_members[neighbor],it_modulelink->second);
		add_penalty += theta(add_linkw-add_inlinks) - theta(it_modulelink->second-mod_inlinks[toM]); // Change at module
		// No change at module's neighbor
	      }
	      else{ // Change associated to connection between fromM and toM
		add_linkw -= wtoM;
		add_networkLength += logChoose(add_Nmem*remove_Nmem,add_linkw) - logChoose(mod_members[toM]*remove_Nmem,add_linkw);
		add_penalty += theta(add_linkw+fromM_weight-remove_inlinks) - theta(add_linkw-remove_inlinks); // Change at module
		add_penalty += theta(add_linkw+fromM_weight-add_inlinks) - theta(add_linkw-mod_inlinks[toM]); // Change at modules's neighbor
	      }
	      it_modulelink++;
	    }
	    else{
	      
	      int neighbor = it_nodelink->first;
	      int add_linkw = it_nodelink->second;
	      if(neighbor != toM){ // New connection to new module's neighbor (and not connection to itself)
		add_networkLength += logChoose(add_Nmem*mod_members[neighbor],add_linkw);
		add_penalty += theta(add_linkw-add_inlinks); // Change at module
		add_penalty += theta(add_linkw-mod_inlinks[neighbor]); // Change at modules's neighbor 
	      }

	      it_nodelink++;
	    }
	    
	  }
	  
	  double deltaScore = delta_modelLength + remove_networkLength + pF*remove_penalty + add_networkLength + pF*add_penalty; 
	  
	  if(deltaScore < best_deltaScore){
	    bestM = toM;
	    best_weight = wtoM;
	    best_deltaScore = deltaScore;
	    best_deltaModelLength = delta_modelLength;
	    best_networkLength = remove_networkLength + add_networkLength;
	    best_penalty = remove_penalty + add_penalty;
	    
	  }
	  
	}
      }
       
      // Option to move node to empty module (if node not already alone)
      if(mod_members[fromM] > static_cast<int>(node[flip]->members.size())){
	if(Nempty > 0){
	  int toM = mod_empty[Nempty-1];
	  int wtoM = 0;
	  double delta_modelLength = Nmem*log(1.0*(Nmod+1)/Nmod)/log(2.0)+(Nmod+1)*log(1.0*Nlinks)/log(2.0);   //Move results in more modules
	  
	  // Calculate change in description length when node is removed 
	  double add_networkLength = 0.0;
	  int add_penalty = 0;
	  int add_Nmem = node[flip]->members.size();
	  int add_inlinks = node[flip]->inlinks;
	  // Change associated to the same module
	  add_networkLength += logChoose(add_Nmem*(add_Nmem-1)/2,add_inlinks);
	  
	  for(map<int,int>::iterator it_nodelink = wNtoM.begin(); it_nodelink != wNtoM.end(); it_nodelink++){	
            
	    int neighbor = it_nodelink->first;
	    int add_linkw = it_nodelink->second;
	    if(neighbor != fromM){
	      add_networkLength += logChoose(add_Nmem*mod_members[neighbor],add_linkw);
	      add_penalty += theta(add_linkw-add_inlinks); // Change at module
	      add_penalty += theta(add_linkw-mod_inlinks[neighbor]); // Change at modules's neighbor 
	    }
	    else{
	      add_networkLength += logChoose(add_Nmem*remove_Nmem,add_linkw);
	      add_penalty += theta(add_linkw-add_inlinks);
	      add_penalty += theta(add_linkw-remove_Nmem);
	    }
	  }
	  
	  double deltaScore = delta_modelLength + remove_networkLength + pF*remove_penalty + add_networkLength + pF*add_penalty; 

	  if(deltaScore < best_deltaScore){

	    bestM = toM;
	    best_weight = wtoM;
	    best_deltaScore = deltaScore;
	    best_deltaModelLength = delta_modelLength;
	    best_networkLength = remove_networkLength + add_networkLength;
	    best_penalty = remove_penalty + add_penalty;
	    
	  }
	}
      }

      // Make best possible move
      if(bestM != fromM){
	
	//Update empty module vector
	if(mod_members[bestM] == 0){
	  Nmod++;
	  Nempty--;
	}
	if(mod_members[fromM] == static_cast<int>(node[flip]->members.size())){
	  mod_empty[Nempty] = fromM;
	  Nmod--;
	  Nempty++;
	}
	
	modelLength += best_deltaModelLength;
	networkLength += best_networkLength;
	penalty += best_penalty;
	codeLength = modelLength + networkLength;
	score = codeLength + penalty;
	
	mod_inlinks[fromM] -= node[flip]->inlinks + fromM_weight;
	mod_members[fromM] -= node[flip]->members.size();
	mod_inlinks[bestM] += node[flip]->inlinks + best_weight;
	mod_members[bestM] += node[flip]->members.size();
	
	// Update links associated with fromM module 
	for(map<int,int>::iterator it_nodelink = wNtoM.begin(); it_nodelink != wNtoM.end(); it_nodelink++){
	  int neighbor = it_nodelink->first;
	  if(neighbor != fromM){ // Module is not connected to itself with a link
	    int weight = it_nodelink->second;
	    map<int,int>::iterator it_modulelink = mod_links[fromM].find(neighbor);
	    map<int,int>::iterator it_neighbormodulelink = mod_links[neighbor].find(fromM);
	    
	    if(it_modulelink->second == weight){ // Erase if the node has the only connections between the modules
	      mod_links[fromM].erase(it_modulelink);
	      mod_links[neighbor].erase(it_neighbormodulelink);
	    }
	    else{ // Otherwise update the number of connections
	      it_modulelink->second -= weight;
	      it_neighbormodulelink->second -= weight;
	    }
	  }
	}
	
	// Update links associated with bestM module 
	for(map<int,int>::iterator it_nodelink = wNtoM.begin(); it_nodelink != wNtoM.end(); it_nodelink++){
	  int neighbor = it_nodelink->first;
	  if(neighbor != bestM){ // Module is not connected to itself with a link
	    int weight = it_nodelink->second;
	    map<int,int>::iterator it_modulelink = mod_links[bestM].find(neighbor);
	    if(it_modulelink == mod_links[bestM].end()){ // Modules are not already connected
	      mod_links[bestM].insert(make_pair(neighbor,weight));
	      mod_links[neighbor].insert(make_pair(bestM,weight));   
	    }
	    else{ // Otherwise update the number of connections
	      it_modulelink->second += weight;  
	      map<int,int>::iterator it_neighbormodulelink = mod_links[neighbor].find(bestM);
	      it_neighbormodulelink->second += weight;
	    }
	  }
	}
	
	node[flip]->index = bestM;
	moved = true;
	
      }
      
    }
    
  }
}
  

void Greedy::initiate(void){
  
  genLogTable(static_cast<int>(Nnode*(Nnode-1)/2+10)); // generate look-up table
  
  calibrate();
  
}


void Greedy::calibrate(void){
  
  // !!! create mod_empty > Nmod ???
  mod_empty.clear();
  mod_empty = vector<int>(Nmod);
  Nempty = 0;
  
  int size = mod_links.size();
  for(int i=0;i<size;i++)
    mod_links[i].clear();
  mod_links.clear();
  mod_inlinks.clear();
  mod_links = vector<map<int,int> >(Nmod);
  mod_inlinks = vector<int>(Nmod);
  mod_members = vector<int>(Nmod);
  
  for(int i=0;i<Nmod;i++){
    
    int nlinks = node[i]->links.size();
    for(int j=0;j<nlinks;j++)
      mod_links[i].insert(make_pair(node[i]->links[j].first,node[i]->links[j].second));

    mod_inlinks[i] = node[i]->inlinks;

    mod_members[i] = node[i]->members.size();
    node[i]->index = i;
  }

  modelLength = Nmem*log(1.0*Nmod)/log(2.0) + 0.5*Nmod*(Nmod+1)*log(1.0*Nlinks)/log(2.0);
  
  networkLength = 0.0;
  penalty = 0;
  
  for(int i=0;i<Nmod;i++){
    
    for(map<int,int>::iterator it = mod_links[i].begin(); it != mod_links[i].end(); it++){

      if(i < it->first)
	networkLength += logChoose(mod_members[i]*mod_members[it->first],it->second);
      
      penalty += theta(it->second-mod_inlinks[i]);
    }
    
    networkLength += logChoose(mod_members[i]*(mod_members[i]-1)/2,mod_inlinks[i]);
    
  }
  
  codeLength = modelLength + networkLength; 
  
  score = codeLength + pF*penalty;

}

void Greedy::prepare(bool sort){

  Nmod = 0;
  modWnode.clear();
  
  if(sort){
    
    multimap<int,int> Msize;
    for(int i=0;i<Nnode;i++){
      if(mod_members[i] > 0){
	Nmod++;
	Msize.insert(make_pair(mod_members[i],i));
      }
    }
    
    for(multimap<int,int>::reverse_iterator it = Msize.rbegin(); it != Msize.rend(); it++)
      modWnode.push_back(it->second);
    
  }
  else{
    
    for(int i=0;i<Nnode;i++){
      if(mod_members[i] > 0){
	Nmod++;
	modWnode.push_back(i);
      }
    }
    
  }
  
}

void Greedy::level(Node ***node_tmp, bool sort){
  
  prepare(sort);
  
  (*node_tmp) = new Node*[Nmod];
  
  vector<int> nodeInMod = vector<int>(Nnode);
  for(int i=0;i<Nmod;i++){
    (*node_tmp)[i] = new Node();
    (*node_tmp)[i]->index = i;
    (*node_tmp)[i]->inlinks = mod_inlinks[modWnode[i]];
    nodeInMod[modWnode[i]] = i;
  }

  // Update links
  for(int i=0;i<Nmod;i++)
    for(map<int,int>::iterator it = mod_links[modWnode[i]].begin(); it != mod_links[modWnode[i]].end(); it++)
      (*node_tmp)[i]->links.push_back(make_pair(nodeInMod[it->first],it->second));
  
  // Update members
  for(int i=0;i<Nnode;i++)
    copy(node[i]->members.begin(),node[i]->members.end(),back_inserter((*node_tmp)[nodeInMod[node[i]->index]]->members));


  // Option to move to empty module
  mod_empty.clear();
  Nempty = 0;
  for(int i=0;i<Nnode;i++){
    delete node[i];
  }
  delete [] node;
  
  Nnode = Nmod;
  node = (*node_tmp);
  
  calibrate();
  
}

double Greedy::logChoose(int n,int k){
  
  return logFac[n] - logFac[k] - logFac[n-k];
  
}

void Greedy::genLogTable(int maxsize){
  
  logFac = vector<double>(maxsize);
  logFac[0] = 0.0;
  for(int i=1;i<maxsize;i++)
    logFac[i] = logFac[i-1]+log(1.0*i)/log(2.0);
  
}

int Greedy::theta(int pen){

  if(pen < 0)
    return 0;
  else
    return pen;

}

void Greedy::determMove(int *moveTo){

  for(int i=0;i<Nnode;i++){
    
    int fromM = i;
    int bestM = moveTo[i];
    
    if(fromM != bestM){
      
      // Create map with module links
      map<int,int> wNtoM;
      map<int,int>::iterator it_M;
      int nlinks = node[fromM]->links.size();
      

      for(int j=0;j<nlinks;j++){
	int nb_M = node[node[fromM]->links[j].first]->index;
	int nb_w = node[fromM]->links[j].second;
	
	it_M = wNtoM.find(nb_M);
	if (it_M != wNtoM.end())
	  it_M->second += nb_w;
	else
	  wNtoM.insert(make_pair(nb_M,nb_w));
      }
    
      // Calculate exit weight to own module
      int fromM_weight = 0;
      it_M = wNtoM.find(fromM);
      if (it_M != wNtoM.end()){
	fromM_weight = it_M->second;
      }
      
      // Calculate change in description length when node is removed from old module
      double remove_networkLength = 0.0;
      int remove_penalty = 0;
      int remove_Nmem = mod_members[fromM] - node[fromM]->members.size();
      int remove_inlinks = mod_inlinks[fromM] - node[fromM]->inlinks - fromM_weight;
      // Change associated to the same module
      remove_networkLength += logChoose(remove_Nmem*(remove_Nmem-1)/2,remove_inlinks) - logChoose(mod_members[fromM]*(mod_members[fromM]-1)/2,mod_inlinks[fromM]);
      
      // Change associated to neighboring modules
      map<int,int>::iterator it_nodelink =  wNtoM.begin();
      if(it_nodelink != wNtoM.end() && it_nodelink->first == fromM)
	it_nodelink++;

      for(map<int,int>::iterator it_modulelink = mod_links[fromM].begin(); it_modulelink != mod_links[fromM].end(); it_modulelink++){
	int neighbor = it_modulelink->first;
	int remove_linkw = it_modulelink->second;
	
	if(it_nodelink != wNtoM.end()){
	  if(it_nodelink->first == neighbor){
	    remove_linkw -= it_nodelink->second; // New weight between neighbors
	    remove_penalty += theta(remove_linkw-mod_inlinks[neighbor]) - theta(it_modulelink->second-mod_inlinks[neighbor]); // Change at modules's neighbor
	    it_nodelink++;
	    if(it_nodelink != wNtoM.end() && it_nodelink->first == fromM)
	      it_nodelink++;
	  }
	}
	remove_networkLength += logChoose(remove_Nmem*mod_members[neighbor],remove_linkw) - logChoose(mod_members[fromM]*mod_members[neighbor],it_modulelink->second);
	remove_penalty += theta(remove_linkw-remove_inlinks) - theta(it_modulelink->second-mod_inlinks[fromM]); // Change at module
      }

      int toM = bestM;
      int wtoM = 0;
      
      // Calculate change in description length when node is added to new module
      it_M = wNtoM.find(bestM);
      if(it_M != wNtoM.end())
	wtoM = it_M->second; 

      
      double delta_modelLength = 0.0;
      if(mod_members[fromM] == static_cast<int>(node[fromM]->members.size())) // If move results in fewer modules
	delta_modelLength = Nmem*log(1.0*(Nmod-1)/Nmod)/log(2.0)-Nmod*log(1.0*Nlinks)/log(2.0);
      else if(mod_members[toM] == 0) // If move results in more modules
	delta_modelLength = Nmem*log(1.0*(Nmod+1)/Nmod)/log(2.0)+(Nmod+1)*log(1.0*Nlinks)/log(2.0);

      // Calculate change in description length when node is removed 
      double add_networkLength = 0.0;
      int add_penalty = 0;
      int add_Nmem = mod_members[toM] + node[fromM]->members.size();
      int add_inlinks = mod_inlinks[toM] + node[fromM]->inlinks + wtoM;
      // Change associated to the same module
      add_networkLength += logChoose(add_Nmem*(add_Nmem-1)/2,add_inlinks) - logChoose(mod_members[toM]*(mod_members[toM]-1)/2,mod_inlinks[toM]);
      
      it_nodelink =  wNtoM.begin();
      map<int,int>::iterator it_modulelink = mod_links[toM].begin();	
      while(it_nodelink != wNtoM.end() || it_modulelink != mod_links[toM].end()){
	
	bool newLink = false;
	bool oldLink = false;
	
	// Check if new module and node have link 
	if(it_nodelink == wNtoM.end()){
	  oldLink = true;
	}
	else if(it_modulelink == mod_links[toM].end()){
	  newLink = true;
	}
	else{
	  if(it_nodelink->first == it_modulelink->first){
	    oldLink = true;
	    newLink = true;
	  }
	  else if(it_nodelink->first > it_modulelink->first){   
	    oldLink = true;
	    }
	  else{
	    newLink = true;
	  }
	}
	
	if(oldLink & newLink){
	  int neighbor = it_modulelink->first;
	  int add_linkw = it_modulelink->second + it_nodelink->second;
	  if(neighbor != fromM){
	    add_networkLength += logChoose(add_Nmem*mod_members[neighbor],add_linkw) - logChoose(mod_members[toM]*mod_members[neighbor],it_modulelink->second);
	    add_penalty += theta(add_linkw-add_inlinks) - theta(it_modulelink->second-mod_inlinks[toM]); // Change at module
	    add_penalty += theta(add_linkw-mod_inlinks[neighbor]) - theta(it_modulelink->second-mod_inlinks[neighbor]); // Change at modules's neighbor
	  }
	  else{ // Change associated to connection between fromM and toM
	    add_linkw -= wtoM;
	    add_networkLength += logChoose(add_Nmem*remove_Nmem,add_linkw) - logChoose(mod_members[toM]*remove_Nmem,add_linkw-fromM_weight);
	      add_penalty += theta(add_linkw-remove_inlinks) - theta(add_linkw-fromM_weight-remove_inlinks); // Change at module
	      add_penalty += theta(add_linkw-add_inlinks) - theta(add_linkw-fromM_weight-mod_inlinks[toM]); // Change at modules's neighbor
	  }
	  
	  it_nodelink++;
	  it_modulelink++;
	  
	}
	else if(oldLink){
	  int neighbor = it_modulelink->first;
	  int add_linkw = it_modulelink->second;
	  if(neighbor != fromM){
	    add_networkLength += logChoose(add_Nmem*mod_members[neighbor],add_linkw) - logChoose(mod_members[toM]*mod_members[neighbor],it_modulelink->second);
	    add_penalty += theta(add_linkw-add_inlinks) - theta(it_modulelink->second-mod_inlinks[toM]); // Change at module
	    // No change at module's neighbor
	  }
	  else{ // Change associated to connection between fromM and toM
	    add_linkw -= wtoM;
	    add_networkLength += logChoose(add_Nmem*remove_Nmem,add_linkw) - logChoose(mod_members[toM]*remove_Nmem,add_linkw);
	    add_penalty += theta(add_linkw+fromM_weight-remove_inlinks) - theta(add_linkw-remove_inlinks); // Change at module
	      add_penalty += theta(add_linkw+fromM_weight-add_inlinks) - theta(add_linkw-mod_inlinks[toM]); // Change at modules's neighbor
	  }
	  it_modulelink++;
	}
	else{
	  
	  int neighbor = it_nodelink->first;
	  int add_linkw = it_nodelink->second;
	  if(neighbor != toM){ // New connection to new module's neighbor (and not connection to itself)
	    add_networkLength += logChoose(add_Nmem*mod_members[neighbor],add_linkw);
	    add_penalty += theta(add_linkw-add_inlinks); // Change at module
	    add_penalty += theta(add_linkw-mod_inlinks[neighbor]); // Change at modules's neighbor 
	  }
	  else{
	    
	    
	  }
	  it_nodelink++;
	}	
      }
      
      int best_weight = wtoM;
      double best_deltaModelLength = delta_modelLength;
      double best_networkLength = remove_networkLength + add_networkLength;
      int best_penalty = remove_penalty + add_penalty;
            
      //Update empty module vector
      if(mod_members[bestM] == 0){
	Nmod++;
	Nempty--;
      }
      if(mod_members[fromM] == static_cast<int>(node[fromM]->members.size())){
 	mod_empty[Nempty] = fromM;
	Nmod--;
 	Nempty++;
      }
      
      modelLength += best_deltaModelLength;
      networkLength += best_networkLength;
      penalty += best_penalty;
      codeLength = modelLength + networkLength;
      score = codeLength + pF*penalty;

      mod_inlinks[fromM] -= node[fromM]->inlinks + fromM_weight;
      mod_members[fromM] -= node[fromM]->members.size();
      mod_inlinks[bestM] += node[fromM]->inlinks + best_weight;
      mod_members[bestM] += node[fromM]->members.size();
      
      // Update links associated with fromM module 
      for(map<int,int>::iterator it_nodelink = wNtoM.begin(); it_nodelink != wNtoM.end(); it_nodelink++){
	int neighbor = it_nodelink->first;
	if(neighbor != fromM){ // Module is not connected to itself with a link
	  int weight = it_nodelink->second;
	  map<int,int>::iterator it_modulelink = mod_links[fromM].find(neighbor);
	  map<int,int>::iterator it_neighbormodulelink = mod_links[neighbor].find(fromM);
	
	  if(it_modulelink->second == weight){ // Erase if the node has the only connections between the modules
	    mod_links[fromM].erase(it_modulelink);
	    mod_links[neighbor].erase(it_neighbormodulelink);
	  }
	  else{ // Otherwise update the number of connections
	    it_modulelink->second -= weight;
	    it_neighbormodulelink->second -= weight;
	  }
	}
      }

      // Update links associated with bestM module 
      for(map<int,int>::iterator it_nodelink = wNtoM.begin(); it_nodelink != wNtoM.end(); it_nodelink++){
	int neighbor = it_nodelink->first;
	if(neighbor != bestM){ // Module is not connected to itself with a link
	  int weight = it_nodelink->second;
	  map<int,int>::iterator it_modulelink = mod_links[bestM].find(neighbor);
	  if(it_modulelink == mod_links[bestM].end()){ // Modules are not already connected
	    mod_links[bestM].insert(make_pair(neighbor,weight));
	    mod_links[neighbor].insert(make_pair(bestM,weight));   
	  }
	  else{ // Otherwise update the number of connections
	    it_modulelink->second += weight;  
	    map<int,int>::iterator it_neighbormodulelink = mod_links[neighbor].find(bestM);
	    it_neighbormodulelink->second += weight;
	  }
	}
      }
      
      node[fromM]->index = bestM;

    }
  }
  
}
