##CAG.jar

This folder includes a runnable jar file which computes different agreement measures, and 3 sets of example input files. It basically compares two groupings of the same set of datapoints and gives a score between -1/0 to 1 to their agreement, 1 being the highest. For the jar file use -h to see its parameters.

The groupings should be formatted as a simple text file where each line is one group, that lists nodes in that group. Network, if given, should be .pairs or .gml [other formats might work, not tested].

Examples:
### Disjoint Clusterings:
#### Only Structure Independent Clustering Agreement
		$./cag.jar examples/NMIexample_V.list examples/NMIexample_U1.list 
		ARI'     
		0.703    
	
#### Also Structure Dependent Variations
		$./cag.jar examples/NMIexample_V.list examples/NMIexample_U1.list -g examples/NMIexample.gml 
		ARI'     , AR_{x^2}^{\Sigma d}   , AR_{x^2}^{\xi}   
		0.703    , 0.762                 , 0.787            
	
#### Print All Agreement Variations
		$./cag.jar examples/NMIexample_V.list examples/NMIexample_U1.list -all  
		ARI'     , RI'   , Jaccard  , F_{\beta=1.0}  , VI   , NMI_{+}   
		0.703    , 0.851 , 0.64     , 0.78           , 0.842, 0.785 

### Overlapping Clusterings
####  Only Structure Independent Clustering Agreement
#####  Main Measure
			$./cag.jar examples/NMIexample_V.list examples/NMIexample_U1.list -o 
			ARI_\delta'    
			0.703    
#####  Print All Agreement Variations
			$./cag.jar examples/NMIexample_V.list examples/NMIexample_U1.list -all  -o 
			ARI_\delta'    , RI_\delta'  , I_{norm}'  , I_\sqrt{tr}'  , \omega  , A\omega  
			0.703          , 0.851       , 0.705      , 0.84          , 0.836   , 0.66  
#####  Also print nmi variations: 
			$./cag.jar examples/NMIexample_V.list examples/NMIexample_U1.list -all  -o +nmiVars 
			ARI_\delta'    , RI_\delta'  , I_{norm}'  , I_\sqrt{tr}'  , NMI'  , NMI''  , \omega  , A\omega  
			0.703          , 0.851       , 0.705      , 0.84          , 0.61  , 0.538  , 0.836   , 0.66 
#####   This one uses external code, which may need to be recompiled on the machine
	
####  Also Structure Dependent Variations
#####   Main Measure
		$./cag.jar examples/NMIexample_V.list examples/NMIexample_U1.list -g examples/NMIexample.gml -o 
		ARI_\delta'    , _{\bot}ARI_\delta'  , _{+}ARI_\delta'  
		0.703          , 0.806               , 0.793  
	
#####   Print All Variations
		$./cag.jar examples/NMIexample_V.list examples/NMIexample_U1.list -g examples/NMIexample.gml -all  -o 
		ARI_\delta'    , RI_\delta'  , I_{norm}'  , I_\sqrt{tr}'  , \omega  , A\omega  , _{\bot}ARI_\delta'  , _{\bot}RI_\delta'  , _{\bot}I_{norm}'  , _{\bot}I_\sqrt{tr}'  , _{+}ARI_\delta'  , _{+}RI_\delta'  , _{+}I_{norm}'  , _{+}I_\sqrt{tr}'  
		0.703          , 0.851       , 0.705      , 0.84          , 0.836   , 0.66     , 0.806               , 0.906              , 0.768             , 0.902                , 0.793            , 0.922           , 0.839          , 0.866             

		
		
###   Other Parameters:
#### print the details
		$./cag.jar NMIexample_V.list NMIexample_U1.list +details
		ARI' : AR_{x^2}^{|\cap|}     
		0.7031888798037612    
		
		$./cag.jar NMIexample_V.list NMIexample_U1.list +details -all
		ARI' : AR_{x^2}^{|\cap|}     , RI' : R_{x^2}^{|\cap|}   , Jaccard  , F_{\beta=1.0}       , VI : R_{xlogx}^{|\cap|}   , NMI_{+} : AR_{xlogx}^{|\cap|}   
		0.7031888798037612           , 0.8512396694214877       , 0.64     , 0.7804878048780487  , 0.8423282765538793        , 0.7847044314890324              
		       
#### choose between implementations
		 $./cag.jar NMIexample_V.list NMIexample_U1.list +details -etaBased
		ARI' : AR_{x^2}^{|\cap|}     
		0.7031888798037612    
		$./cag.jar NMIexample_V.list NMIexample_U1.list +details -algCont
		ARI_{o}'              
		0.7031888798037612    
		$./cag.jar NMIexample_V.list NMIexample_U1.list +details -algDelta
		ARI_\delta'           
		0.7031888798037612  
####  use the exact formulations which ignores the pairs of same nodes
		$./cag.jar NMIexample_V.list NMIexample_U1.list +details -algDelta -exact
		ARI_\delta            
		0.6597938144329897 
		$./cag.jar NMIexample_V.list NMIexample_U1.list +details -originals -exact
		ARI                   
		0.6597938144329897    
		
### A note on choosing the implementations) 
	The defaults are recommended implementations. 
	In case of disjoint clusters, the etaBased implementation is better than the originals. 
	Especially in case the clusters are non-covering. They are equivalent otherwise, except for ARI exact. 
	In case of overlapping, the algDelta is the only valid implementation. It also works for disjoint clusters.
	However since it is computationally more expensive, if clusters are disjoint, the etaBased is recommended.   
		
####  Different implementations are mostly interchangeable for disjoint covering clusters	
		$./cag.jar NMIexample_V.list NMIexample_U1.list +details -all -originals 
		ARI                   , Jaccard  , F_{\beta=1.0}       , VI                  , NMI_\Sigma          
		0.6597938144329897    , 0.64     , 0.7804878048780487  , 0.8423282765538793  , 0.7847044314890322 
		$./cag.jar NMIexample_V.list NMIexample_U1.list +details -all -exact -algDelta 
		ARI_\delta            , RI_\delta           , I_{norm}            , I_\sqrt{tr}         , Jaccard  , F_{\beta=1.0}       , VI^a                , NMI_{+}^a           
		0.6597938144329897    , 0.8363636363636364  , 0.6666666666666667  , 0.8395939690034018  , 0.64     , 0.7804878048780487  , 0.8423282765538793  , 0.7847044314890322  
		$./cag.jar NMIexample_V.list NMIexample_U1.list +details -all -exact -algCont
		ARI_{o}               , RI_{o}              , Jaccard  , F_{\beta=1.0}       , VI_{o}              , AVI_{o}             
		0.6597938144329897    , 0.8363636363636364  , 0.64     , 0.7804878048780487  , 0.8423282765538794  , 0.6691821805086787  
		$./cag.jar NMIexample_V.list NMIexample_U1.list +details -all  -algCont
		ARI_{o}'              , RI_{o}'             , Jaccard  , F_{\beta=1.0}       , VI_{o}'             , AVI_{o}'            
				0.7031888798037612    , 0.8512396694214877  , 0.64     , 0.7804878048780487  , 0.8423282765538794  , 0.7847044314890324 
		$./cag.jar NMIexample_V.list NMIexample_U1.list +details -all  -etaBased
		ARI' : AR_{x^2}^{|\cap|}     , RI' : R_{x^2}^{|\cap|}   , Jaccard  , F_{\beta=1.0}       , VI : R_{xlogx}^{|\cap|}   , NMI_{+} : AR_{xlogx}^{|\cap|}   
		0.7031888798037612           , 0.8512396694214877       , 0.64     , 0.7804878048780487  , 0.8423282765538793        , 0.7847044314890324              
		$./cag.jar NMIexample_V.list NMIexample_U1.list +details -all  -algDelta
		ARI_\delta'           , RI_\delta'          , I_{norm}'           , I_\sqrt{tr}'        , Jaccard  , F_{\beta=1.0}       , VI^a                , NMI_{+}^a           
		0.7031888798037612    , 0.8512396694214877  , 0.7047095797180738  , 0.8395939690034018  , 0.64     , 0.7804878048780487  , 0.8423282765538793  , 0.7847044314890322  
		
		
