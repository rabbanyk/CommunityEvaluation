// agmfast.cpp : Defines the entry point for the console application.
//
#include "stdafx.h"
#include "agmfast.h"
#include "agm.h"
#include "agmfit.h"
#include <omp.h>


int main(int argc, char* argv[]) {
	Env = TEnv(argc, argv, TNotify::StdNotify);
	Env.PrepArgs(TStr::Fmt("ragm. build: %s, %s. Time: %s", __TIME__, __DATE__, TExeTm::GetCurTm()));
	TExeTm ExeTm;
	Try
  TStr OutFPrx = Env.GetIfArgPrefixStr("-o:", "test.txt", "output file name for detected community affiliation");
	const TStr InFNm = Env.GetIfArgPrefixStr("-i:", "oregon1_010331.txt", "input edgelist file name");
  int OptComs = Env.GetIfArgPrefixInt("-c:", 100, "the number of communities to detect (-1:detect automatically)");
  const int MinComs = Env.GetIfArgPrefixInt("-mc:", 5, "minimum number of communities to try for cross-validation");
  const int MaxComs = Env.GetIfArgPrefixInt("-xc:", 100, "maximum number of communities to try for cross-validation");
  const int DivComs = Env.GetIfArgPrefixInt("-nc:", 5, "how many numbers to try in cross-validation for the number of communities");
  const int NumThreads = Env.GetIfArgPrefixInt("-nt:", 4, "number of threads for parallelization");
  const double StepAlpha = Env.GetIfArgPrefixFlt("-sa:", 0.05, "alpha for stepsize");
  const double StepBeta = Env.GetIfArgPrefixFlt("-sb:", 0.3, "beta for stepsize");

  omp_set_num_threads(NumThreads);
  PUNGraph G;
  if (InFNm.IsStrIn(".ungraph")) {
    TFIn GFIn(InFNm);
    G = TUNGraph::Load(GFIn);
  } else {
    PUNGraph FullG = TSnap::LoadEdgeList<PUNGraph>(InFNm);
    G = TSnap::GetMxWcc(FullG);
  }
  TSnap::DelSelfEdges(G);
	printf("Graph: %d Nodes %d Edges\n", G->GetNodes(), G->GetEdges());
  
	TVec<TIntV> EstCmtyVV;
	TExeTm RunTm;
	TAGMFast RAGM(G, 10, 10);
  
  if (OptComs == -1) {
    printf("finding number of communities\n");
    OptComs = RAGM.FindComsByCV(NumThreads, MaxComs, MinComs, DivComs, OutFPrx, StepAlpha, StepBeta);
  }

  RAGM.NeighborComInit(OptComs);
  if (NumThreads == 1 || G->GetEdges() < 1000 || G->GetNodes() < 300) {
    RAGM.MLEGradAscent(0.0001, 1000 * G->GetNodes(), "", StepAlpha, StepBeta);
  } else {
    RAGM.MLEGradAscentParallel(0.0001, 1000, NumThreads, "", StepAlpha, StepBeta);
  }
	RAGM.GetCmtyVV(EstCmtyVV);
  FILE* F = fopen(OutFPrx.CStr(), "wt"); 
  for (int c = 0; c < EstCmtyVV.Len(); c++) {
    for (int u = 0; u < EstCmtyVV[c].Len(); u++) {
      fprintf(F, "%d", EstCmtyVV[c][u].Val);
      if (u < EstCmtyVV[c].Len() - 1) { fprintf(F, "\t"); }
    }
    fprintf(F, "\n");
  }
  fclose(F);

	Catch

  printf("\nrun time: %s (%s)\n", ExeTm.GetTmStr(), TSecTm::GetCurTm().GetTmStr().CStr());

  return 0;
}
