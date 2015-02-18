// parser.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <imdb.h>

int main(int argc, char* argv[]) {
	Env = TEnv(argc, argv, TNotify::StdNotify);
	Env.PrepArgs(TStr::Fmt("realcom. build: %s, %s. Time: %s", __TIME__, __DATE__, TExeTm::GetCurTm()));
	TExeTm ExeTm;
	Try
  TStr GFPrx = Env.GetIfArgPrefixStr("-g:", "", "Graph data prefix");
	TStr InFPrx = Env.GetIfArgPrefixStr("-d:", "ppilc", "Community data set prefix");
	TStr SettingFNm = Env.GetIfArgPrefixStr("-s:", "Settings.ini", "Setting file name");
	TStr CmtySettingNm = Env.GetIfArgPrefixStr("-:", "FileNames.ini", "Setting file name for community affiliation files");
	TStrStrH SettingsH, CmtyFNmH;
	TYUtil::ReadStrStrH(SettingFNm, SettingsH, ssfTabSep);
	TYUtil::ReadStrStrH(CmtySettingNm, CmtyFNmH, ssfTabSep);
	IAssert(SettingsH.IsKey("DATAPATH"));
	IAssert(SettingsH.IsKey("OUTPATH"));

	TStr DataPath = SettingsH.GetDat("DATAPATH");
	TStr OutDir = SettingsH.GetDat("OUTPATH");

	printf("Setting loaded\n");
  if (GFPrx.Empty()) { GFPrx = InFPrx; }
	TFIn GFIn(DataPath + GFPrx + ".ungraph");
	TStr CFNm = DataPath + CmtyFNmH.GetDat(InFPrx);
	TFIn CFIn(CFNm);
	PUNGraph G = TUNGraph::Load(GFIn);
	THash<TInt, TIntV> CmtyVH;
	if (CFNm.GetFExt() == ".cmtyVV") {
		TVec<TIntV> CmtyVV(CFIn);
		CmtyVH.Gen(CmtyVV.Len());
		for (int c = 0; c < CmtyVV.Len(); c++) {
			CmtyVH.AddDat(c, CmtyVV[c]);
		}
	} else {
		CmtyVH.Load(CFIn);
	}
	printf("graph and communities are loaded\n");
  TStr StatFNm = OutDir + SettingsH.GetDat("STATFILE");
  TStrFltH StatH;
	TCommStruct CA(G, CmtyVH);
  CA.Init();
  CA.CutTooSmallCommunity(3);
  CA.LoadQ(InFPrx + ".QV");
  CA.LoadAvgDistForCommunity(InFPrx + ".AvgDistV");
  CA.LoadAvgDistForCommunityForRewired(InFPrx + ".AvgDistRewiredV");
  CA.PlotMemSize(InFPrx, OutDir);
  CA.PlotGroupStat(InFPrx, OutDir);
  CA.PlotNodeStat(InFPrx, OutDir);
	CA.PlotGroupSize(InFPrx, OutDir);
  CA.PlotMaxOlpFrac(InFPrx, OutDir);
  CA.GetMemSzStats(StatH);

  
  FILE* F = fopen(StatFNm.CStr(), "at");
  fprintf(F, "%s\t%d\t%d\t%d\t", InFPrx.CStr(), (int) StatH.GetDat("Nodes").Val, (int) StatH.GetDat("Edges").Val, (int) StatH.GetDat("Coms").Val);
  fprintf(F, "%.2f\t%.2f\t", StatH.GetDat("CCF").Val, StatH.GetDat("EffDiam").Val);
  fprintf(F, "%.1f\t%.1f\t%.3f\t", StatH.GetDat("MeanSz").Val, StatH.GetDat("MedianSz").Val, StatH.GetDat("FracHigh2MeanSz").Val);
  fprintf(F, "%.1f\t%.1f\t%.3f\t", StatH.GetDat("MeanMem").Val, StatH.GetDat("MedianMem").Val, StatH.GetDat("FracHigh2MeanMem").Val);
  fprintf(F, "\n");
  fclose(F);
  Catch

  printf("\nrun time: %s (%s)\n", ExeTm.GetTmStr(), TSecTm::GetCurTm().GetTmStr().CStr());

  return 0;
}


