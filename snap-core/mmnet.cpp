//
//  mmnet.cpp
//  snap-core
//
//  Created by Sameep Bagadia on 11/18/14.
//  Copyright (c) 2014 infolab. All rights reserved.
//


//typedef TSVNet* PSVNet;
void GetPageRankMM(const PSVNet& Graph, TVec<TIntFltH>& PRankHV, const double& C, const double& Eps, const int& MaxIter) {
  const int NNodes = Graph->GetNodes();
  int NTypeCnt = Graph->GetNTypeCnt();
  int ETypeCnt = Graph->GetETypeCnt();
  //const double OneOver = 1.0/double(NNodes);
  for (int i = 0; i < NTypeCnt; i++) {
    TIntFltH PRankH;
    PRankH.Gen(Graph->GetNodes(i));
    PRankHV.Add(PRankH);
  }
  for (int NType = 0; NType < NTypeCnt; NType++) {
    for (TSVNet::TNodeI NI = Graph->BegNI(NType); NI < Graph->EndNI(NType); NI++) {
      PRankHV[NType].AddDat(NI.GetId(), 1.0/NNodes);
    }
  }
  TFltV TmpV(NNodes);
  for (int iter = 0; iter < MaxIter; iter++) {
    int j = 0;
    for (int NType = 0; NType < NTypeCnt; NType++) {
      for (TSVNet::TNodeI NI = Graph->BegNI(NType); NI < Graph->EndNI(NType); NI++, j++) {
        TmpV[j] = 0;
        for (int EType = 0; EType < ETypeCnt; EType++) {
          for (int e = 0; e < NI.GetInDeg(EType); e++) {
            const int InEId = NI.GetInEId(e, EType);
            const int InNType = Graph->GetSrcNType(EType);
            const int InNId = Graph->GetSrcNId(InEId, EType);
            const int OutDeg = Graph->GetNI(InNId, InNType).GetOutDeg();
            TmpV[j] += PRankHV[InNType].GetDat(InNId) / OutDeg;
          }
        }
        TmpV[j] =  C*TmpV[j]; // Berkhin (the correct way of doing it)
        //TmpV[j] =  C*TmpV[j] + (1.0-C)*OneOver; // iGraph
      }
    }
    double diff=0, sum=0, NewVal;
    for (int i = 0; i < TmpV.Len(); i++) { sum += TmpV[i]; }
    const double Leaked = (1.0-sum) / double(NNodes);
    j = 0;
    for (int NType = 0; NType < NTypeCnt; NType++) {
      for (int i = 0; i < PRankHV[NType].Len(); i++, j++) {
        NewVal = TmpV[j] + Leaked; // Berkhin
        //NewVal = TmpV[j] / sum;  // iGraph
        diff += fabs(NewVal-PRankHV[NType][i]);
        PRankHV[NType][i] = NewVal;
      }
    }
    if (diff < Eps) { break; }
  }
}





#ifdef _OPENMP
// Page Rank -- there are two different implementations (uncomment the desired 2 lines):
//   Berkhin -- (the correct way) see Algorithm 1 of P. Berkhin, A Survey on PageRank Computing, Internet Mathematics, 2005
//   iGraph -- iGraph implementation(which treats leaked PageRank in a funny way)
// This is a parallel, optimized version.

void GetPageRankMMMP2(const PSVNet& Graph, TVec<TIntFltH>& PRankHV, const double& C, const double& Eps, const int& MaxIter) {
  const int NNodes = Graph->GetNodes();
  int NTypeCnt = Graph->GetNTypeCnt();
  int ETypeCnt = Graph->GetETypeCnt();
  //TIntV NV;
  TVec<TSVNet::TNodeI> NV;
  PRankHV = TVec<TIntFltH>();
  //const double OneOver = 1.0/double(NNodes);
  for (int i = 0; i < NTypeCnt; i++) {
    TIntFltH PRankH;
    PRankH.Gen(Graph->GetNodes(i));
    PRankHV.Add(PRankH);
  }
  //PRankH.Gen(NNodes);
  //int MxId = -1;
  //time_t t = time(0);
  //printf("%s", ctime(&t));
  for (int NType = 0; NType < NTypeCnt; NType++) {
    for (TSVNet::TNodeI NI = Graph->BegNI(NType); NI < Graph->EndNI(NType); NI++) {
      NV.Add(NI);
      PRankHV[NType].AddDat(NI.GetId(), 1.0/NNodes);
    }
  }
  /*
  for (typename PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
    NV.Add(NI);
    PRankH.AddDat(NI.GetId(), 1.0/NNodes);
    int Id = NI.GetId();
    if (Id > MxId) {
      MxId = Id;
    }
    //IAssert(NI.GetId() == PRankH.GetKey(PRankH.Len()-1));
  }*/
  //t = time(0);
  //printf("%s", ctime(&t));
  TVec<TFltV> PRankVV;
  TVec<TIntV> OutDegVV;
  for (int NType = 0; NType < NTypeCnt; NType++) {
    TFltV PRankV(Graph->GetMxNId(NType)+1);
    TIntV OutDegV(Graph->GetMxNId(NType)+1);
    PRankVV.Add(PRankV);
    OutDegVV.Add(OutDegV);
  }
  /*
  TFltV PRankV(MxId+1);
  TIntV OutDegV(MxId+1);
  */
  //for (typename PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
#pragma omp parallel for schedule(dynamic,10000)
  for (int j = 0; j < NNodes; j++) {
    TSVNet::TNodeI NI = NV[j];
    int NId = NI.GetId();
    int NType = NI.GetNType();
    PRankVV[NType][NId] = 1.0/NNodes;
    OutDegVV[NType][NId] = NI.GetOutDeg();
    
    /*
    typename PGraph::TObj::TNodeI NI = NV[j];
    int Id = NI.GetId();
    PRankV[Id] = 1.0/NNodes;
    OutDegV[Id] = NI.GetOutDeg();*/
  }
  
  TFltV TmpV(NNodes);
  
  //int hcount1 = 0;
  //int hcount2 = 0;
  for (int iter = 0; iter < MaxIter; iter++) {
    //time_t t = time(0);
    //printf("%s%d\n", ctime(&t),iter);
    //int j = 0;
    //for (typename PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++, j++) {
#pragma omp parallel for schedule(dynamic,10000)
    for (int j = 0; j < NNodes; j++) {
      TSVNet::TNodeI NI = NV[j];
      TFlt Tmp = 0;
      for (int EType = 0; EType < ETypeCnt; EType++) {
        for (int e = 0; e < NI.GetInDeg(EType); e++) {
          const int InEId = NI.GetInEId(e, EType);
          const int InNType = Graph->GetSrcNType(EType);
          const int InNId = Graph->GetSrcNId(InEId, EType);
          const int OutDeg = OutDegVV[InNType][InNId];
          Tmp += PRankVV[InNType][InNId] / OutDeg;
        }
      }
      TmpV[j] = C*Tmp;
    
      
      
      /*
      //typename PGraph::TObj::TNodeI NI = Graph->GetNI(NV[j]);
      typename PGraph::TObj::TNodeI NI = NV[j];
      TFlt Tmp = 0;
      //TmpV[j] = 0;
      for (int e = 0; e < NI.GetInDeg(); e++) {
        const int InNId = NI.GetInNId(e);
        //hcount1++;
        //const int OutDeg = Graph->GetNI(InNId).GetOutDeg();
        const int OutDeg = OutDegV[InNId];
        //if (OutDeg != OutDegV[InNId]) {
        //printf("*** ERROR *** InNId %d, OutDeg %d, OutDegV %d\n",
        //InNId, OutDeg, OutDegV[InNId].Val);
        //}
        if (OutDeg > 0) {
          //hcount2++;
          //TmpV[j] += PRankH.GetDat(InNId) / OutDeg;
          //TmpV[j] += PRankV[InNId] / OutDeg;
          Tmp += PRankV[InNId] / OutDeg;
        }
      }
      TmpV[j] =  C*Tmp; // Berkhin (the correct way of doing it)
      //TmpV[j] =  C*TmpV[j]; // Berkhin (the correct way of doing it)
      ////TmpV[j] =  C*TmpV[j] + (1.0-C)*OneOver; // iGraph*/
    }
    double sum = 0;
#pragma omp parallel for reduction(+:sum) schedule(dynamic,10000)
    for (int i = 0; i < TmpV.Len(); i++) { sum += TmpV[i]; }
    const double Leaked = (1.0-sum) / double(NNodes);
    
    double diff = 0;
    //#pragma omp parallel for reduction(+:diff) schedule(dynamic,10000)
    //for (int i = 0; i < PRankH.Len(); i++) { // re-instert leaked PageRank
    //double NewVal = TmpV[i] + Leaked; // Berkhin
    ////NewVal = TmpV[i] / sum;  // iGraph
    //diff += fabs(NewVal-PRankH[i]);
    //PRankH[i] = NewVal;
    //}
#pragma omp parallel for reduction(+:diff) schedule(dynamic,10000)
    for (int i = 0; i < NNodes; i++) {
      TSVNet::TNodeI NI = NV[i];
      double NewVal = TmpV[i] + Leaked; // Berkhin
      //NewVal = TmpV[i] / sum;  // iGraph
      int NId = NI.GetId();
      int NType = NI.GetNType();
      diff += fabs(NewVal-PRankVV[NType][NId]);
      PRankVV[NType][NId] = NewVal;
    }
    //printf("counts %d %d\n", hcount1, hcount2);
    if (diff < Eps) { break; }
  }
  
#pragma omp parallel for schedule(dynamic,10000)
  for (int i = 0; i < NNodes; i++) {
    TSVNet::TNodeI NI = NV[i];
    int NId = NI.GetId();
    int NType = NI.GetNType();
    //PRankH.AddDat(NI.GetId(), PRankV[NI.GetId()]);
    //PRankHV[NType][i] = PRankVV[NType][NId];
    PRankHV[NType].AddDat(NId, PRankVV[NType][NId]);
  }
}

#endif

void GetBfsLevelMM(const PSVNet& Graph, TVec<TIntIntH>& BfsLevelHV, const int& StartNId, const int& StartNType) {
  int NTypeCnt = Graph->GetNTypeCnt();
  int ETypeCnt = Graph->GetETypeCnt();
  
  BfsLevelHV = TVec<TIntIntH>();
  
  
  for (int i = 0; i < NTypeCnt; i++) {
    TIntIntH BfsLevelH;
    BfsLevelH.Gen(Graph->GetNodes(i));
    BfsLevelHV.Add(BfsLevelH);
  }
  for (int NType = 0; NType < NTypeCnt; NType++) {
    for (TSVNet::TNodeI NI = Graph->BegNI(NType); NI < Graph->EndNI(NType); NI++) {
      BfsLevelHV[NType].AddDat(NI.GetId(), -1);
    }
  }
  int LevelCnt = 0;
  TVec<TIntPr> CurrLevel;
  TVec<TIntPr> NextLevel;
  TIntPr StartNode(StartNType, StartNId);
  CurrLevel.Add(StartNode);
  BfsLevelHV[StartNode.Val1].AddDat(StartNode.Val2, 0);
  while (CurrLevel.Len() > 0) {
    LevelCnt++;
    //printf("level count: %d\n", LevelCnt);
    NextLevel.Clr();
    for (int i = 0; i < CurrLevel.Len(); i++) {
      int NId = CurrLevel[i].Val2;
      int NType = CurrLevel[i].Val1;
      //printf("NId = %d, NType = %d\n", NId, NType);
      TSVNet::TNodeI NI = Graph->GetNI(NId, NType);
      for (int EType = 0; EType < ETypeCnt; EType++) {
        //printf("EType = %d\n", EType);
        for (int e = 0; e < NI.GetOutDeg(EType); e++) {
          const int OutEId = NI.GetOutEId(e, EType);
          const int OutNType = Graph->GetDstNType(EType);
          const int OutNId = Graph->GetDstNId(OutEId, EType);
          //printf("OutEId = %d, OutNType = %d, OutNId = %d\n", OutEId, OutNType, OutNId);
          if (BfsLevelHV[OutNType].GetDat(OutNId) == -1) {
            BfsLevelHV[OutNType].AddDat(OutNId, LevelCnt);
            TIntPr OutNode(OutNType, OutNId);
            NextLevel.Add(OutNode);
          }
        }
      }
    }
    CurrLevel = NextLevel;
  }
}

#ifdef _OPENMP
void GetBfsLevelMMMP(const PSVNet& Graph, TVec<TIntIntH>& BfsLevelHV, const int& StartNId, const int& StartNType) {
  int NTypeCnt = Graph->GetNTypeCnt();
  int ETypeCnt = Graph->GetETypeCnt();
  
  BfsLevelHV = TVec<TIntIntH>();
  
  
  for (int i = 0; i < NTypeCnt; i++) {
    TIntIntH BfsLevelH;
    BfsLevelH.Gen(Graph->GetNodes(i));
    BfsLevelHV.Add(BfsLevelH);
  }
  for (int NType = 0; NType < NTypeCnt; NType++) {
    for (TSVNet::TNodeI NI = Graph->BegNI(NType); NI < Graph->EndNI(NType); NI++) {
      BfsLevelHV[NType].AddDat(NI.GetId(), -1);
    }
  }
  int LevelCnt = 0;
  TVec<TIntPr> CurrLevel;
  TVec<TIntPr> NextLevel;
  TIntPr StartNode(StartNType, StartNId);
  CurrLevel.Add(StartNode);
  BfsLevelHV[StartNode.Val1].AddDat(StartNode.Val2, 0);
  while (CurrLevel.Len() > 0) {
    //printf("level count: %d\n", LevelCnt);
    NextLevel.Clr();
#pragma omp parallel for schedule(dynamic,1000)
    for (int i = 0; i < CurrLevel.Len(); i++) {
      int NId = CurrLevel[i].Val2;
      int NType = CurrLevel[i].Val1;
      //printf("NId = %d, NType = %d\n", NId, NType);
      TSVNet::TNodeI NI = Graph->GetNI(NId, NType);
      for (int EType = 0; EType < ETypeCnt; EType++) {
        //printf("EType = %d\n", EType);
        for (int e = 0; e < NI.GetOutDeg(EType); e++) {
          const int OutEId = NI.GetOutEId(e, EType);
          const int OutNType = Graph->GetDstNType(EType);
          const int OutNId = Graph->GetDstNId(OutEId, EType);
          //printf("OutEId = %d, OutNType = %d, OutNId = %d\n", OutEId, OutNType, OutNId);
          if (BfsLevelHV[OutNType].GetDat(OutNId) == -1) {
            TIntPr OutNode(OutNType, OutNId);
            NextLevel.AddMP(OutNode);
          }
        }
      }
    }
    LevelCnt++;
    NextLevel.Merge();
    for (int i = 0; i < NextLevel.Len(); i++) {
      BfsLevelHV[NextLevel[i].Val1].AddDat(NextLevel[i].Val2, LevelCnt);
    }
    CurrLevel = NextLevel;
  }
}


#endif
