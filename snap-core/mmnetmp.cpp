//
//  mmnetmp.cpp
//  snap-core
//
//  Created by Sameep Bagadia on 2/18/15.
//  Copyright (c) 2015 infolab. All rights reserved.
//


void GetPageRankMM(const PSVNetMP& Graph, TVec<TIntFltH>& PRankHV, const double& C, const double& Eps, const int& MaxIter) {
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
    for (TSVNetMP::TNodeI NI = Graph->BegNI(NType); NI < Graph->EndNI(NType); NI++) {
      PRankHV[NType].AddDat(NI.GetId(), 1.0/NNodes);
    }
  }
  TFltV TmpV(NNodes);
  for (int iter = 0; iter < MaxIter; iter++) {
    int j = 0;
    for (int NType = 0; NType < NTypeCnt; NType++) {
      for (TSVNetMP::TNodeI NI = Graph->BegNI(NType); NI < Graph->EndNI(NType); NI++, j++) {
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
/*
 void GetSimRankMM(const PSVNetMP& Graph, THash<TPair<TIntPr, TIntPr>, TFlt>& SRankH, const double& C, const double& Eps, const int& MaxIter) {
 const int NNodes = Graph->GetNodes();
 int NTypeCnt = Graph->GetNTypeCnt();
 int ETypeCnt = Graph->GetETypeCnt();
 for (int NType = 0; NType < NTypeCnt; NType++) {
 for (TSVNetMP::TNodeI NI = Graph->BegNI(NType); NI < Graph->EndNI(NType); NI++, j++) {
 SRankH.AddDat(TPair<TIntPr, TIntPr>(TPair<TInt, TInt>(NType, NI.GetId()), TPair<TInt, TInt>(NType, NI.GetId())), 1.0);
 }
 }
 THash<TPair<TIntPr, TIntPr>, TFlt>& SRankHOld;
 for (int iter = 0; iter < MaxIter; iter++) {
 //copy the entire hashmap and clear the original
 SRankHOld = SRankH;
 SRankH.Clr();
 
 //for each of the present values, add to the next set of values
 
 }
 }
 */

PSVNetMP TSVNetMP::GetSubGraph(TIntV NTypeV, TIntV ETypeV) {
  // if nodetype vector is empty, initialize it to proper node types
  
  TVec<bool> NTypeBool;
  TVec<bool> ETypeBool;
  
  int NTypeCnt = GetNTypeCnt();
  int ETypeCnt = GetETypeCnt();
  
  //initialize the vectors
  for (int NType = 0; NType < NTypeCnt; NType++) {
    NTypeBool.Add(false);
  }
  for (int EType = 0; EType < ETypeCnt; EType++) {
    ETypeBool.Add(false);
  }
  
  if (NTypeV.Empty()) {
    for (int i = 0; i < ETypeV.Len(); i++) {
      int EType = ETypeV[i];
      ETypeBool[EType] = true;
      NTypeBool[GetSrcNType(EType)] = true;
      NTypeBool[GetDstNType(EType)] = true;
    }
  }
  // if edgetype vector is empty, initialize it to proper edge types
  else if (ETypeV.Empty()) {
    for (int i = 0; i < NTypeV.Len(); i++) {
      NTypeBool[NTypeV[i]] = true;
    }
    for (int i = 0; i < ETypeCnt; i++) {
      if (NTypeBool[GetSrcNType(i)] && NTypeBool[GetDstNType(i)]) {
        ETypeBool[i] = true;
      }
    }
  }
  else {
    for (int i = 0; i < NTypeV.Len(); i++) {
      NTypeBool[NTypeV[i]] = true;
    }
    for (int i = 0; i < ETypeV.Len(); i++) {
      ETypeBool[ETypeV[i]] = true;
    }
  }
  /*
   // sort and remove duplicates in NTypeV and ETYpeV
   NTypeV.Merge();
   ETypeV.Merge();
   */
  
  PSVNetMP Graph = new TSVNetMP();
  //adding ntypes
  for (int NType = 0; NType < GetNTypeCnt(); NType++) {
    Graph->AddNType();
  }
  //adding etypes
  for (int EType = 0; EType < GetETypeCnt(); EType++) {
    Graph->AddEType(GetSrcNType(EType), GetDstNType(EType));
  }
  //adding nodes
  for (int NType = 0; NType < NTypeCnt ; NType++) {
    if (NTypeBool[NType]) {
      for (TSVNetMP::TNodeI NI = BegNI(NType); NI < EndNI(NType); NI++) {
        int NId = NI.GetId();
        Graph->AddNode(NType, NId);
      }
    }
  }
  //adding edges
  for (int EType = 0; EType < ETypeCnt; EType++) {
    if (ETypeBool[EType]) {
      for (THashMP<TInt, TEdge>::TIter it = EdgeHV[EType].BegI(); it < EdgeHV[EType].EndI(); it++) {
        int SrcNId = it.GetDat().GetSrcNId();
        int DstNId = it.GetDat().GetDstNId();
        int EId = it.GetDat().GetId();
        Graph->AddEdge(SrcNId, DstNId, EType, EId);
      }
    }
  }
  return Graph;
  
}

#ifdef _OPENMP
PSVNetMP TSVNetMP::GetSubGraphMP(TIntV NTypeV, TIntV ETypeV) {
  // if nodetype vector is empty, initialize it to proper node types
  
  float start = omp_get_wtime();
  
  TVec<bool> NTypeBool;
  TVec<bool> ETypeBool;
  
  int NTypeCnt = GetNTypeCnt();
  int ETypeCnt = GetETypeCnt();
  
  //initialize the vectors
  for (int NType = 0; NType < NTypeCnt; NType++) {
    NTypeBool.Add(false);
  }
  for (int EType = 0; EType < ETypeCnt; EType++) {
    ETypeBool.Add(false);
  }
  
  if (NTypeV.Empty()) {
    for (int i = 0; i < ETypeV.Len(); i++) {
      int EType = ETypeV[i];
      ETypeBool[EType] = true;
      NTypeBool[GetSrcNType(EType)] = true;
      NTypeBool[GetDstNType(EType)] = true;
    }
  }
  // if edgetype vector is empty, initialize it to proper edge types
  else if (ETypeV.Empty()) {
    for (int i = 0; i < NTypeV.Len(); i++) {
      NTypeBool[NTypeV[i]] = true;
    }
    for (int i = 0; i < ETypeCnt; i++) {
      if (NTypeBool[GetSrcNType(i)] && NTypeBool[GetDstNType(i)]) {
        ETypeBool[i] = true;
      }
    }
  }
  else {
    for (int i = 0; i < NTypeV.Len(); i++) {
      NTypeBool[NTypeV[i]] = true;
    }
    for (int i = 0; i < ETypeV.Len(); i++) {
      ETypeBool[ETypeV[i]] = true;
    }
  }
  /*
   // sort and remove duplicates in NTypeV and ETYpeV
   NTypeV.Merge();
   ETypeV.Merge();
   */
  
  PSVNetMP Graph = new TSVNetMP();
  //adding ntypes
  for (int NType = 0; NType < GetNTypeCnt(); NType++) {
    Graph->AddNType();
    if (NTypeBool[NType]) {
      Graph->ReserveNodes(NType, MxNIdV[NType]);
      Graph->SetMxNId(NType, MxNIdV[NType]);
    }
  }
  //adding etypes
  for (int EType = 0; EType < GetETypeCnt(); EType++) {
    Graph->AddEType(GetSrcNType(EType), GetDstNType(EType));
    if (ETypeBool[EType]) {
      Graph->ReserveEdges(EType, MxEIdV[EType]);
      Graph->SetMxEId(EType, MxEIdV[EType]);
    }
  }
  
  TimeV[0] += omp_get_wtime() - start;
  start = omp_get_wtime();
  
  //adding nodes
  TIntV NTypeVec;
  //TIntV NIdVec;
  TIntV NodeKeyIdVec;
  TIntV KeyIdVec;
  //#pragma omp parallel for schedule(dynamic, 1)
  for (int NType = 0; NType < NTypeCnt; NType++) {
    if (NTypeBool[NType]) {
      int MaxKeys = NodeHV[NType].GetMxKeyIds();
      for (int KeyId = NodeHV[NType].FFirstKeyId(); KeyId < MaxKeys; KeyId++) {
        if (NodeHV[NType].IsKeyId2(KeyId)) {
          TInt NId = NodeHV[NType].GetKey(KeyId);
          TInt NodeKeyId = Graph->AddNodeToHash(NId, NType, ETypeCnt);
          NTypeVec.Add(NType);
          //NIdVec.Add(NId);
          NodeKeyIdVec.Add(NodeKeyId);
          KeyIdVec.Add(KeyId);
        }
      }
      Graph->SetNCnt(NType, GetNodes(NType));
    }
  }
  TimeV[1] += omp_get_wtime() - start;
  start = omp_get_wtime();
  //adding edges to nodes
  TInt NodeCnt = NTypeVec.Len();
  //TIntV test(NodeCnt);
  //#pragma omp parallel for schedule(dynamic,10000)
  for (int i = 0; i < NodeCnt; i++) {
    TInt NType = NTypeVec[i];
    TInt KeyId = KeyIdVec[i];
    TVec<TIntV> InEIdVV = NodeHV[NType][KeyId].GetInEIdVV();
    TVec<TIntV> OutEIdVV = NodeHV[NType][KeyId].GetOutEIdVV();
    //test[i] = InEIdVV.Len() + OutEIdVV.Len();
    //printf("%d %d\n", InEIdVV.Len(), OutEIdVV.Len());
    Graph->AddEdgesToNode(NType, NodeKeyIdVec[i], InEIdVV, OutEIdVV);
  }
  
  /*
  for (int NType = 0; NType < NTypeCnt ; NType++) {
    if (NTypeBool[NType]) {
      int MaxKeys = NodeHV[NType].GetMxKeyIds();
      #pragma omp parallel for schedule(dynamic,10000)
      for (int KeyId = NodeHV[NType].FFirstKeyId(); KeyId < MaxKeys; KeyId++) {
        if (NodeHV[NType].IsKeyId2(KeyId)) {
          TInt NId = NodeHV[NType].GetKey(KeyId);
          TVec<TIntV> InEIdVV = NodeHV[NType][KeyId].GetInEIdVV();
          TVec<TIntV> OutEIdVV = NodeHV[NType][KeyId].GetOutEIdVV();
          Graph->AddEdgesToNode(NId, NType, InEIdVV, OutEIdVV);
        }
      }
      Graph->SetNCnt(NType, GetNodes(NType));
    }
  }*/
  TimeV[2] += omp_get_wtime() - start;
  start = omp_get_wtime();
  //adding edges
  for (int EType = 0; EType < ETypeCnt; EType++) {
    if (ETypeBool[EType]) {
      int MaxKeys = EdgeHV[EType].GetMxKeyIds();
      //#pragma omp parallel for schedule(dynamic,10000)
      for (int KeyId = EdgeHV[EType].FFirstKeyId(); KeyId < MaxKeys; KeyId++) {
        if (EdgeHV[EType].IsKeyId2(KeyId)) {
          TInt EId = EdgeHV[EType].GetKey(KeyId);
          TInt SrcNId = EdgeHV[EType][KeyId].GetSrcNId();
          TInt DstNId = EdgeHV[EType][KeyId].GetDstNId();
          Graph->AddEdgeToHash(EId, SrcNId, DstNId, EType);
        }
      }
      Graph->SetECnt(EType, GetEdges(EType));
    }
  }
  TimeV[3] += omp_get_wtime() - start;
  
  //int test_sum = 0;
  //for (int i = 0; i < NodeCnt; i++) {test_sum += test[i];}
  //printf("%d\n", test_sum);
  
  return Graph;
  
}
#endif


PNEANet TSVNetMP::GetSubGraphTNEANet(TIntV NTypeV, TIntV ETypeV, TIntIntH& Offsets) {
  // if nodetype vector is empty, initialize it to proper node types
  
  TVec<bool> NTypeBool;
  TVec<bool> ETypeBool;
  
  int NTypeCnt = GetNTypeCnt();
  int ETypeCnt = GetETypeCnt();
  
  //initialize the vectors
  for (int NType = 0; NType < NTypeCnt; NType++) {
    NTypeBool.Add(false);
  }
  for (int EType = 0; EType < ETypeCnt; EType++) {
    ETypeBool.Add(false);
  }
  
  if (NTypeV.Empty()) {
    for (int i = 0; i < ETypeV.Len(); i++) {
      int EType = ETypeV[i];
      ETypeBool[EType] = true;
      NTypeBool[GetSrcNType(EType)] = true;
      NTypeBool[GetDstNType(EType)] = true;
    }
  }
  // if edgetype vector is empty, initialize it to proper edge types
  else if (ETypeV.Empty()) {
    for (int i = 0; i < NTypeV.Len(); i++) {
      NTypeBool[NTypeV[i]] = true;
    }
    for (int i = 0; i < ETypeCnt; i++) {
      if (NTypeBool[GetSrcNType(i)] && NTypeBool[GetDstNType(i)]) {
        ETypeBool[i] = true;
      }
    }
  }
  else {
    for (int i = 0; i < NTypeV.Len(); i++) {
      NTypeBool[NTypeV[i]] = true;
    }
    for (int i = 0; i < ETypeV.Len(); i++) {
      ETypeBool[ETypeV[i]] = true;
    }
  }
  
  // sort and remove duplicates in NTypeV and ETYpeV
  //NTypeV.Merge();
  //ETypeV.Merge();
  
  
  PNEANet Graph = new TNEANet();
  //adding nodes
  TInt offset = 0;
  TIntV OffsetV;
  for (int NType = 0; NType < GetNTypeCnt(); NType++) {
    if(NTypeBool[NType]){
      Offsets.AddDat(NType, offset);
      OffsetV.Add(offset);
      for (TSVNetMP::TNodeI NI = BegNI(NType); NI < EndNI(NType); NI++) {
        int NId = NI.GetId();
        Graph->AddNode(NId + offset);
      }
      offset += MxNIdV[NType];
    }
    else {
      OffsetV.Add(-1);
    }
  }
  
  //adding edges
  for (int EType = 0; EType < ETypeCnt; EType++) {
    if (ETypeBool[EType]) {
      int SrcNType = GetSrcNType(EType);
      int DstNType = GetDstNType(EType);
      int SrcOffset = OffsetV[SrcNType];
      int DstOffset = OffsetV[DstNType];
      for (THashMP<TInt, TEdge>::TIter it = EdgeHV[EType].BegI(); it < EdgeHV[EType].EndI(); it++) {
        int SrcNId = it.GetDat().GetSrcNId();
        int DstNId = it.GetDat().GetDstNId();
        //int EId = it.GetDat().GetId();
        int NewSrcNId = SrcNId + SrcOffset;
        int NewDstNId = DstNId + DstOffset;
        Graph->AddEdge(NewSrcNId, NewDstNId);
      }
    }
  }
  return Graph;
}


#ifdef _OPENMP
// Page Rank -- there are two different implementations (uncomment the desired 2 lines):
//   Berkhin -- (the correct way) see Algorithm 1 of P. Berkhin, A Survey on PageRank Computing, Internet Mathematics, 2005
//   iGraph -- iGraph implementation(which treats leaked PageRank in a funny way)
// This is a parallel, optimized version.

void GetPageRankMMMP2(const PSVNetMP& Graph, TVec<TIntFltH>& PRankHV, const double& C, const double& Eps, const int& MaxIter) {
  const int NNodes = Graph->GetNodes();
  int NTypeCnt = Graph->GetNTypeCnt();
  int ETypeCnt = Graph->GetETypeCnt();
  //TIntV NV;
  TVec<TSVNetMP::TNodeI> NV;
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
    for (TSVNetMP::TNodeI NI = Graph->BegNI(NType); NI < Graph->EndNI(NType); NI++) {
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
    TSVNetMP::TNodeI NI = NV[j];
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
      TSVNetMP::TNodeI NI = NV[j];
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
      TSVNetMP::TNodeI NI = NV[i];
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
    TSVNetMP::TNodeI NI = NV[i];
    int NId = NI.GetId();
    int NType = NI.GetNType();
    //PRankH.AddDat(NI.GetId(), PRankV[NI.GetId()]);
    //PRankHV[NType][i] = PRankVV[NType][NId];
    PRankHV[NType].AddDat(NId, PRankVV[NType][NId]);
  }
}

#endif

void GetBfsLevelMM(const PSVNetMP& Graph, TVec<TIntIntH>& BfsLevelHV, const int& StartNId, const int& StartNType) {
  //printf("BFS called\n");
  int NTypeCnt = Graph->GetNTypeCnt();
  int ETypeCnt = Graph->GetETypeCnt();
  
  BfsLevelHV = TVec<TIntIntH>();
  
  
  for (int i = 0; i < NTypeCnt; i++) {
    TIntIntH BfsLevelH;
    BfsLevelH.Gen(Graph->GetNodes(i));
    BfsLevelHV.Add(BfsLevelH);
  }
  for (int NType = 0; NType < NTypeCnt; NType++) {
    for (TSVNetMP::TNodeI NI = Graph->BegNI(NType); NI < Graph->EndNI(NType); NI++) {
      BfsLevelHV[NType].AddDat(NI.GetId(), -1);
    }
  }
  //printf("reached here1\n");
  int LevelCnt = 0;
  TVec<TIntPr> CurrLevel;
  TVec<TIntPr> NextLevel;
  TIntPr StartNode(StartNType, StartNId);
  CurrLevel.Add(StartNode);
  BfsLevelHV[StartNode.Val1].AddDat(StartNode.Val2, 0);
  //printf("reached here2\n");
  while (CurrLevel.Len() > 0) {
    //printf("bfs level: %d\n", LevelCnt);
    LevelCnt++;
    //printf("level count: %d\n", LevelCnt);
    NextLevel.Clr();
    for (int i = 0; i < CurrLevel.Len(); i++) {
      int NId = CurrLevel[i].Val2;
      int NType = CurrLevel[i].Val1;
      //printf("NId = %d, NType = %d\n", NId, NType);
      TSVNetMP::TNodeI NI = Graph->GetNI(NId, NType);
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
void GetBfsLevelMMMP(const PSVNetMP& Graph, TVec<TIntIntH>& BfsLevelHV, const int& StartNId, const int& StartNType) {
  int NTypeCnt = Graph->GetNTypeCnt();
  int ETypeCnt = Graph->GetETypeCnt();
  const int NNodes = Graph->GetNodes();
  
  BfsLevelHV = TVec<TIntIntH>();
  
  
  for (int i = 0; i < NTypeCnt; i++) {
    TIntIntH BfsLevelH;
    BfsLevelH.Gen(Graph->GetNodes(i));
    BfsLevelHV.Add(BfsLevelH);
  }
  for (int NType = 0; NType < NTypeCnt; NType++) {
    for (TSVNetMP::TNodeI NI = Graph->BegNI(NType); NI < Graph->EndNI(NType); NI++) {
      BfsLevelHV[NType].AddDat(NI.GetId(), -1);
    }
  }
  int LevelCnt = 0;
  TVec<TIntPr> CurrLevel;
  TVec<TIntPr> NextLevel(NNodes);
  TIntPr StartNode(StartNType, StartNId);
  CurrLevel.Add(StartNode);
  BfsLevelHV[StartNode.Val1].AddDat(StartNode.Val2, 0);
  while (CurrLevel.Len() > 0) {
    //printf("level count: %d\n", LevelCnt);
    TVec<TIntPr> NextLevel(NNodes);
    NextLevel.Reduce(0);
#pragma omp parallel for schedule(dynamic,1000)
    for (int i = 0; i < CurrLevel.Len(); i++) {
      int NId = CurrLevel[i].Val2;
      int NType = CurrLevel[i].Val1;
      //printf("NId = %d, NType = %d\n", NId, NType);
      TSVNetMP::TNodeI NI = Graph->GetNI(NId, NType);
      for (int EType = 0; EType < ETypeCnt; EType++) {
        //printf("EType = %d\n", EType);
        for (int e = 0; e < NI.GetOutDeg(EType); e++) {
          const int OutEId = NI.GetOutEId(e, EType);
          const int OutNType = Graph->GetDstNType(EType);
          const int OutNId = Graph->GetDstNId(OutEId, EType);
          //printf("OutEId = %d, OutNType = %d, OutNId = %d\n", OutEId, OutNType, OutNId);
          //TIntV test(100);
          if (BfsLevelHV[OutNType].GetDat(OutNId) == -1) {
            TIntPr OutNode(OutNType, OutNId);
            printf("Reached before MP2\n");
            //test.AddAtm(OutNId);
            NextLevel.AddAtm(OutNode);
            printf("Reached after MP2\n");
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

int GetBfsLevelMMMP2(const PSVNetMP& Graph, TVec<TIntV >& BfsLevelVV, const int& StartNId, const int& StartNType) {
  int NTypeCnt = Graph->GetNTypeCnt();
  int ETypeCnt = Graph->GetETypeCnt();
  const int NNodes = Graph->GetNodes();
  //int MxNId = Graph->GetMxNId();
  int NonNodeDepth = 2147483647; // INT_MAX
  int InfDepth = 2147483646; // INT_MAX - 1
  //ShortestDists.Gen(MxNId);
  BfsLevelVV.Gen(NTypeCnt);
  for (int NType = 0; NType < NTypeCnt; NType++) {
    BfsLevelVV[NType].Gen(Graph->GetMxNId(NType));
  }
  for (int NType = 0; NType < NTypeCnt; NType++) {
    int MxNId = Graph->GetMxNId(NType);
#pragma omp parallel for schedule(dynamic,10000)
    for (int NId = 0; NId < MxNId; NId++) {
      if (Graph->IsNode(NId, NType)) { BfsLevelVV[NType][NId] = InfDepth; }
      else { BfsLevelVV[NType][NId] = NonNodeDepth; }
    }
  }
  
  /*
   TIntV NIdVec1(NNodes, 0); // ensure enough capacity
   TIntV NIdVec2(NNodes, 0); // ensure enough capacity
   TIntV NTypeVec1(NNodes, 0); // ensure enough capacity
   TIntV NTypeVec2(NNodes, 0); // ensure enough capacity
   */
  TIntPrV Vec1(NNodes);
  TIntPrV Vec2(NNodes); // ensure enough capacity
  
  BfsLevelVV[StartNType][StartNId] = 0;
  /*
   TIntV* PCurNIdV = &NIdVec1;
   TIntV* PCurNTypeV = &NTypeVec1;
   PCurNIdV->Add(StartNId);
   PCurNTypeV->Add(StartNType);
   TIntV* PNextNIdV = &NIdVec2;
   TIntV* PNextNTypeV = &NTypeVec2;
   */
  TIntPrV* PCurV = &Vec1;
  TIntPr StartNode(StartNId, StartNType);
  PCurV->Add(StartNode);
  TIntPrV* PNextV = &Vec2;
  int Depth = 0; // current depth
  
  while (!PCurV->Empty()) {
    Depth++; // increase depth
    //printf("Starting depth : %d\n", Depth);
#pragma omp parallel for schedule(dynamic,10000)
    for (int i = 0; i < PCurV->Len(); i++) {
      int NId = PCurV->GetVal(i).Val1;
      int NType = PCurV->GetVal(i).Val2;
      TSVNetMP::TNodeI NI = Graph->GetNI(NId, NType);
      for (int EType = 0; EType < ETypeCnt; EType++) {
        for (int e = 0; e < NI.GetOutDeg(EType); e++) {
          
          const int EId = NI.GetOutEId(e, EType);
          const int OutNId = Graph->GetDstNId(EId, EType);
          const int OutNType = Graph->GetDstNType(EType);
          //printf("OutNid : %d, OutNType: %d\n", OutNId, OutNType);
          //const int OutNId = NI.GetOutNId(e);
          if (__sync_bool_compare_and_swap(&(BfsLevelVV[OutNType][OutNId].Val), InfDepth, Depth)) {
            //printf("before AddAtm\n");
            TIntPr out_node(OutNId, OutNType);
            PNextV->AddAtm(out_node);
            //PNextNTypeV->AddAtm(OutNType);
            //printf("After AddAtm\n");
          }
        }
      }
      
    }
    //      #pragma omp parallel for schedule(dynamic,10000)
    //      for (int NId = 0; NId < MxNId; NId++) {
    //        if (ShortestDists[NId] == InfDepth) {
    //          typename PGraph::TObj::TNodeI NI = Graph->GetNI(NId);
    //          for (int e = 0; e < NI.GetInDeg(); e++) {
    //            const int InNId = NI.GetInNId(e);
    //            if (ShortestDists[InNId] < Depth) {
    //              ShortestDists[NId] = Depth;
    //              PNextV->AddAtm(NId);
    //              break;
    //            }
    //          }
    //        }
    //      }
    // swap pointer, no copying
    TIntPrV* Tmp = PCurV;
    PCurV = PNextV;
    PNextV = Tmp;
    
    //TIntV* TmpNType = PCurNTypeV;
    //PCurNTypeV = PNextNTypeV;
    //PNextNTypeV = TmpNType;
    // clear next
    PNextV->Reduce(0); // reduce length, does not initialize new array
    //PNextNTypeV->Reduce(0);
  }
  return Depth-1;
}


#endif
