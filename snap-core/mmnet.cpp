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
