//
//  mmnetmp.h
//  snap-core
//
//  Created by Sameep Bagadia on 2/18/15.
//  Copyright (c) 2015 infolab. All rights reserved.
//

#ifndef snap_core_mmnetmp_h
#define snap_core_mmnetmp_h

class TSVNetMP {
public:
  class TNode {
  private:
    TInt Id, Type;
    TVec<TIntV > InEIdVV, OutEIdVV;
  public:
    TNode() : Id(-1), Type(-1), InEIdVV(), OutEIdVV() { }
    TNode(const int& NId, const int& NType, const int& ETypeCnt) : Id(NId), Type(NType), InEIdVV(), OutEIdVV() {
      for (int i = 0; i < ETypeCnt; i++) {
        InEIdVV.Add(TIntV());
        OutEIdVV.Add(TIntV());
      }
    }
    TNode(const TNode& Node) : Id(Node.Id), Type(Node.Type), InEIdVV(Node.InEIdVV), OutEIdVV(Node.OutEIdVV) { }
    //  TNode(TSIn& SIn) : Id(SIn), InEIdV(SIn), OutEIdV(SIn) { }
    //  void Save(TSOut& SOut) const { Id.Save(SOut); InEIdV.Save(SOut); OutEIdV.Save(SOut); }
    int GetId() const { return Id; }
    int GetNType() const { return Type; }
    int GetDeg() const { return GetInDeg() + GetOutDeg(); }
    int GetInDeg(int EType) {return InEIdVV[EType].Len();}
    int GetInDeg() const {
      TInt InDeg = 0;
      for (int i = 0; i < InEIdVV.Len(); i++) {
        InDeg += InEIdVV[i].Len();
      }
      return InDeg;
    }
    int GetOutDeg(int EType) {return OutEIdVV[EType].Len();}
    int GetOutDeg() const {
      TInt OutDeg = 0;
      for (int i = 0; i < OutEIdVV.Len(); i++) {
        OutDeg += OutEIdVV[i].Len();
      }
      return OutDeg;
    }
    void AddEType() {
      InEIdVV.Add(TIntV());
      OutEIdVV.Add(TIntV());
    }
    void AddOutEdge(int EId, int EType) {
      OutEIdVV[EType].Add(EId);
    }
    void AddInEdge(int EId, int EType) {
      InEIdVV[EType].Add(EId);
    }
    int GetInEId(const int& EdgeN, const int& EType) const { return InEIdVV[EType][EdgeN]; }
    int GetOutEId(const int& EdgeN, const int& EType) const { return OutEIdVV[EType][EdgeN]; }
    
    void SetInEIdVV(TVec<TIntV>& InVV) {
      for (int i = 0; i < InVV.Len(); i++) {
        if (InVV[i].Len() > 0) {
          InEIdVV[i].MoveFrom(InVV[i]);
        }
      }
    }
    
    void SetOutEIdVV(TVec<TIntV>& OutVV) {
      for (int i = 0; i < OutVV.Len(); i++) {
        if (OutVV[i].Len() > 0) {
          OutEIdVV[i].MoveFrom(OutVV[i]);
        }
      }
    }
    
    TVec<TIntV> GetInEIdVV() { return InEIdVV; }
    TVec<TIntV> GetOutEIdVV() { return OutEIdVV; }
    
    // int GetNbrEId(const int& EdgeN) const { return EdgeN<GetOutDeg()?GetOutEId(EdgeN):GetInEId(EdgeN-GetOutDeg()); }
    // bool IsInEId(const int& EId) const { return InEIdV.SearchBin(EId) != -1; }
    // bool IsOutEId(const int& EId) const { return OutEIdV.SearchBin(EId) != -1; }
    friend class TMMNet;
  };
  
  class TEdge {
  private:
    TInt Id, SrcNId, DstNId, Type;
  public:
    TEdge() : Id(-1), SrcNId(-1), DstNId(-1), Type(-1) { }
    TEdge(const int& EId, const int& SourceNId, const int& DestNId, const int& EType) : Id(EId), SrcNId(SourceNId), DstNId(DestNId), Type(EType) { }
    TEdge(const TEdge& Edge) : Id(Edge.Id), SrcNId(Edge.SrcNId), DstNId(Edge.DstNId), Type(Edge.Type) { }
    // TEdge(TSIn& SIn) : Id(SIn), SrcNId(SIn), DstNId(SIn) { }
    // void Save(TSOut& SOut) const { Id.Save(SOut); SrcNId.Save(SOut); DstNId.Save(SOut); }
    int GetId() const { return Id; }
    int GetSrcNId() const { return SrcNId; }
    int GetDstNId() const { return DstNId; }
    int GetType() const { return Type; }
    friend class TMMNet;
  };
  
  class TNodeI {
  private:
    typedef THashMP<TInt, TNode>::TIter THashIter;
    THashIter NodeHI;
    TSVNetMP *Net;
  public:
    TNodeI() : NodeHI(), Net(NULL) { }
    TNodeI(const THashIter& NodeHIter, const TSVNetMP* NetPt) : NodeHI(NodeHIter), Net((TSVNetMP *) NetPt) { }
    TNodeI(const TNodeI& NodeI) : NodeHI(NodeI.NodeHI), Net(NodeI.Net) { }
    TNodeI& operator = (const TNodeI& NodeI) { NodeHI=NodeI.NodeHI; Net=NodeI.Net; return *this; }
    /// Increment iterator.
    TNodeI& operator++ (int) { NodeHI++;  return *this; }
    bool operator < (const TNodeI& NodeI) const { return NodeHI < NodeI.NodeHI; }
    bool operator == (const TNodeI& NodeI) const { return NodeHI == NodeI.NodeHI; }
    
    /// Returns ID of the current node.
    int GetId() const { return NodeHI.GetDat().GetId(); }
    /// Returns Node Type of the current node.
    int GetNType() const { return NodeHI.GetDat().GetNType(); }
    /// Returns degree of the current node.
    int GetDeg() const { return NodeHI.GetDat().GetDeg(); }
    /// Returns in-degree of the current node.
    int GetInDeg() const { return NodeHI.GetDat().GetInDeg(); }
    /// Returns out-degree of the current node.
    int GetOutDeg() const { return NodeHI.GetDat().GetOutDeg(); }
    /// Returns in-degree of particular edge type of the current node.
    int GetInDeg(int EType)  { return NodeHI.GetDat().GetInDeg(EType); }
    /// Returns out-degree of particular edge type of the current node.
    int GetOutDeg(int EType)  { return NodeHI.GetDat().GetOutDeg(EType); }
    /// Returns ID of EdgeN-th in-edge of particular edge type.
    int GetInEId(const int& EdgeN, const int& EType) const { return NodeHI.GetDat().GetInEId(EdgeN, EType); }
    /// Returns ID of EdgeN-th out-edge of particular edge type.
    int GetOutEId(const int& EdgeN, const int& EType) const { return NodeHI.GetDat().GetOutEId(EdgeN, EType); }
    
    
    /*/// Returns ID of NodeN-th out-node (the node the current node points to). ##TNodeNet::TNodeI::GetOutNId
     int GetOutNId(const int& NodeN) const { return NodeHI.GetDat().GetOutNId(NodeN); }
     /// Returns ID of NodeN-th neighboring node. ##TNodeNet::TNodeI::GetNbrNId
     int GetNbrNId(const int& NodeN) const { return NodeHI.GetDat().GetNbrNId(NodeN); }
     /// Tests whether node with ID NId points to the current node.
     bool IsInNId(const int& NId) const { return NodeHI.GetDat().IsInNId(NId); }
     /// Tests whether the current node points to node with ID NId.
     bool IsOutNId(const int& NId) const { return NodeHI.GetDat().IsOutNId(NId); }
     /// Tests whether node with ID NId is a neighbor of the current node.
     bool IsNbrNId(const int& NId) const { return IsOutNId(NId) || IsInNId(NId); }
     const TNodeData& operator () () const { return NodeHI.GetDat().NodeDat; }
     TNodeData& operator () () { return NodeHI.GetDat().GetDat(); }
     const TNodeData& GetDat() const { return NodeHI.GetDat().GetDat(); }
     TNodeData& GetDat() { return NodeHI.GetDat().GetDat(); }
     const TNodeData& GetInNDat(const int& NodeN) const { return Net->GetNDat(GetInNId(NodeN)); }
     TNodeData& GetInNDat(const int& NodeN) { return Net->GetNDat(GetInNId(NodeN)); }
     const TNodeData& GetOutNDat(const int& NodeN) const { return Net->GetNDat(GetOutNId(NodeN)); }
     TNodeData& GetOutNDat(const int& NodeN) { return Net->GetNDat(GetOutNId(NodeN)); }
     const TNodeData& GetNbrNDat(const int& NodeN) const { return Net->GetNDat(GetNbrNId(NodeN)); }
     TNodeData& GetNbrNDat(const int& NodeN) { return Net->GetNDat(GetNbrNId(NodeN)); }
     friend class TNodeNet<TNodeData>;*/
  };
  
  
private:
  
  
  TVec<THashMP<TInt, TNode> > NodeHV;
  TVec<THashMP<TInt, TEdge> > EdgeHV;
  TIntPrV Mapping;
  TIntV MxNIdV;
  TIntV MxEIdV;
  
  
public:
  typedef THashMP<TInt, TNode>::TIter THashIter;
  typedef TSVNetMP* PSVNetMP;
  
  
  TSVNetMP(): NodeHV(), EdgeHV(), MxNIdV(), MxEIdV() { }
  
  TInt AddNType() {
    THashMP<TInt, TNode> NodeH;
    NodeHV.Add(NodeH);
    MxNIdV.Add(0);
    return NodeHV.Len() - 1;
  }
  
  TInt AddEType(const TInt& SrcNType, const TInt& DstNType) {
    Mapping.Add(TIntPr(SrcNType, DstNType));
    for (int i = 0; i < NodeHV.Len(); i++) {
      for (THashMPKeyDatI<TInt, TNode> it = NodeHV[i].BegI(); it != NodeHV[i].EndI(); it++) {
        it.GetDat().AddEType();
      }
    }
    MxEIdV.Add(0);
    EdgeHV.Add(THashMP<TInt, TEdge>());
    return Mapping.Len() - 1;
  }
  
  int AddNode(int NType, int NId = -1) {
    IAssertR(IsNType(NType), TStr::Fmt("Node type %d does not exist", NType));
    if (NId == -1) {
      NId = MxNIdV[NType];
      MxNIdV[NType]++;
    }
    else {
      IAssertR(!IsNode(NId, NType), TStr::Fmt("NodeId %d, NType %d already exists", NId, NType));
      MxNIdV[NType] = (TInt(NId+1) > MxNIdV[NType])?TInt(NId+1):MxNIdV[NType];
    }
    NodeHV[NType].AddDat(NId, TNode(NId, NType, Mapping.Len()));
    return NId;
  }
  
  int AddEdge(int SrcNId, int DstNId, int EType, int EId = -1) {
    IAssertR(IsEType(EType), TStr::Fmt("Edge type %d does not exist", EType)); //
    int SrcNType = Mapping[EType].Val1;
    int DstNType = Mapping[EType].Val2;
    if (! IsNode(SrcNId, SrcNType)) {
      AddNode(SrcNType, SrcNId);
    }
    if (! IsNode(DstNId, DstNType)) {
      AddNode(DstNType, DstNId);
    }
    //IAssertR(IsNode(SrcNId, SrcNType), TStr::Fmt("Node id %d, Node type %d does not exist", SrcNId, SrcNType));
    //IAssertR(IsNode(DstNId, DstNType), TStr::Fmt("Node id %d, Node type %d does not exist", DstNId, DstNType));
    if (EId == -1) {
      EId = MxEIdV[EType];
      MxEIdV[EType]++;
    }
    else {
      IAssertR(!IsEdge(EId, EType), TStr::Fmt("EdgeId %d, EType %d already exists", EId, EType)); //
      MxEIdV[EType] = (TInt(EId+1) > MxEIdV[EType])?TInt(EId+1):MxEIdV[EType];
    }
    EdgeHV[EType].AddDat(EId, TEdge(EId, SrcNId, DstNId, EType));
    NodeHV[SrcNType].GetDat(SrcNId).AddOutEdge(EId, EType);
    NodeHV[DstNType].GetDat(DstNId).AddInEdge(EId, EType);
    return EId;
  }
  
  bool IsNode(const int& NId, const int& NType) { return NodeHV[NType].IsKey(NId); }
  bool IsNType(const int& NType) { return Mapping.Len() > NType; }
  bool IsEdge(const int& EId, const int& EType) { return EdgeHV[EType].IsKey(EId); }
  bool IsEType(int EType) { return Mapping.Len() > EType; }
  
  int GetNodes() {
    int NodeCnt = 0;
    for (int i = 0; i < NodeHV.Len(); i++) {
      NodeCnt += NodeHV[i].Len();
    }
    return NodeCnt;
  }
  int GetEdges() {
    int EdgeCnt = 0;
    for (int i = 0; i < EdgeHV.Len(); i++) {
      EdgeCnt += EdgeHV[i].Len();
    }
    return EdgeCnt;
  }
  int GetNodes(int NType) { return NodeHV[NType].Len(); }
  int GetEdges(int EType) { return EdgeHV[EType].Len(); }
  int GetNTypeCnt() { return NodeHV.Len(); }
  int GetETypeCnt() { return EdgeHV.Len(); }
  
  TNodeI BegNI(const int& NType) {
    IAssertR(IsNType(NType), TStr::Fmt("Node type %d does not exist", NType));
    return TNodeI(NodeHV[NType].BegI(), this);
  }
  
  /// Returns an iterator referring to the past-the-end node in the network.
  TNodeI EndNI(const int& NType) {
    IAssertR(IsNType(NType), TStr::Fmt("Node type %d does not exist", NType));
    return TNodeI(NodeHV[NType].EndI(), this);
  }
  int GetSrcNType(const int& EType) { return Mapping[EType].Val1; }
  int GetDstNType(const int& EType) { return Mapping[EType].Val2; }
  int GetSrcNId(const int& EId, const int& EType) { return EdgeHV[EType].GetDat(EId).GetSrcNId(); }
  int GetDstNId(const int& EId, const int& EType) { return EdgeHV[EType].GetDat(EId).GetDstNId(); }
  //TNode& GetNI(const int& NId, const int& NType) { return NodeHV[NType].GetDat(NId); }
  /// Returns an iterator referring to the node of ID NId in the network.
  TNodeI GetNI(const int& NId, const int& NType) const { return TNodeI(NodeHV[NType].GetI(NId), this); }
  int GetMxNId(const int& NType) { return MxNIdV[NType] - 1; }
  int GetMxEId(const int& EType) { return MxEIdV[EType] - 1; }
  void SetMxNId(const int& NType, const TInt& MxNId) { MxNIdV[NType] = MxNId; }
  void SetMxEId(const int& EType, const TInt& MxEId) { MxEIdV[EType] = MxEId; }
  void SetNCnt(const TInt& NType, const TInt& L) { NodeHV[NType].SetLen(L); }
  void SetECnt(const TInt& EType, const TInt& L) { EdgeHV[EType].SetLen(L); }
  
  /*
   void GetPageRankMM(const PSVNetMP& Graph, TVec<TIntFltH>& PRankHV, const double& C, const double& Eps, const int& MaxIter);
   */
  
  void ReserveNodes(const int& NType, const int& Nodes) { if(Nodes > 0) { NodeHV[NType].Gen(Nodes); } }
  void ReserveEdges(const int& EType, const int& Edges) { if(Edges > 0) { EdgeHV[EType].Gen(Edges); } }
  //void SetEdges(const int& EType, )
  void AddEdgeToHash(const TInt& EId, const TInt& SrcNId, const TInt& DstNId, const TInt& EType) {
    int EdgeIdx = abs((EId.GetPrimHashCd()) % EdgeHV[EType].GetReservedKeyIds());
    int EdgeKeyId = EdgeHV[EType].AddKey13(EdgeIdx, EId);
    EdgeHV[EType][EdgeKeyId] = TEdge(EId, SrcNId, DstNId, EType);
  }
  void AddNodeWithEdges(const TInt& NId, const TInt& NType, TVec<TIntV>& InEIdVV, TVec<TIntV>& OutEIdVV) {
    int NodeIdx = abs((NId.GetPrimHashCd()) % NodeHV[NType].GetReservedKeyIds());
    int NodeKeyId = NodeHV[NType].AddKey13(NodeIdx, NId);
    int ETypeCnt = InEIdVV.Len();
    NodeHV[NType][NodeKeyId] = TNode(NId, NType, 0);
    //NodeHV[NType][NodeKeyId].SetInEIdVV(InEIdVV);
    //NodeHV[NType][NodeKeyId].SetOutEIdVV(OutEIdVV);
  }
  
  
  PSVNetMP GetSubGraph(TIntV NTypeV, TIntV ETypeV);
  PNEANet GetSubGraphTNEANet(TIntV NTypeV, TIntV ETypeV, TIntIntH& Offsets);
  
  #ifdef _OPENMP
  PSVNetMP GetSubGraphMP(TIntV NTypeV, TIntV ETypeV);
  #endif
  
};

typedef TSVNetMP* PSVNetMP;
void GetPageRankMM(const PSVNetMP& Graph, TVec<TIntFltH>& PRankHV, const double& C, const double& Eps, const int& MaxIter);
void GetBfsLevelMM(const PSVNetMP& Graph, TVec<TIntIntH>& BfsHV, const int& StartNId, const int& StartNType);
int GetBfsLevelMMMP2(const PSVNetMP& Graph, TVec<TIntV >& BfsLevelVV, const int& StartNId, const int& StartNType);
//void GetSimRankMM(const PSVNetMP& Graph, THashMP<TPair<TIntPr, TIntPr>, TFlt>& SRankH, const double& C, const double& Eps, const int& MaxIter);



#ifdef _OPENMP
void GetPageRankMMMP2(const PSVNetMP& Graph, TVec<TIntFltH>& PRankHV, const double& C, const double& Eps, const int& MaxIter);
void GetBfsLevelMMMP(const PSVNetMP& Graph, TVec<TIntIntH>& BfsHV, const int& StartNId, const int& StartNType);
#endif



#endif
