#include "Snap.h"
//using namespace TSnap;
typedef TVec<PTable> PTableV;
typedef TVec<TIntPr > TIntPrV;


//Function to read a table of nodes
PTable AddNodeTable(TTableContext& Context) {
  Schema NodeScm;
  NodeScm.Add(TPair<TStr, TAttrType>("NodeID", atStr));
  char FileName[50];
  int ColCnt = 0;
  printf("Adding Node Table\n");
  printf("Enter filename and number of columns (>= 1) \n");
  scanf("%s %d", FileName, &ColCnt);
  for (TInt i = 1; i < ColCnt; i++) {
    TStr ColName = "Attribute" + i.GetStr();
    NodeScm.Add(TPair<TStr, TAttrType>(ColName, atStr));
  }
  TStr FName(FileName);
  PTable T = TTable::LoadSS(NodeScm, FName, Context);
  return T;
}

// Function to read in a table of edges
PTable AddEdgeTable(TTableContext& Context) {
  Schema EdgeScm;
  EdgeScm.Add(TPair<TStr, TAttrType>("SrcID", atStr));
  EdgeScm.Add(TPair<TStr, TAttrType>("DstID", atStr));
  char FileName[200];
  int ColCnt;
  printf("Adding Edge Table\n");
  printf("Enter filename and number of columns (>= 2)\n");
  scanf("%s %d", FileName, &ColCnt);
  for (TInt i = 1; i < ColCnt-1; i++) {
    TStr ColName = "Attribute" + i.GetStr();
    EdgeScm.Add(TPair<TStr, TAttrType>(ColName, atStr));
  }
  TStr FName(FileName);
  PTable T = TTable::LoadSS(EdgeScm, FName, Context);
  return T;
}

// Function to read in a table of edges in reverse direction as it is undirected graph. So each edge is read in both direction
PTable AddReverseEdgeTable(TTableContext & Context) {
  //TTableContext Context;
  Schema Edges;
  //Just named the first list as dstid and second one as srcid to solve the purpose of reversing
  Edges.Add(TPair<TStr, TAttrType>("DstID", atStr));
  Edges.Add(TPair<TStr, TAttrType>("SrcID", atStr));
  char filename[200];
  int col_count;
  printf("Adding Edge Table\n");
  printf("Enter filename and number of columns (>= 2)\n");
  scanf("%s %d", filename, &col_count);
  for (TInt i = 1; i < col_count-1; i++) {
    TStr ColName = "Attribute" + i.GetStr();
    Edges.Add(TPair<TStr, TAttrType>(ColName, atStr));
  }
  TStr FName(filename);
  PTable T = TTable::LoadSS(Edges, FName, Context);
  return T;
}

//one iteration of pagerank algorithm
void IteratePRank(const PTableV& NTables, const PTableV& ETables, const TIntPrV& Mapping,
                  TVec<TFltV>& PRank, const TVec<TStrIntH >& NodeHash, int NodeCnt) {
  //Make a copy of page rank values
  float d = 0.85;
  TVec<TFltV> OldPR;
  for (int i = 0; i < PRank.Len(); i++) {
    TFltV FltVec;
    for (int j = 0; j < PRank[i].Len(); j++) {
      FltVec.Add(PRank[i][j]);
    }
    OldPR.Add(FltVec);
  }
  
  // update pagerank scores
  int NTblCnt = NTables.Len();
  int ETblCnt = ETables.Len();
  
  // for each node in the graph
  for (int NType = 0; NType < NTblCnt; NType++) {
    printf("Starting node type %d\n", NType);
    int RowCnt = (NTables[NType]->GetNumRows()).Val;
    for (int NodeID = 0; NodeID < (NTables[NType]->GetNumRows()).Val; NodeID++) {
      printf("NType:%d, \tnode %d / %d\n", NType, NodeID, RowCnt);
      TStr CurrNode = NTables[NType]->GetStrVal("NodeID", NodeID);
      PRank[NType][NodeID] = (1 - d) / float(NodeCnt);
      
      //check all the incoming edges to the node
      for (int EType = 0; EType < ETblCnt; EType++) {
        //node type should match the destination node for edge type
        if (Mapping[EType].Val2 == NType) {
          for (int EId = 0; EId < (ETables[EType]->GetNumRows()).Val; EId++) {
            if (CurrNode == ETables[EType]->GetStrVal("DstID", EId)) {
              //add to the pagerank score
              TStr InNode = ETables[EType]->GetStrVal("SrcID", EId);
              int InNodeId = (NodeHash[Mapping[EType].Val1].GetDat(InNode)).Val;
              PRank[NType][NodeID] += (d * OldPR[Mapping[EType].Val1][InNodeId]) / (NTables[Mapping[EType].Val1]->GetIntVal("OutDeg", InNodeId)).Val;
            }
          }
        }
        if (Mapping[EType].Val1 == NType) {
          for (int EId = 0; EId < (ETables[EType]->GetNumRows()).Val; EId++) {
            if (CurrNode == ETables[EType]->GetStrVal("SrcID", EId)) {
              //add to the pagerank score
              TStr InNode = ETables[EType]->GetStrVal("DstID", EId);
              int InNodeId = (NodeHash[Mapping[EType].Val1].GetDat(InNode)).Val;
              PRank[NType][NodeID] += (d * OldPR[Mapping[EType].Val1][InNodeId]) / (NTables[Mapping[EType].Val1]->GetIntVal("OutDeg", InNodeId)).Val;
            }
          }
        }
      }
    }
  }
}

void PrintPRank(const TVec<TFltV> & PRank) {
  for (int i = 0; i < PRank.Len(); i++) {
    for (int j = 0; j < PRank[i].Len(); j++) {
      printf("%lf ", (PRank[i][j]).Val);
    }
    printf("\n");
  }
  printf("\n");
}

void GetPagerankMM(const PTableV& NTables,const PTableV& ETables,const TIntPrV& Mapping, PTableV& Result) {
  TVec<TFltV> PRank;
  int NTblCnt = NTables.Len();
  int ETblCnt = ETables.Len();
  int NodeCnt = 0;
  for (int i = 0; i < NTblCnt; i++) {
    NodeCnt += NTables[i]->GetNumRows();
  }
  printf("%d\n", NodeCnt);
  return;
  float PRInit = 1.0 / float(NodeCnt);
  
  //Initialize pagerank value for each node as 1/N
  for (int i = 0; i < NTblCnt; i++) {
    TFltV FltVec;
    int RowCnt = (NTables[i]->GetNumRows()).Val;
    for (int j = 0; j < RowCnt; j++) {
      FltVec.Add(PRInit);
    }
    PRank.Add(FltVec);
  }
  
  //Create the hashtable from NodeID to index
  TVec<TStrIntH > NodeHash;
  for (int i = 0; i < NTblCnt; i++) {
    TStrIntH hash;
    for (int j = 0; j < NTables[i]->GetNumRows().Val; j++) {
      hash.AddDat(NTables[i]->GetStrVal("NodeID", j), j);
    }
    NodeHash.Add(hash);
  }
  
  //Run pagerank for some iterations
  int IterCnt = 10;
  for (int i = 0; i < IterCnt; i++) {
    printf("Starting pagerank iteration number: %d\n", i);
    IteratePRank(NTables, ETables, Mapping, PRank, NodeHash, NodeCnt);
  }
  
  //Print the final pagerank values
  printf("Final Page rank values:\n");
  PrintPRank(PRank);
  
  //Store the pagerank values in Result
  Result.Clr();
  for (int i = 0; i < NTblCnt; i++) {
    PTable NewPTable = new TTable(*NTables[i]);
    NewPTable->StoreFltCol("Pagerank", PRank[i]);
    Result.Add(NewPTable);
  }
}


// Calculate the outdegrees of each node and stores as a new attribute in the node tables.
void CalcOutDeg(PTableV& NTables, const PTableV& ETables, const TIntPrV& Mapping) {
  printf("Calculating out degrees\n");
  TVec<TStrIntH > OutDeg;
  int NTblCnt = NTables.Len();
  int ETblCnt = ETables.Len();
  
  //Initialize the vector to 0
  for (int i = 0; i < NTblCnt; i++) {
    TStrIntH DegHash;
    int RowCnt = (NTables[i]->GetNumRows()).Val;
    for (int j = 0; j < RowCnt; j++) {
      DegHash.AddDat(NTables[i]->GetStrVal("NodeID", j), 0);
    }
    OutDeg.Add(DegHash);
  }
  // Calculate the outdegrees
  
  for (int i = 0; i < ETblCnt; i++) {
    int SrcTbl = (Mapping[i].Val1).Val;
    int DstTbl = (Mapping[i].Val2).Val;
    printf("Src table = %d, dest table = %d\n", SrcTbl, DstTbl);
    int EdgeCnt = (ETables[i]->GetNumRows()).Val;
    printf("# of edges = %d\n", EdgeCnt);
    for (int j = 0; j < EdgeCnt; j++) {
      TStr Key = ETables[i]->GetStrVal("SrcID", j);
      int tempval = OutDeg[SrcTbl].GetDat(Key);
      //OutDeg[SrcTbl].DelKey(Key);
      OutDeg[SrcTbl].AddDat(Key, tempval+1);
      
      //printf("nodes%d, id = %d, new_val = %d\n", SrcTbl, (ETables[i]->GetIntVal("SrcID", j)).Val, (OutDeg[SrcTbl].GetDat(Key)).Val);
    }
  }
  printf("Outdegrees calculated\n");
  //Store the outdegrees in the tables.
  for (int i = 0; i < NTblCnt; i++) {
    TIntV intvec;
    for (int j = 0; j < (NTables[i]->GetNumRows()).Val; j++) {
      intvec.Add(OutDeg[i].GetDat(NTables[i]->GetStrVal("NodeID", j)));
    }
    NTables[i]->StoreIntCol("OutDeg", intvec);
  }
  printf("Outdegrees stored\n");
}



int main(int argc, char* []) {
  
  int NTblCnt;
  int ETblCnt;
  printf("Enter the number of node tables and edge tables\n");
  scanf("%d %d", &NTblCnt, &ETblCnt);
  
  TVec<PTable> NTables;
  TVec<PTable> ETables;
  TVec<PTable> OutNTables;
  TVec<TPair<TInt, TInt> > Mapping;
  
  TTableContext Context;
  //Read in tables
  for (int i = 0; i < NTblCnt; i++) {
    PTable P = AddNodeTable(Context);
    P->SaveSS("ntable_abc.txt");
    NTables.Add(P);
  }
  NTables[0]->SaveSS("ntabledef.txt");
  for (int i = 0; i < ETblCnt; i++) {
    ETables.Add(AddEdgeTable(Context));
    //ETables.Add(AddReverseEdgeTable(Context));
    printf("Enter the source and destination node table number (index starting from 0)\n");
    int SrcId, DstId;
    scanf("%d %d", &SrcId, &DstId);
    Mapping.Add(TIntPr(SrcId, DstId));
    Mapping.Add(TIntPr(DstId, SrcId));
  }
  ETblCnt *= 2;
  
  /* TEST
   
   
   
   NTables[0]->SaveSS("ntable0.txt");
   
   printf("check: %s\n", NTables[0]->GetStrVal("NodeID", 2).CStr());
   */
  
  //*/////////////
  
  PTableV Result;
  
  //Caculate outdegrees of each node
  CalcOutDeg(NTables, ETables, Mapping);
  
  printf("starting pagerank\n");
  //Call Pagerank
  GetPagerankMM(NTables, ETables, Mapping, Result);
  
  //Print the pagerank values
  for (int i = 0; i < Result.Len(); i++) {
    printf("Node table #%d:\n", i);
    for (int j = 0; j < Result[i]->GetNumRows(); j++) {
      printf("NodeID = %s, Pagerank score = %f\n", (Result[i]->GetStrVal("NodeID", j)).CStr(), (Result[i]->GetFltVal("Pagerank", j)).Val);
    }
    printf("\n");
    TInt idx = i;
    TStr outfile = "Pagerank" + idx.GetStr() + ".tsv";
    Result[i]->SaveSS(outfile);
  }
}

