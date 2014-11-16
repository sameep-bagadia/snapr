#include "Snap.h"
#include <sys/time.h>
#include <sys/resource.h>
#include <stdio.h>
//#include <proc/readproc.h>

typedef TVec<PTable> PTableV;
typedef TVec<TIntPr > TIntPrV;


/**************************Benchmark Functions*****************************/

float getcputime() {
  struct rusage rusage;
  float result;
  
  getrusage(RUSAGE_SELF, &rusage);
  
  result =
  ((float) (rusage.ru_utime.tv_usec + rusage.ru_stime.tv_usec) / 1000000) +
  ((float) (rusage.ru_utime.tv_sec + rusage.ru_stime.tv_sec));
  return result;
}
/*
void getcpumem(float *scpu, float *smem) {
  struct rusage rusage;
  struct proc_t usage;
  
  getrusage(RUSAGE_SELF, &rusage);
  *scpu =
  ((float) (rusage.ru_utime.tv_usec + rusage.ru_stime.tv_usec) / 1000000) +
  ((float) (rusage.ru_utime.tv_sec + rusage.ru_stime.tv_sec));
  
  look_up_our_self(&usage);
  *smem = (float) usage.vsize / 1000000;
}

void getmaxcpumem(float *scpu, float *smem) {
  struct rusage rusage;
  
  getrusage(RUSAGE_SELF, &rusage);
  *scpu =
  ((float) (rusage.ru_utime.tv_usec + rusage.ru_stime.tv_usec) / 1000000) +
  ((float) (rusage.ru_utime.tv_sec + rusage.ru_stime.tv_sec));
  *smem = (float) (rusage.ru_maxrss) / 1000;
}
*/
/***************************************************************************/


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
  char FileName[200];
  int ColCnt;
  int Reverse;
  printf("Adding Edge Table\n");
  printf("Enter filename, number of columns (>= 2), and whether reverse? (reverse = 1, not reverse = 0\n");
  scanf("%s %d %d", FileName, &ColCnt, &Reverse);
  Schema EdgeScm;
  if (Reverse == 1) {
    EdgeScm.Add(TPair<TStr, TAttrType>("DstID", atStr));
    EdgeScm.Add(TPair<TStr, TAttrType>("SrcID", atStr));
  }
  else {
    EdgeScm.Add(TPair<TStr, TAttrType>("SrcID", atStr));
    EdgeScm.Add(TPair<TStr, TAttrType>("DstID", atStr));
  }
  for (TInt i = 1; i < ColCnt-1; i++) {
    TStr ColName = "Attribute" + i.GetStr();
    EdgeScm.Add(TPair<TStr, TAttrType>(ColName, atStr));
  }
  TStr FName(FileName);
  PTable T = TTable::LoadSS(EdgeScm, FName, Context);
  return T;
}


//one iteration of pagerank algorithm
void IteratePRankMP(const PTableV& NTables, const PTableV& ETables, const TIntPrV& Mapping, const TFltV Weights,
                  const TVec<THash<TGroupKey, TPair<TInt, TIntV> > >& GroupMap, TVec<TFltV>& PRank, int NodeCnt) {
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
  int NType, NodeID;
  #pragma omp parallel for collapse(2) schedule(dynamic, 10000)
  for (NType = 0; NType < NTblCnt; NType++) {
    //printf("Starting node type %d\n", NType);
    int RowCnt = (NTables[NType]->GetNumRows()).Val;
    for (NodeID = 0; NodeID < (NTables[NType]->GetNumRows()).Val; NodeID++) {
      //printf("NType:%d, \tnode %d / %d\n", NType, NodeID, RowCnt);
      PRank[NType][NodeID] = (1 - d) / float(NodeCnt);
      TGroupKey NodeKey;
      NodeKey.Val1.Add(NodeID);
      //check all the incoming edges to the node
      
      
      for (int EType = 0; EType < ETblCnt; EType++) {
        //node type should match the destination node for edge type
        if (Mapping[EType].Val2 == NType) {
          //Check if given node has edges
          if (GroupMap[EType].IsKey(NodeKey)) {
            //Get the vector of inedges
            const TIntV & InEdges = GroupMap[EType].GetDat(NodeKey).Val2;
            //Iterate through each inedge to add to pagerank score of the node
            for (int i = 0; i < InEdges.Len(); i++) {
              int EId = InEdges[i];
              int InNodeId = (ETables[EType]->GetIntVal("SrcTTID", EId)).Val;
              PRank[NType][NodeID] += (d * Weights[EType] * OldPR[Mapping[EType].Val1][InNodeId]) / (NTables[Mapping[EType].Val1]->GetFltVal("OutWt", InNodeId)).Val;
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

void GetPagerankMMMP(const PTableV& NTables,const PTableV& ETables,const TIntPrV& Mapping, const TFltV Weights, const TVec<THash<TGroupKey, TPair<TInt, TIntV> > >& GroupMap, PTableV& Result) {
  TVec<TFltV> PRank;
  int NTblCnt = NTables.Len();
  int ETblCnt = ETables.Len();
  int NodeCnt = 0;
  for (int i = 0; i < NTblCnt; i++) {
    NodeCnt += NTables[i]->GetNumRows();
  }
  printf("%d\n", NodeCnt);
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
  
  //Run pagerank for some iterations
  int IterCnt = 10;
  //PrintPRank(PRank);
  for (int i = 0; i < IterCnt; i++) {
    printf("Starting pagerank iteration number: %d\n", i);
    IteratePRankMP(NTables, ETables, Mapping, Weights, GroupMap, PRank, NodeCnt);
    //PrintPRank(PRank);
  }
  /*
  //Print the final pagerank values
  printf("Final Page rank values:\n");
  PrintPRank(PRank);
  */
  //Store the pagerank values in Result
  Result.Clr();
  for (int i = 0; i < NTblCnt; i++) {
    PTable NewPTable = new TTable(*NTables[i]);
    NewPTable->StoreFltCol("Pagerank", PRank[i]);
    Result.Add(NewPTable);
  }
}

void StoreTTableIds(const PTableV& NTables, PTableV& ETables, const TIntPrV& Mapping) {
  int NTblCnt = NTables.Len();
  int ETblCnt = ETables.Len();
  
  //Create the hashtable from NodeID to index
  TVec<TStrIntH > NodeHash;
  for (int i = 0; i < NTblCnt; i++) {
    TStrIntH hash;
    for (int j = 0; j < NTables[i]->GetNumRows().Val; j++) {
      hash.AddDat(NTables[i]->GetStrVal("NodeID", j), j);
    }
    NodeHash.Add(hash);
  }
  
  //Store the TTable Ids in the edge tables
  for (int EType = 0; EType < ETblCnt; EType++) {
    TIntV SrcTTId;
    TIntV DstTTId;
    int SrcType = Mapping[EType].Val1;
    int DstType = Mapping[EType].Val2;
    for (int EId = 0; EId < (ETables[EType]->GetNumRows()).Val; EId++) {
      SrcTTId.Add(NodeHash[SrcType].GetDat(ETables[EType]->GetStrVal("SrcID" ,EId)));
      DstTTId.Add(NodeHash[DstType].GetDat(ETables[EType]->GetStrVal("DstID" ,EId)));
    }
    ETables[EType]->StoreIntCol("SrcTTID", SrcTTId);
    ETables[EType]->StoreIntCol("DstTTID", DstTTId);
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
      TStr Key1 = ETables[i]->GetStrVal("SrcID", j);
      int tempval1 = OutDeg[SrcTbl].GetDat(Key1);
      OutDeg[SrcTbl].AddDat(Key1, tempval1+1);
      //printf("nodes%d, id = %s, new_val = %d\n", SrcTbl, (ETables[i]->GetStrVal("SrcID", j)).CStr(), (OutDeg[SrcTbl].GetDat(Key1)).Val);
      
      TStr Key2 = ETables[i]->GetStrVal("DstID", j);
      int tempval2 = OutDeg[DstTbl].GetDat(Key2);
      OutDeg[DstTbl].AddDat(Key2, tempval2+1);
      
      //printf("nodes%d, id = %s, new_val = %d\n", DstTbl, (ETables[i]->GetStrVal("DstID", j)).CStr(), (OutDeg[DstTbl].GetDat(Key2)).Val);
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

// Calculate the outweights of each node and stores as a new attribute in the node tables.
void CalcOutWt(PTableV& NTables, const PTableV& ETables, const TIntPrV& Mapping, const TFltV& Weights) {
  printf("Calculating out weights\n");
  TVec<TStrFltH > OutWt;
  int NTblCnt = NTables.Len();
  int ETblCnt = ETables.Len();
  
  //Initialize the vector to 0
  for (int i = 0; i < NTblCnt; i++) {
    TStrFltH WtHash;
    int RowCnt = (NTables[i]->GetNumRows()).Val;
    for (int j = 0; j < RowCnt; j++) {
      WtHash.AddDat(NTables[i]->GetStrVal("NodeID", j), 0);
    }
    OutWt.Add(WtHash);
  }
  // Calculate the outweights
  
  for (int i = 0; i < ETblCnt; i++) {
    int SrcTbl = (Mapping[i].Val1).Val;
    int DstTbl = (Mapping[i].Val2).Val;
    printf("Src table = %d, dest table = %d\n", SrcTbl, DstTbl);
    int EdgeCnt = (ETables[i]->GetNumRows()).Val;
    printf("# of edges = %d\n", EdgeCnt);
    for (int j = 0; j < EdgeCnt; j++) {
      TStr Key1 = ETables[i]->GetStrVal("SrcID", j);
      float tempval1 = OutWt[SrcTbl].GetDat(Key1);
      OutWt[SrcTbl].AddDat(Key1, tempval1 + Weights[i]);
      //printf("nodes%d, id = %s, new_val = %d\n", SrcTbl, (ETables[i]->GetStrVal("SrcID", j)).CStr(), (OutDeg[SrcTbl].GetDat(Key1)).Val);
    }
  }
  printf("Outweights calculated\n");
  //Store the outweights in the tables.
  for (int i = 0; i < NTblCnt; i++) {
    TFltV fltvec;
    for (int j = 0; j < (NTables[i]->GetNumRows()).Val; j++) {
      fltvec.Add(OutWt[i].GetDat(NTables[i]->GetStrVal("NodeID", j)));
    }
    NTables[i]->StoreFltCol("OutWt", fltvec);
  }
  printf("Outweights stored\n");
}

// Get the groupmap which has the information for in-edges for each node
void GetGroupMap(PTableV& ETables, TVec<THash<TGroupKey, TPair<TInt, TIntV> > >& GroupMap) {
  int ETblCnt = ETables.Len();
  for (int EType = 0; EType < ETblCnt; EType++) {
    THash<TGroupKey, TPair<TInt, TIntV> > Grouping;
    TStrV GroupBy;
    GroupBy.Add("DstTTID");
    ETables[EType]->GetGroupMapping(GroupBy, "", Grouping);
    GroupMap.Add(Grouping);
  }
}

void PrintBenchmarks(FILE* outfile) {
  float cputime;
  cputime = getcputime();
  /*getcpumem(&scpu, &smem);
  getmaxcpumem(&maxscpu, &maxsmem);*/
  fprintf(outfile, "%f\n", cputime);
}



int main(int argc, char* []) {
  
  FILE * outfile;
  outfile = fopen ("benchmark.txt","w");
  
  //outfile << "Hello World!\n";
  fprintf(outfile, "Hello World!\n");
  PrintBenchmarks(outfile);
  
  int NTblCnt;
  int ETblCnt;
  printf("Enter the number of node tables and edge tables\n");
  scanf("%d %d", &NTblCnt, &ETblCnt);
  
  TVec<PTable> NTables;
  TVec<PTable> ETables;
  TVec<PTable> OutNTables;
  TVec<TPair<TInt, TInt> > Mapping;
  TFltV Weights;
  
  TTableContext Context;
  //Read in tables
  for (int i = 0; i < NTblCnt; i++) {
    PTable P = AddNodeTable(Context);
    P->SaveSS("ntable_abc.txt");
    NTables.Add(P);
  }
  for (int i = 0; i < ETblCnt; i++) {
    ETables.Add(AddEdgeTable(Context));
    printf("Enter the source and destination node table number (index starting from 0), and weight of the edges\n");
    int SrcId;
    int DestId;
    float Wt;
    scanf("%d %d %f", &SrcId, &DestId, &Wt);
    printf("%d %d\n", SrcId, DestId);
    Mapping.Add(TPair<TInt, TInt>(SrcId, DestId));
    Weights.Add(Wt);
  }
  
  fprintf(outfile, "Tables Loaded\n");
  PrintBenchmarks(outfile);
    
  //Caculate outweight of each node
  //CalcOutDeg(NTables, ETables, Mapping);
  CalcOutWt(NTables, ETables, Mapping, Weights);
  
  fprintf(outfile,"Out weights calculated\n");
  PrintBenchmarks(outfile);
  
  // Store TTable ids in the edge files
  StoreTTableIds(NTables, ETables, Mapping);
  fprintf(outfile,"TTable ids stored\n");
  PrintBenchmarks(outfile);
  
  printf("TTable ids stored\n");
  
  // Get Group Mapping
  TVec<THash<TGroupKey, TPair<TInt, TIntV> > > GroupMap;
  GetGroupMap(ETables, GroupMap);
  fprintf(outfile,"Group mappings calculated\n");
  PrintBenchmarks(outfile);
  
  printf("Group mappings calculated\n");
  
  int PRCnt;
  printf("Enter the number of times to run pagerank for 10 iterations\n");
  scanf("%d", &PRCnt);
  
  PTableV Result;
  for (int i = 0; i < PRCnt; i++) {
    printf("starting pagerank no %d\n", i);
    //Call Pagerank
    Result.Clr();
    GetPagerankMMMP(NTables, ETables, Mapping, Weights, GroupMap, Result);
  }
  
  fprintf(outfile,"Pagerank completed\n");
  PrintBenchmarks(outfile);
  
  //Print the pagerank values
  for (int i = 0; i < Result.Len(); i++) {
    /*printf("Node table #%d:\n", i);
    for (int j = 0; j < Result[i]->GetNumRows(); j++) {
    printf("NodeID = %s, Weight = %f, Pagerank score = %f\n", (Result[i]->GetStrVal("NodeID", j)).CStr(), Result[i]->GetFltVal("OutWt", j).Val, (Result[i]->GetFltVal("Pagerank", j)).Val);
    }
    printf("\n");*/
    TInt idx = i;
    TStr outfile = "Pagerank" + idx.GetStr() + ".tsv";
    Result[i]->SaveSS(outfile);
  }
  
  fprintf(outfile,"Tables Saved\n");
  PrintBenchmarks(outfile);
  fclose(outfile);
}

