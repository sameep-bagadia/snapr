#include "Snap.h"
#include <sys/time.h>
#include <sys/resource.h>
#include <iostream>
#include <stdio.h>
#include <proc/readproc.h>

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
              int InNodeId = (NodeHash[Mapping[EType].Val2].GetDat(InNode)).Val;
              PRank[NType][NodeID] += (d * OldPR[Mapping[EType].Val2][InNodeId]) / (NTables[Mapping[EType].Val2]->GetIntVal("OutDeg", InNodeId)).Val;
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
  //PrintPRank(PRank);
  for (int i = 0; i < IterCnt; i++) {
    printf("Starting pagerank iteration number: %d\n", i);
    IteratePRank(NTables, ETables, Mapping, PRank, NodeHash, NodeCnt);
    //PrintPRank(PRank);
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

void PrintBenchmarks(FILE* outfile) {
  float scpu, smem, maxscpu, maxsmem, cputime;
  cputime = getcputime();
  getcpumem(&scpu, &smem);
  getmaxcpumem(&maxscpu, &maxsmem);
  fprintf(outfile, "%f %f %f %f %f\n", cputime, scpu, smem, maxscpu, maxsmem);
}



int main(int argc, char* []) {
  
  //ofstream outfile;
  //outfile.open("benchmark.txt");
  
  FILE * outfile;
  outfile = fopen ("benchmark.txt","w");
  
  //outfile << "Hello World!\n";
  fprintf("Hello World!\n");
  PrintBenchmarks(outfile);
  
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
  for (int i = 0; i < ETblCnt; i++) {
    ETables.Add(AddEdgeTable(Context));
    printf("Enter the source and destination node table number (index starting from 0)\n");
    int SrcId;
    int DestId;
    scanf("%d %d", &SrcId, &DestId);
    printf("%d %d\n", SrcId, DestId);
    Mapping.Add(TPair<TInt, TInt>(SrcId, DestId));
  }

  fprintf(outfile, "Tables Loaded\n");
  PrintBenchmarks(outfile);
  /* TEST
   
   
   
   NTables[0]->SaveSS("ntable0.txt");
   
   printf("check: %s\n", NTables[0]->GetStrVal("NodeID", 2).CStr());
   */
  
  ///////////////
  
  PTableV Result;
  
  //Caculate outdegrees of each node
  CalcOutDeg(NTables, ETables, Mapping);
  
  fprintf(outfile,"Degree calculated\n");
  PrintBenchmarks(outfile);
  
  printf("starting pagerank\n");
  //Call Pagerank
  GetPagerankMM(NTables, ETables, Mapping, Result);
  
  fprintf(outfile,"Pagerank completed\n");
  PrintBenchmarks(outfile);
  
  //Print the pagerank values
  for (int i = 0; i < Result.Len(); i++) {
    printf("Node table #%d:\n", i);
    for (int j = 0; j < Result[i]->GetNumRows(); j++) {
      printf("NodeID = %s, Degree = %d, Pagerank score = %f\n", (Result[i]->GetStrVal("NodeID", j)).CStr(), Result[i]->GetIntVal("OutDeg", j).Val, (Result[i]->GetFltVal("Pagerank", j)).Val);
    }
    printf("\n");
    TInt idx = i;
    TStr outfile = "Pagerank" + idx.GetStr() + ".tsv";
    Result[i]->SaveSS(outfile);
  }
  
  fprintf(outfile,"Tables Saved\n");
  PrintBenchmarks(outfile);
  fclose(outfile);
}

