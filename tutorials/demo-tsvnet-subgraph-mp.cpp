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

void PrintBenchmarks(FILE* outfile) {
  float cputime;
  cputime = getcputime();
  /*getcpumem(&scpu, &smem);
   getmaxcpumem(&maxscpu, &maxsmem);*/
  //fprintf(outfile, "%f\n", cputime);
  fprintf(outfile, "%f \t%f\n", cputime, omp_get_wtime());
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
  //TFltV Weights;
  
  TTableContext Context;
  //Read in tables
  for (int i = 0; i < NTblCnt; i++) {
    PTable P = AddNodeTable(Context);
    P->SaveSS("ntable_abc.txt");
    NTables.Add(P);
  }
  for (int i = 0; i < ETblCnt; i++) {
    ETables.Add(AddEdgeTable(Context));
    //printf("Enter the source and destination node table number (index starting from 0), and weight of the edges\n");
    printf("Enter the source and destination node table number (index starting from 0)\n");
    int SrcId;
    int DestId;
    //float Wt;
    //scanf("%d %d %f", &SrcId, &DestId, &Wt);
    scanf("%d %d", &SrcId, &DestId);
    printf("%d %d\n", SrcId, DestId);
    Mapping.Add(TPair<TInt, TInt>(SrcId, DestId));
    //Weights.Add(Wt);
  }
  
  /*
   char nodes_file[200];
   int NodesCnt;
   printf("Enter the filename containing bfs nodes and number of nodes to take");
   scanf("%s %d", nodes_file, &NodesCnt);
   
   FILE* infile = fopen (nodes_file,"r");
   char NodeTypeStr[100];
   char NodeIdStr[200];
   TIntV NodeTypeV;
   TStrV NodeIdV;
   for (int i = 0; i < NodesCnt; i++) {
   fscanf(infile, "%s %s", NodeTypeStr, NodeIdStr);
   int NodeTypeInt;
   TStr NodeIdStr2 = TStr(NodeIdStr);
   if (NodeTypeStr[0] == 'c') { NodeTypeInt = 0; }
   else if (NodeTypeStr[0] == 'l') { NodeTypeInt = 1; }
   else if (NodeTypeStr[0] == 'p') { NodeTypeInt = 2; }
   else if (NodeTypeStr[0] == 't') { NodeTypeInt = 3; }
   else if (NodeTypeStr[0] == 'u') { NodeTypeInt = 4; }
   else { printf("ERROR! Node Type not found!\n"); }
   NodeTypeV.Add(NodeTypeInt);
   NodeIdV.Add(NodeIdStr2);
   
   }
   */
  
  fprintf(outfile, "Tables Loaded\n");
  PrintBenchmarks(outfile);
  printf("Tables Loaded\n");
  
  //Convert to Graph
  printf("Converting to Graph\n");
  TSVNetMP Graph;
  TIntV NTypeV;
  TIntV ETypeV;
  TVec<THash<TStr, TInt> > NodesHV;
  //Adding node types
  for (int i = 0; i < NTblCnt; i++) {
    NTypeV.Add(Graph.AddNType());
  }
  printf("a\n");
  //Adding edge types
  for (int i = 0; i < ETblCnt; i++) {
    ETypeV.Add(Graph.AddEType(NTypeV[Mapping[i].Val1], NTypeV[Mapping[i].Val2]));
  }
  printf("a\n");
  //Adding nodes
  for (int i = 0; i < NTblCnt; i++) {
    THash<TStr, TInt> NodesH;
    NodesHV.Add(NodesH);
    int NType = NTypeV[i];
    Graph.ReserveNodes(i, NTables[i]->GetNumRows());
    for (int j = 0; j < NTables[i]->GetNumRows(); j++) {
      NodesHV[i].AddDat(NTables[i]->GetStrVal("NodeID", j), Graph.AddNode(NType));
    }
    Graph.SetNCnt(i, NTables[i]->GetNumRows());
  }
  //printf("a\n");
  //Adding edges
  for (int i = 0; i < ETblCnt; i++) {
    int EType = ETypeV[i];
    int SrcNType = NTypeV[Mapping[i].Val1];
    int DstNType = NTypeV[Mapping[i].Val2];
    Graph.ReserveEdges(i, ETables[i]->GetNumRows());
    for (int j = 0; j < ETables[i]->GetNumRows(); j++) {
      int SrcNId = NodesHV[SrcNType].GetDat(ETables[i]->GetStrVal("SrcID", j));
      int DstNId = NodesHV[DstNType].GetDat(ETables[i]->GetStrVal("DstID", j));
      Graph.AddEdge(SrcNId, DstNId, EType);
    }
    Graph.SetECnt(i, ETables[i]->GetNumRows());
  }
  //printf("a\n");
  
  TIntV EdgeVec1;
  //EdgeVec.Add(0);
  //EdgeVec.Add(2);
  TIntV NodeVec1;
  NodeVec1.Add(2);
  NodeVec1.Add(4);
  //TIntIntH Offsets1;
  
  TIntV EdgeVec2;
  //EdgeVec.Add(0);
  //EdgeVec.Add(2);
  TIntV NodeVec2;
  NodeVec2.Add(2);
  NodeVec2.Add(3);
  NodeVec2.Add(4);
  //TIntIntH Offsets2;
  
  
  TIntV EdgeVec3;
  //EdgeVec.Add(0);
  //EdgeVec.Add(2);
  TIntV NodeVec3;
  NodeVec3.Add(0);
  NodeVec3.Add(2);
  NodeVec3.Add(3);
  NodeVec3.Add(4);
  //TIntIntH Offsets3;
  
  int iter_count = 100;
  fprintf(outfile, "Converted to Graph\n");
  PrintBenchmarks(outfile);
  printf("Converted to Graph\n");
  printf("Starting subgraph\n");
  
  PSVNetMP Graph1;
  PSVNetMP Graph2;
  PSVNetMP Graph3;
  
  
  for (int iter = 0; iter < iter_count; iter++) {
    Graph1 = Graph.GetSubGraphMP(NodeVec1, EdgeVec1);
  }
  
  fprintf(outfile, "subgraph completed\n");
  PrintBenchmarks(outfile);
  printf("subgraph completed\n");
  
  for (int iter = 0; iter < iter_count; iter++) {
    Graph2 = Graph.GetSubGraphMP(NodeVec2, EdgeVec2);
  }
  
  fprintf(outfile, "subgraph completed\n");
  PrintBenchmarks(outfile);
  printf("subgraph completed\n");
  
  
  for (int iter = 0; iter < iter_count; iter++) {
    Graph3 = Graph.GetSubGraphMP(NodeVec3, EdgeVec3);
  }
  
  fprintf(outfile, "subgraph completed\n");
  PrintBenchmarks(outfile);
  printf("subgraph completed\n");
  
  printf("Original: Nodes = %d, edges = %d\n", Graph.GetNodes(), Graph.GetEdges());
  printf("Subgraph1: Nodes = %d, edges = %d\n", Graph1->GetNodes(), Graph1->GetEdges());
  printf("Subgraph2: Nodes = %d, edges = %d\n", Graph2->GetNodes(), Graph2->GetEdges());
  printf("Subgraph3: Nodes = %d, edges = %d\n", Graph3->GetNodes(), Graph3->GetEdges());
  printf("%d\n", Graph.GetSrcNId(2, 0));
  
  
  //Store bfs output
  /*
   printf("Storing Bfs output\n");
   FILE * outpr;
   outpr = fopen ("output_bfs.txt","w");
   for (int i = 0; i < NTblCnt; i++) {
   for (int j = 0; j < NTables[i]->GetNumRows(); j++) {
   fprintf(outpr, "%d \t%s \t%d\n", i, (NTables[i]->GetStrVal("NodeID", j)).CStr(), int(BfsLevelHV[i].GetDat(NodesHV[i].GetDat(NTables[i]->GetStrVal("NodeID", j)))));
   }
   }*/
  
  fprintf(outfile, "Output stored\n");
  PrintBenchmarks(outfile);
  printf("Output stored\n");
  
  
  fclose(outfile);
  //fclose(infile);
  //fclose(outpr);
}


