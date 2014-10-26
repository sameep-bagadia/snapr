#include "Snap.h"
//using namespace TSnap;


//Function to read a table of nodes
PTable AddNodeTable(TTableContext & Context) {
    //TTableContext Context;
    Schema Nodes;
    Nodes.Add(TPair<TStr, TAttrType>("NodeID", atStr));
    char filename[50];
    int col_count;
    printf("Adding Node Table\n");
    printf("Enter filename and number of columns (>= 1) \n");
    scanf("%s %d", filename, &col_count);
    for (TInt i = 1; i < col_count; i++) {
        TStr col_name = "Attribute" + i.GetStr();
        Nodes.Add(TPair<TStr, TAttrType>(col_name, atStr));
    }
    TStr fname(filename);
    PTable T = TTable::LoadSS(Nodes, fname, Context);
    return T;
}

// Function to read in a table of edges
PTable AddEdgeTable(TTableContext & Context) {
    //TTableContext Context;
    Schema Edges;
    Edges.Add(TPair<TStr, TAttrType>("SrcID", atStr));
    Edges.Add(TPair<TStr, TAttrType>("DstID", atStr));
    char filename[200];
    int col_count;
    printf("Adding Edge Table\n");
    printf("Enter filename and number of columns (>= 2)\n");
    scanf("%s %d", filename, &col_count);
    for (TInt i = 1; i < col_count-1; i++) {
        TStr col_name = "Attribute" + i.GetStr();
        Edges.Add(TPair<TStr, TAttrType>(col_name, atStr));
    }
    TStr fname(filename);
    PTable T = TTable::LoadSS(Edges, fname, Context);
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
        TStr col_name = "Attribute" + i.GetStr();
        Edges.Add(TPair<TStr, TAttrType>(col_name, atStr));
    }
    TStr fname(filename);
    PTable T = TTable::LoadSS(Edges, fname, Context);
    return T;
}

//one iteration of pagerank algorithm
void PRank_iteration(const TVec<PTable>& NTables,const TVec<PTable>& ETables,const TVec<TPair<TInt, TInt> > & Mapping, TVec<TFltV> & PRank, const TVec<TStrIntH > &NodeHash, int total_nodes) {
    //Make a copy of page rank values
    float d = 0.85;
    TVec<TFltV> OldPR;
    for (int i = 0; i < PRank.Len(); i++) {
        TFltV fltvec;
        for (int j = 0; j < PRank[i].Len(); j++) {
            fltvec.Add(PRank[i][j]);
        }
        OldPR.Add(fltvec);
    }
    
    // update pagerank scores
    int NTableCount = NTables.Len();
    int ETableCount = ETables.Len();
    
    // for each node in the graph
    for (int SNType = 0; SNType < NTableCount; SNType++) {
        printf("Starting node type %d\n", SNType);
        int num_rows = (NTables[SNType]->GetNumRows()).Val;
        for (int NodeID = 0; NodeID < (NTables[SNType]->GetNumRows()).Val; NodeID++) {
            printf("SNType:%d, \tnode %d / %d\n", SNType, NodeID, num_rows);
            TStr curr_node = NTables[SNType]->GetStrVal("NodeID", NodeID);
            PRank[SNType][NodeID] = (1 - d) / float(total_nodes);
            
            //check all the incoming edges to the node
            for (int EType = 0; EType < ETableCount; EType++) {
                //node type should match the destination node for edge type
                if (Mapping[EType].Val2 == SNType) {
                    for (int EID = 0; EID < (ETables[EType]->GetNumRows()).Val; EID++) {
                        if (curr_node == ETables[EType]->GetStrVal("DstID", EID)) {
                            //add to the pagerank score
                            TStr in_node = ETables[EType]->GetStrVal("SrcID", EID);
                            int in_node_index = (NodeHash[Mapping[EType].Val1].GetDat(in_node)).Val;
                            PRank[SNType][NodeID] += (d * OldPR[Mapping[EType].Val1][in_node_index]) / (NTables[Mapping[EType].Val1]->GetIntVal("OutDeg", in_node_index)).Val;
                        }
                    }
                }
            }
        }
    }
}

void print_PRank(const TVec<TFltV> & PRank) {
    for (int i = 0; i < PRank.Len(); i++) {
        for (int j = 0; j < PRank[i].Len(); j++) {
            printf("%lf ", (PRank[i][j]).Val);
        }
        printf("\n");
    }
    printf("\n");
}

void Pagerank_MM(const TVec<PTable>& NTables,const TVec<PTable>& ETables,const TVec<TPair<TInt, TInt> > & Mapping, TVec<PTable>& Result) {
    TVec<TFltV> PRank;
    int NTableCount = NTables.Len();
    int ETableCount = ETables.Len();
    int total_nodes = 0;
    for (int i = 0; i < NTableCount; i++) {
        total_nodes += NTables[i]->GetNumRows();
    }
    float pr_init = 1.0 / float(total_nodes);
    
    //Initialize pagerank value for each node as 1/N
    for (int i = 0; i < NTableCount; i++) {
        TFltV fltvec;
        int nrows = (NTables[i]->GetNumRows()).Val;
        for (int j = 0; j < nrows; j++) {
            fltvec.Add(pr_init);
        }
        PRank.Add(fltvec);
    }
    
    //Create the hashtable from NodeID to index
    TVec<TStrIntH > NodeHash;
    for (int i = 0; i < NTableCount; i++) {
        TStrIntH hash;
        for (int j = 0; j < NTables[i]->GetNumRows().Val; j++) {
            hash.AddDat(NTables[i]->GetStrVal("NodeID", j), j);
        }
        NodeHash.Add(hash);
    }
    
    //Run pagerank for some iterations
    int numIterations = 10;
    for (int i = 0; i < numIterations; i++) {
        printf("Starting pagerank iteration number: %d\n", i);
        PRank_iteration(NTables, ETables, Mapping, PRank, NodeHash, total_nodes);
    }
    
    //Print the final pagerank values
    printf("Final Page rank values:\n");
    print_PRank(PRank);
    
    //Store the pagerank values in Result
    Result.Clr();
    for (int i = 0; i < NTableCount; i++) {
        PTable new_ptable = new TTable(*NTables[i]);
        new_ptable->StoreFltCol("Pagerank", PRank[i]);
        Result.Add(new_ptable);
    }
}


// Calculate the outdegrees of each node and stores as a new attribute in the node tables.
void Calc_Outdegrees(TVec<PTable>& NTables, const TVec<PTable>& ETables, const TVec<TPair<TInt, TInt> > & Mapping) {
    printf("Calculating out degrees\n");
    TVec<TStrIntH > OutDeg;
    int NTableCount = NTables.Len();
    int ETableCount = ETables.Len();
    
    //Initialize the vector to 0
    for (int i = 0; i < NTableCount; i++) {
        TStrIntH strhash;
        int nrows = (NTables[i]->GetNumRows()).Val;
        for (int j = 0; j < nrows; j++) {
            strhash.AddDat(NTables[i]->GetStrVal("NodeID", j), 0);
        }
        OutDeg.Add(strhash);
    }
    // Calculate the outdegrees
    
    for (int i = 0; i < ETableCount; i++) {
        int SrcTbl = (Mapping[i].Val1).Val;
        int DstTbl = (Mapping[i].Val2).Val;
        printf("Src table = %d, dest table = %d\n", SrcTbl, DstTbl);
        int NumEdges = (ETables[i]->GetNumRows()).Val;
        printf("# of edges = %d\n", NumEdges);
        for (int j = 0; j < NumEdges; j++) {
            TStr key = ETables[i]->GetStrVal("SrcID", j);
            int tempval = OutDeg[SrcTbl].GetDat(key);
            //OutDeg[SrcTbl].DelKey(key);
            OutDeg[SrcTbl].AddDat(key, tempval+1);
            
            //printf("nodes%d, id = %d, new_val = %d\n", SrcTbl, (ETables[i]->GetIntVal("SrcID", j)).Val, (OutDeg[SrcTbl].GetDat(key)).Val);
        }
    }
    printf("Outdegrees calculated\n");
    //Store the outdegrees in the tables.
    for (int i = 0; i < NTableCount; i++) {
        TIntV intvec;
        for (int j = 0; j < (NTables[i]->GetNumRows()).Val; j++) {
            intvec.Add(OutDeg[i].GetDat(NTables[i]->GetStrVal("NodeID", j)));
        }
        NTables[i]->StoreIntCol("OutDeg", intvec);
    }
    printf("Outdegrees stored\n");
}



int main(int argc, char* []) {
    
    int NTableCount;
    int ETableCount;
    printf("Enter the number of node tables and edge tables\n");
    scanf("%d %d", &NTableCount, &ETableCount);
    
    TVec<PTable> NTables;
    TVec<PTable> ETables;
    TVec<PTable> OutNTables;
    TVec<TPair<TInt, TInt> > Mapping;
    
    TTableContext Context;
    //Read in tables
    for (int i = 0; i < NTableCount; i++) {
        PTable P = AddNodeTable(Context);
        P->SaveSS("ntable_abc.txt");
        NTables.Add(P);
    }
    NTables[0]->SaveSS("ntabledef.txt");
    for (int i = 0; i < ETableCount; i++) {
        ETables.Add(AddEdgeTable(Context));
        ETables.Add(AddReverseEdgeTable(Context));
        printf("Enter the source and destination node table number (index starting from 0)\n");
        int srcid, dstid;
        scanf("%d %d", &srcid, &dstid);
        Mapping.Add(TPair<TInt, TInt>(srcid, dstid));
        Mapping.Add(TPair<TInt, TInt>(dstid, srcid));
    }
    ETableCount *= 2;
    
    /* TEST
     
     
     
    NTables[0]->SaveSS("ntable0.txt");
    
    printf("check: %s\n", NTables[0]->GetStrVal("NodeID", 2).CStr());
    */
    
    ///////////////
    
    TVec<PTable> Result;
    
    //Caculate outdegrees of each node
    Calc_Outdegrees(NTables, ETables, Mapping);
    
    printf("starting pagerank\n");
    //Call Pagerank
    Pagerank_MM(NTables, ETables, Mapping, Result);
    
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

