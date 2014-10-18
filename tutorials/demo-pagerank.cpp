#include "Snap.h"
//using namespace TSnap;


//Function to read a table of nodes
PTable AddNodeTable() {
    TTableContext Context;
    Schema Nodes;
    Nodes.Add(TPair<TStr, TAttrType>("NodeID", atInt));
    char filename[50];
    printf("Adding Node Table\n");
    printf("Enter filename\n");
    scanf("%s", filename);
    TStr fname(filename);
    PTable T = TTable::LoadSS(Nodes, fname, Context);
    return T;
}

// Function to read in a table of edges
PTable AddEdgeTable() {
    TTableContext Context;
    Schema Edges;
    Edges.Add(TPair<TStr, TAttrType>("SrcID", atInt));
    Edges.Add(TPair<TStr, TAttrType>("DstID", atInt));
    char filename[50];
    printf("Adding Edge Table\n");
    printf("Enter filename\n");
    scanf("%s", filename);
    TStr fname(filename);
    PTable T = TTable::LoadSS(Edges, fname, Context);
    return T;
}

//one iteration of pagerank algorithm
void PRank_iteration(const TVec<PTable>& NTables,const TVec<PTable>& ETables,const TVec<TPair<TInt, TInt> > & Mapping, TVec<TFltV> & PRank, const TVec<TIntIntH > &NodeHash) {
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
        for (int NodeID = 0; NodeID < (NTables[SNType]->GetNumRows()).Val; NodeID++) {
            int curr_node = NTables[SNType]->GetIntVal("NodeID", NodeID);
            PRank[SNType][NodeID] = 1 - d;
            
            //check all the incoming edges to the node
            for (int EType = 0; EType < ETableCount; EType++) {
                
                //node type should match the destination node for edge type
                if (Mapping[EType].Val2 == SNType) {
                    for (int EID = 0; EID < (ETables[EType]->GetNumRows()).Val; EID++) {
                        if (curr_node == (ETables[EType]->GetIntVal("DstID", EID)).Val) {
                            //add to the pagerank score
                            int in_node = (ETables[EType]->GetIntVal("SrcID", EID)).Val;
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
    
    //Initialize pagerank value for each node as 1
    for (int i = 0; i < NTableCount; i++) {
        TFltV fltvec;
        int nrows = (NTables[i]->GetNumRows()).Val;
        for (int j = 0; j < nrows; j++) {
            fltvec.Add(1.0);
        }
        PRank.Add(fltvec);
    }
    
    //Create the hashtable from NodeID to index
    TVec<TIntIntH > NodeHash;
    for (int i = 0; i < NTableCount; i++) {
        TIntIntH hash;
        for (int j = 0; j < NTables[i]->GetNumRows().Val; j++) {
            hash.AddDat(NTables[i]->GetIntVal("NodeID", j), j);
        }
        NodeHash.Add(hash);
    }
    
    //Run pagerank for some iterations
    int numIterations = 100;
    for (int i = 0; i < numIterations; i++) {
        PRank_iteration(NTables, ETables, Mapping, PRank, NodeHash);
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
    TVec<TIntIntH > OutDeg;
    int NTableCount = NTables.Len();
    int ETableCount = ETables.Len();
    
    //Initialize the vector to 0
    for (int i = 0; i < NTableCount; i++) {
        TIntIntH inthash;
        int nrows = (NTables[i]->GetNumRows()).Val;
        for (int j = 0; j < nrows; j++) {
            inthash.AddDat((NTables[i]->GetIntVal("NodeID", j)).Val, 0);
        }
        OutDeg.Add(inthash);
    }
    printf("Initializiation of outdegrees to 0 done\n");
    // Calculate the outdegrees
    
    for (int i = 0; i < ETableCount; i++) {
        int SrcTbl = (Mapping[i].Val1).Val;
        int DstTbl = (Mapping[i].Val2).Val;
        printf("Src table = %d, dest table = %d\n", SrcTbl, DstTbl);
        int NumEdges = (ETables[i]->GetNumRows()).Val;
        printf("# of edges = %d\n", NumEdges);
        for (int j = 0; j < NumEdges; j++) {
            int key = (ETables[i]->GetIntVal("SrcID", j)).Val;
            int tempval = OutDeg[SrcTbl].GetDat(key);
            OutDeg[SrcTbl].DelKey(key);
            OutDeg[SrcTbl].AddDat(key, tempval+1);
            
            printf("nodes%d, id = %d, new_val = %d\n", SrcTbl, (ETables[i]->GetIntVal("SrcID", j)).Val, (OutDeg[SrcTbl].GetDat(key)).Val);
        }
    }
    printf("Outdegrees calculated\n");
    //Store the outdegrees in the tables.
    for (int i = 0; i < NTableCount; i++) {
        TIntV intvec;
        for (int j = 0; j < (NTables[i]->GetNumRows()).Val; j++) {
            intvec.Add(OutDeg[i].GetDat(NTables[i]->GetIntVal("NodeID", j)));
        }
        NTables[i]->StoreIntCol("OutDeg", intvec);
    }
    printf("Outdegrees stored\n");
}



int main(int argc, char* argv[]) {
    
    int NTableCount;
    int ETableCount;
    printf("Enter the number of node tables and edge tables\n");
    scanf("%d %d", &NTableCount, &ETableCount);
    
    TVec<PTable> NTables;
    TVec<PTable> ETables;
    TVec<PTable> OutNTables;
    TVec<TPair<TInt, TInt> > Mapping;
    
    //Read in tables
    for (int i = 0; i < NTableCount; i++) {
        NTables.Add(AddNodeTable());
    }
    for (int i = 0; i < ETableCount; i++) {
        ETables.Add(AddEdgeTable());
        printf("Enter the source and destination node table number (index starting from 0)\n");
        int srcid, dstid;
        scanf("%d %d", &srcid, &dstid);
        Mapping.Add(TPair<TInt, TInt>(srcid, dstid));
    }
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
            printf("NodeID = %d, Pagerank score = %f\n", (Result[i]->GetIntVal("NodeID", j)).Val, (Result[i]->GetFltVal("Pagerank", j)).Val);
        }
        printf("\n");
    }
}

