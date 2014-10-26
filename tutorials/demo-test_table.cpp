#include "Snap.h"

int main(int argc, char* []) {
    TTableContext Context;
    Schema Nodes;
    Nodes.Add(TPair<TStr, TAttrType>("Src", atStr));
    Nodes.Add(TPair<TStr, TAttrType>("Dst", atStr));
    PTable T = TTable::LoadSS(Nodes, "../../FlickrTables/photosCLEF/photo_tagger_edges_dupremoved.tsv", Context);
    T->SaveSS("test.tsv");
}

