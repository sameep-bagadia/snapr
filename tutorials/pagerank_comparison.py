import snap

G = snap.LoadEdgeList(snap.PNGraph, "input_udpg.txt")
PRankH = snap.TIntFltH()
snap.GetPageRank(G, PRankH)
print G.GetNodes()
print G.GetEdges()
for id in PRankH:
    print id, G.GetNI(id).GetOutDeg() , PRankH[id]