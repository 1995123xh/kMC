function out = outedge(G,nidx)
%to solve the problem of no 'outedges' function before R2018b
edges=[];
n=neighbors(G, nidx);
for nei = n
    edges=[edges, findedge(G, nidx,n)];
end
out=edges;
