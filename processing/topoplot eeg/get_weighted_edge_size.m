function Esize= get_weighted_edge_size(Edges,IsWeighted)

channel_pairs = length(Edges);
Esize = ones(1,channel_pairs)*2;

if IsWeighted
    %Esize=mapminmax(Edges,1,10);

    cmax = 5;
    cmin = 1;
    m = 5;

    index = fix((Edges-cmin)/(cmax-cmin)*m)+1;
    index(index<1) = 1;
    index(index>m) = m;
    Esize = index;
end