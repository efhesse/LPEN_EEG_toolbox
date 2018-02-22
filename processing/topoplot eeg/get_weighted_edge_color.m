function Ecolor= get_weighted_edge_color(Edges,IsWeighted,colors,Ecaxis)

if  IsWeighted==1     
    Ecolor=Edges;
else
    Ecolor = ones(size(Edges,1),1);
end
if min(Ecolor)==max(Ecolor)
    Ecolor=repmat(colors,[length(Edges) 1]);
else
    if ~isempty(Ecaxis)
        cmin = Ecaxis(1);
        cmax = Ecaxis(2);

        m = size(colors,1);
        index = fix((Ecolor-cmin)/(cmax-cmin)*m)+1;
        index(index<1) = 1;
        index(index>m) = m;
        Ecolor = colors(index,:,:);
    else
        Ecolor=mapminmax(Ecolor',1,size(colors,1));
        Ecolor=colors(round(Ecolor),:,:);
    end
end