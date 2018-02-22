function topoplot_connectivity_EEG(mat,filename_to_save,channel_position_file,is_weighted,nodecolor,edgecolors,ecaxis,exts)

%load matrices
load(channel_position_file)
poslocs = chanlocs;

channel_nr = size(mat,1);
channel_pairs = channel_nr * channel_nr - channel_nr;

%create variables
%channel pair indices
[ds.chanPairs(:,1), ds.chanPairs(:,2)]=ind2sub(size(mat(:,:)),find(mat(:,:)));
Edges = mat(find(mat(:,:)));

%node size
Nsize = ones(1,channel_nr);
%for future use
%Nsize=mapminmax(Nsize',1,20);

%node color
%Ncolor = repmat([0 0 1],[channel_nr 1]);
Ncolor = repmat(nodecolor,[channel_nr 1]);
%for future use -> determine node color

%edge size
Esize= get_weighted_edge_size(Edges,is_weighted);

%edge color
color_weight = [];
if isstr(edgecolors)
%if is_weighted == 1 && isstr(edgecolors)
    colors=colormap(edgecolors);
    color_weight = 1;
else
    colors = edgecolors;
    color_weight = 0;
end

if isempty(ecaxis)
    ecaxis_vals = mat(find(mat(:) > 0));
    ecaxis_min = min(ecaxis_vals);
    ecaxis_max = max(ecaxis_vals);
    ecaxis = [ecaxis_min ecaxis_max];
end

Ecaxis = ecaxis;

Ecolor= get_weighted_edge_color(Edges,color_weight,colors,Ecaxis);
%Ecolor= get_weighted_edge_color(Edges,is_weighted,colors,Ecaxis);

%alpha for color of edges
%for future use
Ealpha = ones(channel_pairs,2);

is_undirected = 1;
%h = figure;
topoplot_connect(ds,poslocs,(Nsize+1).*20,Ncolor,is_undirected,Ecolor,Esize,Ealpha,Ecaxis);
colorbar
%savefig
for e = 1 : length(exts)
    saveas(h, [filename_to_save '.' exts{e}],exts{e})   
end    
   
