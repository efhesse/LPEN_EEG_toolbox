function [position] = topoplot_connect(displayStruct,loc_file,EMARKERSIZE,c,IsUnDirected,Ecolor,Esize,Ealpha,Ecaxis)
% NOTE: The cartoon head is drawn using code in topoplot.m file from EEGLAB
% v6.01b (http://sccn.ucsd.edu/eeglab/). 
%
% Reference:
% Delorme A & Makeig S (2004) EEGLAB: an open source toolbox for analysis
% of single-trial EEG dynamics. Journal of Neuroscience Methods 134:9-21 
%
% Usage:
%
% >> topoplot_connect(ds, EEG.chanlocs);
%
% *ds* is the display strcture with the folowing fields:
%
% * *ds.chanPairs* (required) - N x 2 matrix, with N being the number of 
% connected channel pairs. For example, ds.chanPairs = [7, 12; 13 20]; 
% specifies two connected channel pairs (7, 12) and (13, 20).
% * *ds.connectStrength* (optional) - N x 1 matrix, a vector specifying
% connection strengths. If unspecified, then the connections will be
% rendered in a color at the center of the current colormap.
% * *ds.connectStrengthLimits* (optional) - 1 x 2 matrix specifying minimum
% and maximum values of connection strengths to display. If it is not 
% specified, then the minimum and maximum values from ds.connectStrength 
% are used.
%
% EEG.chanlocs is a structure specifying channel locations (or an locs
% filename)
%
% For comments and/or suggestions, please send me an email at
% praneeth@mit.edu
%
global fig
BACKCOLOR = [0.83 0.82 0.78];  % EEGLAB standard
rmax = 0.5;             % actual head radius - Don't change this!
%rmax = 0.64;             % actual head radius - Don't change this!
CIRCGRID   = 201;       % number of angles to use in drawing circles
AXHEADFAC = 0.9;        % head to axes scaling factor
% AXHEADFAC = 1.1;        % head to axes scaling factor
HEADCOLOR = [0 0 0];    % default head color (black)
%EMARKER = '.';          % mark electrode locations with small disks
EFSIZE = get(0,'DefaultAxesFontSize');
%EMARKER=ChanName;
%ECOLOR = [0 0 0];       % default electrode color = black
%EMARKERSIZE = [];       % default depends on number of electrodes, set in code
%EMARKERLINEWIDTH = 2;   % default edge linewidth for emarkers
HLINEWIDTH = 1.7;         % default linewidth for head, nose, ears
HEADRINGWIDTH    = .007;% width of the cartoon head ring

%
%%%%%%%%%%%%%%%%%%%% Read the channel location information %%%%%%%%%%%%%%%%%%%%%%%%
%
if ischar(loc_file)
    [tmpeloc, labels, Th Rd indices] = readlocs( loc_file,'filetype','loc');
elseif isstruct(loc_file) % a locs struct
    [tmpeloc, labels, Th Rd indices] = readlocs( loc_file );
    % Note: Th and Rd correspond to indices channels-with-coordinates only
else
    error('loc_file must be a EEG.locs struct or locs filename');
end
Th = pi/180*Th;                              % convert degrees to radians
plotchans = indices;
allchansind = 1:length(Th);

%%%%%%%%%%%%%%%%%%% remove infinite and NaN values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

[x,y]     = pol2cart(Th,Rd);  % transform electrode locations from polar to cartesian coordinates
plotchans = abs(plotchans);   % reverse indicated channel polarities
Rd        = Rd(plotchans);
x         = x(plotchans);
y         = y(plotchans);
plotrad = min(1.0,max(Rd)*1.02);            % default: just outside the outermost electrode location
 plotrad = max(plotrad,0.5);                 % default: plot out to the 0.5 head boundary
headrad = rmax;
pltchans = find(Rd <= plotrad); % plot channels inside plotting circle
x     = x(pltchans);
y     = y(pltchans);
squeezefac = rmax/plotrad;
x    = x*squeezefac;
y    = y*squeezefac;

%
%%%%%%%%%%%%%%%%%%%%%%% Draw blank head %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
cla
hold on
set(gca,'Xlim',[-rmax rmax]*AXHEADFAC,'Ylim',[-rmax rmax]*AXHEADFAC)

%
%%%%%%%%%%%%%%%%%%% Plot filled ring to mask jagged grid boundary %%%%%%%%%%%%%%%%%%%%%%%%%%%
%
hwidth = HEADRINGWIDTH;                   % width of head ring
hin  = squeezefac*headrad*(1- hwidth/2);  % inner head ring radius
circ = linspace(0,2*pi,CIRCGRID);
rx = sin(circ);
ry = cos(circ);

labels    = labels(plotchans); % remove labels for electrodes without locations
labels    = strvcat(labels); % make a label string matrix

%%%%%%%%%%%%%%%%%%%%%%%%% Plot cartoon head, ears, nose %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
headx = [[rx(:)' rx(1) ]*(hin+hwidth)  [rx(:)' rx(1)]*hin];
heady = [[ry(:)' ry(1) ]*(hin+hwidth)  [ry(:)' ry(1)]*hin];

patch(headx,heady,ones(size(headx)),HEADCOLOR,'edgecolor',HEADCOLOR); hold on

%%%%%%%%%%%%%%%%%%% Plot ears and nose %%%%%%%%%%%%%%%%%%%%%%%%%%%
%
base  = rmax-.0046;
basex = 0.18*rmax;                   % nose width
tip   = 1.15*rmax;
tiphw = .04*rmax;                    % nose tip half width
tipr  = .01*rmax;                    % nose tip rounding
q = .04; % ear lengthening
EarX  = [.497-.005  .510  .518  .5299 .5419  .54    .547   .532   .510   .489-.005]; % rmax = 0.5
EarY  = [q+.0555 q+.0775 q+.0783 q+.0746 q+.0555 -.0055 -.0932 -.1313 -.1384 -.1199];
sf    = headrad/plotrad;

plot3([basex;tiphw;0;-tiphw;-basex]*sf,[base;tip-tipr;tip;tip-tipr;base]*sf,...
    2*ones(size([basex;tiphw;0;-tiphw;-basex])),...
    'Color',HEADCOLOR,'LineWidth',HLINEWIDTH);                 % plot nose
plot3(EarX*sf,EarY*sf,2*ones(size(EarX)),'color',HEADCOLOR,'LineWidth',HLINEWIDTH)    % plot left ear
plot3(-EarX*sf,EarY*sf,2*ones(size(EarY)),'color',HEADCOLOR,'LineWidth',HLINEWIDTH)   % plot right ear


% %%%%%%%%%%%%%%%%%%% Show electrode information %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
plotax = gca;
axis square                                           % make plotax square
axis off

pos = get(gca,'position');
set(plotax,'position',pos);

xlm = get(gca,'xlim');
set(plotax,'xlim',xlm);

ylm = get(gca,'ylim');
set(plotax,'ylim',ylm);                               % copy position and axis limits again


% if isempty(EMARKERSIZE)
%     EMARKERSIZE = 10;
%     if length(y)>=32
%         EMARKERSIZE = 8;
%     elseif length(y)>=48
%         EMARKERSIZE = 6;
%     elseif length(y)>=64
%         EMARKERSIZE = 5;
%     elseif length(y)>=80
%         EMARKERSIZE = 4;
%     elseif length(y)>=100
%         EMARKERSIZE = 3;
%     elseif length(y)>=128
%         EMARKERSIZE = 3;
%     elseif length(y)>=160
%         EMARKERSIZE = 3;
%     end
% end
%
%%%%%%%%%%%%%%%%%%%%%%%% Mark electrode locations only %%%%%%%%%%%%%%%%%%%%%%%%%%
%
ELECTRODE_HEIGHT = -1;  % z value for plotting electrode information (above the surf)
scatter3(y,x,ones(size(x))*ELECTRODE_HEIGHT,EMARKERSIZE,c,'fill');
position=[y' ,x'];
%   for i = 1:size(labels,1)
%     hh(i) = text(double(y(i)+0.01),double(x(i)),...
%         ELECTRODE_HEIGHT,labels(i,:),'HorizontalAlignment','left',...
% 	'VerticalAlignment','middle','Color', ECOLOR,'userdata', num2str(allchansind(i)), ...
% 	'FontSize',EFSIZE, 'buttondownfcn', ...
% 	    ['tmpstr = get(gco, ''userdata'');'...
% 	     'set(gco, ''userdata'', get(gco, ''string''));' ...
% 	     'set(gco, ''string'', tmpstr); clear tmpstr;'] );
%   end
%

% praneeth - starts
numChanPairs = size(displayStruct.chanPairs, 1);
% cM = colormap;
% if ~isfield(displayStruct, 'connectStrength')
%     cmapPos = ceil(size(cM, 1)/1)*ones(size(displayStruct.chanPairs, 1), 1);
% else
%     if ~isfield(displayStruct, 'connectStrengthLimits')
%         displayStruct.connectStrengthLimits = [min(displayStruct.connectStrength), max(displayStruct.connectStrength)];
%     end
%     xp = displayStruct.connectStrengthLimits(1);
%     yp = displayStruct.connectStrengthLimits(2);
%     displayStruct.connectStrength(displayStruct.connectStrength < xp) = xp;
%     displayStruct.connectStrength(displayStruct.connectStrength > yp) = yp;
%     if xp == yp
%         cmapPos = ceil(size(cM, 1)/2)*ones(size(displayStruct.chanPairs, 1), 1);
%     else
%         cmapPos = round((displayStruct.connectStrength - xp)/(yp - xp)*(size(cM, 1) - 1) + 1);
%     end
% end
 %cM(cmapPos(kk), :)
%Esize=mapminmax(Esize,0.1,3);
for kk = 1:numChanPairs
    if IsUnDirected==1
    patchline(y(displayStruct.chanPairs(kk, :)), x(displayStruct.chanPairs(kk, :)), [ELECTRODE_HEIGHT, ELECTRODE_HEIGHT], 'LineWidth',Esize(kk), 'EdgeColor',Ecolor(kk,:),'facealpha',Ealpha(kk,1),'edgealpha',Ealpha(kk,2));
    else
%     arrSize=mapminmax(Esize,0.1,0.2);
arrSize=Esize/20;
    plot_arrow(y(displayStruct.chanPairs(kk,1)),x(displayStruct.chanPairs(kk,1)),y(displayStruct.chanPairs(kk,2)),x(displayStruct.chanPairs(kk,2)),'linewidth',Esize(kk),'headwidth',arrSize(kk)/1.5,'headheight',arrSize(kk) ,'EdgeColor',Ecolor(kk,:),'facecolor',Ecolor(kk,:),'facealpha',Ealpha(kk,1),'edgealpha',Ealpha(kk,2));
    end
end

% praneeth - ends
% set(gcf, 'color', BACKCOLOR);
lighting phong
%MY COMMENT
% %Add menue
% l=findall(fig,'type','uimenu');
% delete(l);
% l= uimenu(fig,'Label','window');
% uimenu(l,'Label','Title ...','Callback',{@title});
% uimenu(l,'Label','Load Files ...','Callback',@loadf);
% uimenu(l,'Label','BCT ...','Callback',@bct);
% uimenu(l,'Label','Exit','separator','on','callback',@exitp,'accelerator','z');
% mh1 = uimenu(fig,'Label','View');
% mh2 = uimenu(fig,'Label','Help');
colorbar
if ~isempty(Ecaxis)
    caxis(Ecaxis)
end
return;
end

function handles = plot_arrow( x1,y1,x2,y2,varargin )
%
% plot_arrow - plots an arrow to the current plot
%
% format:   handles = plot_arrow( x1,y1,x2,y2 [,options...] )
%
% input:    x1,y1   - starting point
%           x2,y2   - end point
%           options - come as pairs of "property","value" as defined for "line" and "patch"
%                     controls, see matlab help for listing of these properties.
%                     note that not all properties where added, one might add them at the end of this file.
%                     
%                     additional options are:
%                     'headwidth':  relative to complete arrow size, default value is 0.07
%                     'headheight': relative to complete arrow size, default value is 0.15
%                     (encoded are maximal values if pixels, for the case that the arrow is very long)
%
% output:   handles - handles of the graphical elements building the arrow
%
% Example:  plot_arrow( -1,-1,15,12,'linewidth',2,'color',[0.5 0.5 0.5],'facecolor',[0.5 0.5 0.5] );
%           plot_arrow( 0,0,5,4,'linewidth',2,'headwidth',0.25,'headheight',0.33 );
%           plot_arrow;   % will launch demo

% =============================================
% for debug - demo - can be erased
% =============================================
% if (nargin==0)
%     figure;
%     axis;
%     set( gca,'nextplot','add' );
%     for x = 0:0.3:2*pi
%         color = [rand rand rand];
%         h = plot_arrow( 1,1,50*rand*cos(x),50*rand*sin(x),...
%             'color',color,'facecolor',color,'edgecolor',color );
%         set( h,'linewidth',2 );
%     end
%     hold off;
%     return
% end
% =============================================
% end of for debug
% =============================================


% =============================================
% constants (can be edited)
% =============================================
alpha       = 0.15;   % head length
beta        = 0.07;   % head width
max_length  = 22;
max_width   = 10;

% =============================================
% check if head properties are given
% =============================================
% if ratio is always fixed, this section can be removed!
if ~isempty( varargin )
    for c = 1:floor(length(varargin)/2)
        try
            switch lower(varargin{c*2-1})
                % head properties - do nothing, since handled above already
            case 'headheight',alpha = max( min( varargin{c*2},1 ),0.01 );
            case 'headwidth', beta = max( min( varargin{c*2},1 ),0.01 );
            end
        catch
            fprintf( 'unrecognized property or value for: %s\n',varargin{c*2-1} );
        end
    end
end

% =============================================
% calculate the arrow head coordinates
% =============================================
den         = x2 - x1 + eps;                                % make sure no devision by zero occurs
teta        = atan( (y2-y1)/den ) + pi*(x2<x1) - pi/2;      % angle of arrow
cs          = cos(teta);                                    % rotation matrix
ss          = sin(teta);
R           = [cs -ss;ss cs];
line_length = sqrt( (y2-y1)^2 + (x2-x1)^2 );                % sizes
head_length = min( line_length*alpha,max_length );
head_width  = min( line_length*beta,max_length );
x0          = x2*cs + y2*ss;                                % build head coordinats
y0          = -x2*ss + y2*cs;
coords      = R*[x0 x0+head_width/2 x0-head_width/2; y0 y0-head_length y0-head_length];

% =============================================
% plot arrow  (= line + patch of a triangle)
% =============================================
h1          = patchline( [x1,x2],[y1,y2],[3,3]);
h2          = patch( coords(1,:),coords(2,:),[0 0 0] );
    
% =============================================
% return handles
% =============================================
handles = [h1 h2];

% =============================================
% check if styling is required 
% =============================================
% if no styling, this section can be removed!
if ~isempty( varargin )
    for c = 1:floor(length(varargin)/2)
        try
            switch lower(varargin{c*2-1})

             % only patch properties    
            case 'facelighting',set( h2,'FaceLighting',varargin{c*2} );
            case 'edgelighting',set( h2,'EdgeLighting',varargin{c*2} );

            % shared properties    
            case 'linestyle', set( handles,'LineStyle',varargin{c*2} );
            case 'linewidth', set( handles,'LineWidth',varargin{c*2} );
            case 'parent',    set( handles,'parent',varargin{c*2} );
            case 'edgealpha', set(handles,'EdgeAlpha',varargin{c*2});
            case 'facealpha', set(handles,'FaceAlpha',varargin{c*2});
            case 'facecolor', set( handles,'FaceColor',varargin{c*2} );
            case 'edgecolor', set( handles,'EdgeColor',varargin{c*2} );               


            % head properties - do nothing, since handled above already
            case 'headwidth',;
            case 'headheight',;
                
            end
        catch
            fprintf( 'unrecognized property or value for: %s\n',varargin{c*2-1} );
        end
    end
end
set(gcf,'color',[0.83 0.82 0.78]);
end
