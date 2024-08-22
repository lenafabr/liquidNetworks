% watch movie of output from example1
% then extract a network object and polygon information
%
% requires networktools code from : 
% https://github.com/lenafabr/networktools
%
% and spatialgraph2D code from FileExchange:
% https://www.mathworks.com/matlabcentral/fileexchange/73630-spatialgraph2d

addpath('../../networktools/')
addpath('../../spatialgraph2D/')

%% read and visualize from file with many snapshots
data = dlmread('example1.dump.out','',0,1);
% get node positions and edge connections
options=struct();
options.plot=false;
[xs,ys,es,ecs] = viewSnapshots(data,options);
%%
options=struct();
options.step=20;
options.dt=1e-3*10; % timestep times snapshot frequency
options.axis=false;
viewSnapshots2(xs,ys,es,ecs,options)

%% Extract network and polygons

% which frame:
tc=10;

% make network
NT = NetworkObj;
% set node positions and connections
NT.nodepos = [xs(tc,:)' ys(tc,:)'];
for ec = 1:ecs(tc)
    NT.edgenodes(ec,:) = [es{ec}(tc,1) es{ec}(tc,2)];
end

% setup network, reset edgelengths, remove double edges and only
% keeplargest component in case disconnected
NT.setupNetwork(true);
NT.removeDoubleEdges;
NT.keepLargestConnComp;

%% make graph, used for spatialgraph2D
G = NT.makeGraph;
% spatial graph 2d object:
obj = spatialgraph2D(G,NT.nodepos(:,1)',NT.nodepos(:,2)');
% obtain polyshape array of all constituent polygons:
pgon = polyshape(obj);
