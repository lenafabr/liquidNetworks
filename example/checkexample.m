% watch movie of output from example1

addpath('../../networktools/')

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