function [xs,ys,es,ecs] = viewSnapshots(data,options)
% processes data file from dlmread of minnet snapshots

% INPUT: 
% data - dlmread of dump from minnet
% options - optional parameters/etc

%%
opt=struct();

% by default, plot results
opt.plot=true;
% by default, plot nodes
opt.plotnodes=true;
% default pause time
opt.pause=.01;
% default step
opt.step = 1;

% if nonempty, view these specific frames
opt.frames = [];

% how far to go in data file
opt.cutoff = NaN;

% getmovie
opt.getMovie = false;
% framerate
opt.rate = 10;
% default movie name
opt.movieName = 'minnetMovie';
% default show axis
opt.axis = true;
% default range
opt.axlim = [];
% show node names
opt.labels = false;
% show title
opt.title = true;
% default time step
opt.dt = 1e-3;

% highlight pinned nodes/fixed nodes, whatever you put in this cell
% arraay
opt.pins = {};
% highlight boundary edges
opt.ebnd = {};
% edge line width
opt.lw = 2;
opt.ecolor=[0 0 0 1];
opt.ebndcolor=[.7 .3 0 1];


% copy over input options
if (exist('options','var'))
    opt = copyStruct(options,opt);
end

%%
if isnan(opt.cutoff)
    cutoff=size(data,1);
else
    cutoff=opt.cutoff;
end


%%
nn=0;
ne=0;
fm=0;
throwaway=0;
net = {};
xs=[];ys=[];
% xs=zeros(1000,10);
% ys=zeros(1000,10);

for lc=1:cutoff
%     [lc,nn,ne]
    % if counted nnode already
    if (nn==0)&&(ne==0)
        fm=fm+1;
        nn=data(lc,1);
        ncs(fm) = nn; continue
    elseif ne==0
        ne=data(lc,1);
        ecs(fm) = ne; continue
    end
    
    % parsing node data
    if nn>0
        nc=data(lc,1);
        xs(fm,nc) = data(lc,2);
        ys(fm,nc) = data(lc,3);
        nn=nn-1;
    % parsing edge data
    else
        es{ne}(fm,1) = data(lc,1);
        es{ne}(fm,2) = data(lc,2);
        ne=ne-1;
    end
end

% if either of the node or edge counter didn't reach zero, then the last
% batch of data is incomplete, don't use it when going through data with
% movie
if (nn~=0)||(ne~=0)
    throwaway = 1;
end

% set axis limits
% if opt.axis
    axlim=[min(xs,[],'all') max(xs,[],'all') min(ys,[],'all') max(ys,[],'all')];
% else
%     axlim=[min(xs,[],'all') max(xs,[],'all') min(ys,[],'all') max(ys,[],'all')];
% end

if ~isempty(opt.axlim)
    axlim = opt.axlim;
end


%% use this method to plot 
if (isempty(opt.frames) || opt.frames(1)==1)
    if (opt.plot)
        
        if (opt.plotnodes)
            scatter(xs(1,1:ncs(1)),ys(1,1:ncs(1)),20,'Filled');
            hold on
            scatter(0,0,50,'w','Filled')
        else
            scatter(0,0,50,'w','Filled')
            hold on
        end


        for ec = 1:ecs(1)
            if (~isempty(es{ec}(1,:)))
                plot([xs(1,es{ec}(1,1)) xs(1,es{ec}(1,2))], ...
                    [ys(1,es{ec}(1,1)) ys(1,es{ec}(1,2))], '-', 'LineWidth',opt.lw,'Color',opt.ecolor)
            end
        end

        % if user specified pins to highlight:
        if ~isempty(opt.pins)
            scatter(xs(1,opt.pins{1}),ys(1,opt.pins{1}),30,[.8 .2 .1],'LineWidth',opt.lw)
            scatter(0,0,60,'w','Filled')
        end

        if ~isempty(opt.ebnd)
            for ec = opt.ebnd{1}
                ec
                ecs(1)-ec+1
                plot([xs(1,es{ecs(1)-ec+1}(1,1)) xs(1,es{ecs(1)-ec+1}(1,2))], ...
                    [ys(1,es{ecs(1)-ec+1}(1,1)) ys(1,es{ecs(1)-ec+1}(1,2))], '-', 'LineWidth',opt.lw,'Color',opt.ebndcolor)
            end
        end
        
        % plot node labels
        if (opt.labels)
            for nc=1:ncs(1)
                xytxt = [xs(1,nc),ys(1,nc)]+[0 .001];
                text(xytxt(1),xytxt(2),num2str(nc),'FontSize',16)
            end
        end



        hold off
        axis('equal')
        if (opt.axis)
            axis(axlim)
        else
            axis(axlim)
            axis off
        end
        if (opt.title)
%             title(['t = ', num2str(1)])
            title(['t = ', num2str(0*opt.dt), 's'],'Interpreter','latex')
        end
        set(gcf,'color','w');
        set(gca,'DefaultTextInterpreter','latex','FontSize',30)
        %%
        pause(opt.pause)
        
        % grab frame if saving movie
        if opt.getMovie
            F(1) = getframe(gcf);
        end
        %%
    end
end

if isempty(opt.frames)
    viewframes = (1+opt.step):opt.step:(size(xs,1)-throwaway);
else
    viewframes = opt.frames;
end
for tc=viewframes
    if (opt.plot)
    % plot nodes
    if (opt.plotnodes)
        scatter(xs(tc,:),ys(tc,:),20,'Filled');
        hold on
        scatter(0,0,50,'w','Filled')
    else
        scatter(0,0,50,'w','Filled')
        hold on
    end

    for ec = 1:ecs(tc)
        if tc<=size(es{ec},1)
            if (~isempty(es{ec}(tc,:)))
% debugging
%                 if ~isempty(find(es{ec}(tc,:)==416))
%                     find(opt.ebnd{tc}==ecs(tc)-ec+1)
%                     es{ec}(tc,:)
%                 elseif ~isempty(find(es{ec}(tc,:)==1))
%                     find(opt.ebnd{tc}==ecs(tc)-ec+1)
%                     es{ec}(tc,:)
%                 elseif ~isempty(find(es{ec}(tc,:)==30))
%                     find(opt.ebnd{tc}==ecs(tc)-ec+1)
%                     es{ec}(tc,:)
%                 end

                    
                % plot edges connecting nodes
                plot([xs(tc,es{ec}(tc,1)) xs(tc,es{ec}(tc,2))], ...
                    [ys(tc,es{ec}(tc,1)) ys(tc,es{ec}(tc,2))], '-', 'LineWidth',opt.lw,'Color',opt.ecolor)
            end
        end
    end

    % if user specified pins to highlight:
    if ~isempty(opt.pins)
        scatter(xs(tc,opt.pins{tc}),ys(tc,opt.pins{tc}),30,[.8 .2 .1],'LineWidth',opt.lw)
%         scatter(xs(tc,opt.pins{tc}),ys(tc,opt.pins{tc}),40+tc,'LineWidth',2)
        scatter(0,0,60,'w','Filled')
    end

    if ~isempty(opt.ebnd)
        for ec = opt.ebnd{tc}
            plot([xs(tc,es{ecs(tc)-ec+1}(tc,1)) xs(tc,es{ecs(tc)-ec+1}(tc,2))], ...
                [ys(tc,es{ecs(tc)-ec+1}(tc,1)) ys(tc,es{ecs(tc)-ec+1}(tc,2))], '-', 'LineWidth',opt.lw,'Color',opt.ebndcolor)
        end
    end

    % plot node labels
    if (opt.labels)
        for nc=1:ncs(tc)
            xytxt = [xs(tc,nc),ys(tc,nc)]+[0 .001];
            text(xytxt(1),xytxt(2),num2str(nc),'FontSize',16)
        end
    end


    hold off
    axis equal
    if (opt.axis)
        axis(axlim)
    else
        axis(axlim)
        axis off
    end
    if (opt.title)
%         title(['t = ', num2str(tc)])
        title(['t = ', num2str((tc-1)*opt.dt), 's'],'Interpreter','latex')
    end
    set(gcf,'color','w');
    set(gca,'DefaultTextInterpreter','latex','FontSize',30)
    pause(opt.pause)
    
    % grab frame if saving movie
    if (opt.getMovie)
        if exist('F','var')==0
            F(1) = getframe(gcf);

        else
            F(end+1)=getframe(gcf);
        end
    end

    end
    
end




%%

if opt.getMovie
    writerObj = VideoWriter(['~/Desktop/',opt.movieName,'.avi']);
    writerObj.FrameRate = opt.rate;
    % set the image per second
    % open the video writer
    open(writerObj);
    % write the frames to the video
    for i=1:length(F)
        % convert the image to a frame
        frame = F(i) ;
        writeVideo(writerObj, frame);
    end
    % close the writer object
    close(writerObj);
end


end