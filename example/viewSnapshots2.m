function viewSnapshots2(xs,ys,es,ecs,options)
% input is already processed data from non-plotting call

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
opt.movieName = ['minnetMovie',char(datetime("today"))];
% default show axis
opt.axis = true;
% default range
opt.axlim = [];
% show node names
opt.labels = false;
% show title
opt.title = true;
% move title?
opt.movetitle = false;
% default time step
opt.dt = 1e-3;
% show scale bar, user supply desired length
opt.scalebar=0;
% left text
opt.left='';
% top text
opt.top='';

% highlight pinned nodes/fixed nodes, whatever you put in this cell
% arraay
opt.pins = {};
% highlight boundary edges
opt.ebnd = {};
% edge line width
opt.lw = 2;
opt.ecolor=[0 0 0 1];
opt.ebndcolor=[.7 .3 0 1];

% plot box
opt.plotbox=[];

% copy over input options
if (exist('options','var'))
    opt = copyStruct(options,opt);
end

% set axis limits
axlim=[min(xs,[],'all') max(xs,[],'all') min(ys,[],'all') max(ys,[],'all')];
if ~isempty(opt.axlim)
    axlim = opt.axlim;
end

maxy=max(ys,[],'all');
minx=min(xs,[],'all');

% scale bar
if opt.scalebar
    mx=max(xs,[],'all');
    my=min(ys,[],'all');
    x1=mx-opt.scalebar;
    xscale=linspace(x1,mx);
    yscale=(my+opt.scalebar/5)*ones(size(xscale));
end

%% plot!
n = max(max(cellfun(@(x)size(x,1),es)));
temp=cellfun(@(x)[x; zeros(n-size(x,1),2)], es, 'uni', 0);
aa=cell2mat(temp);

if isempty(opt.frames)
    viewframes = 1:opt.step:(size(xs,1));
else
    viewframes = opt.frames;
end
%%
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
    edges=nonzeros(aa(tc,:));
    edgex=reshape(xs(tc,edges),2,[]);
    edgey=reshape(ys(tc,edges),2,[]);
    plot(edgex,edgey,'-','LineWidth',opt.lw,'Color',opt.ecolor)
    % if user specified pins to highlight:
    if ~isempty(opt.pins)
        xx=xs(tc,opt.pins{tc});
        yy=ys(tc,opt.pins{tc});
        inds1=find(xx==0);inds2=find(yy==0);inds=union(inds1,inds2);
        length(inds)
        length(xx)
        inds=setdiff(1:length(xx),inds);
        xx=xx(inds);yy=yy(inds);
        scatter(xx,yy,15,[.8 .2 .1],'LineWidth',opt.lw)
%         scatter(xs(tc,opt.pins{tc}),ys(tc,opt.pins{tc}),40+tc,'LineWidth',2)
%         scatter(0,0,60,'w','Filled')
    end
    
    if ~isempty(opt.ebnd)
        bndedges=nonzeros(aa(tc,opt.ebnd{tc}));
        bndedgex=reshape(xs(tc,bndedges),2,[]);
        bndedgey=reshape(ys(tc,bndedges),2,[]);
        plot(bndedgex,bndedgey,'-', 'LineWidth',opt.lw,'Color',opt.ebndcolor)
    end

    % plot box
    if ~isempty(opt.plotbox)
        xmin=opt.plotbox(1);xmax=opt.plotbox(2);
        ymin=opt.plotbox(3);ymax=opt.plotbox(4);
        plot([xmin xmax],[ymin ymin],'r-','Linewidth',2)
        plot([xmin xmax],[ymax ymax],'r-','Linewidth',2)
        plot([xmin xmin],[ymin ymax],'r-','Linewidth',2)
        plot([xmax xmax],[ymin ymax],'r-','Linewidth',2)
    end

    % plot node labels
    if (opt.labels)
        for nc=1:ncs(tc)
            xytxt = [xs(tc,nc),ys(tc,nc)]+[0 .001];
            text(xytxt(1),xytxt(2),num2str(nc),'FontSize',16)
        end
    end
    
    % scale bar
    if (opt.scalebar)
        plot(xscale,yscale,'k','linewidth',8)
        text(mean(xscale),mean(yscale)+opt.scalebar/5,[num2str(opt.scalebar),'$\mu$m'],...
            'Interpreter','latex','FontSize',40,'HorizontalAlignment','center')
    end

    % add text to top and left
    if ~isempty(opt.top)
        text(0,maxy+1,opt.top,'Interpreter','latex','fontsize',40,'HorizontalAlignment','center')
    end
    if ~isempty(opt.left)
        tleft=text(minx-1,0,opt.left,'Interpreter','latex','fontsize',40,'HorizontalAlignment','center');
        tleft.Rotation=90;
    end

    % axes
    hold off
    axis equal
    if (opt.axis)
        axis(axlim)
    else
        axis(axlim)
        axis off
    end

    % title
    if (opt.title)
        title(['t = ', num2str((tc-viewframes(1))*opt.dt), 's'],'Interpreter','latex')
        ax = gca;
        ax.TitleFontSizeMultiplier = 1.75;
        if opt.movetitle
            ax = gca;
            ax.TitleHorizontalAlignment = 'left';
        end
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

%% write movie to file
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