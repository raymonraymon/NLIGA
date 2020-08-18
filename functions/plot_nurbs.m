function plot_nurbs(nurbs, ctrl, param)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% plot nurbs curve, surface or volume
% 
% acceptable forms:
%       plotnurbs(nurbs, ctrl, param)
% 
% inputs:
%       nurbs: nurbs curve, surface or volume
%       ctrl:  a boolean value, if ctrl=1, plot control net of the nurbs
%       param: a boolean value, if param=1, plot parameter net of the nurbs
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

% default values
subnum = 100;
hold_flag = ishold;

if (iscell (nurbs.knots))                                       % plot nurbs volume or surface
    if (size (nurbs.knots,2) == 2)                              % plot a nurbs surface
        plot_nurbs_srf(nurbs, param);
    elseif (size (nurbs.knots,2) == 3)                          % plot the boundaries of a nurbs volume
        volume_boundaries = nrbextract (nurbs);
        
        hold on;   
        for i = 1:6
            plot_nurbs_srf (volume_boundaries(i), param);
        end        
    end
else                                                            % plot nurbs curve
    order = nurbs.order;
    p = nrbeval (nurbs, linspace (nurbs.knots(order), nurbs.knots(end-order+1), subnum));
    bc = getboundingbox(p);
    lc = sqrt((bc(1,1)-bc(1,2))^2+(bc(2,1)-bc(2,2))^2+(bc(3,1)-bc(3,2))^2);
    rc = lc/400;
    wc = rc*2;
    
    edgenum = subnum-1;
    edge = zeros(edgenum,6);
    for i = 2:size(p,2)
        edge(i-1,:) = reshape([p(:,i-1);p(:,i)],1,6);
    end
    
    hold on;
    for i = 1:size(edge,1)
        p1 = edge(i,1:3);
        p2 = edge(i,4:6);
        drawCylinder(p1,p2,rc,'b');
    end
    
    if param
        hold on;
        order = nurbs.order;
        p = nrbeval (nurbs, unique (nurbs.knots(order:end-order+1)));

        for i = 1:size(p,2)
            drawSphere(p(:,i), wc, 'g');
        end
    end
end

% calculate the bounding box of the nurbs entity
points = getctrlpoints(nurbs);
bb = getboundingbox(points);

length = sqrt((bb(1,1)-bb(1,2))^2+(bb(2,1)-bb(2,2))^2+(bb(3,1)-bb(3,2))^2);
r = length/130;
w = r/3;

% plot control net of the nurbs
if ctrl
    [center, edge] = getctrlnet(points);
    hold on;
    for i = 1:size(center,1)
        drawSphere(center(i,:), r, 'r');
    end

    for i = 1:size(edge,1)
        p1 = edge(i,1:3);
        p2 = edge(i,4:6);
        drawCylinder(p1,p2,w);
    end
end

if (~hold_flag)
    hold off;
end

% axis equal;
set(gca, 'XLim',[bb(1,1)-2*r bb(1,2)+2*r]);
set(gca, 'YLim',[bb(2,1)-2*r bb(2,2)+2*r]);
set(gca, 'ZLim',[bb(3,1)-2*r bb(3,2)+2*r]);
set(gcf,'color','w');
set(gcf, 'renderer', 'opengl')
light;
view([30 30]);


end