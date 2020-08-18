function plot_nurbs_srf(nurbs, param)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% plot nurbs surface
% just for calling by function plotnurbs(nurbs, param)
% 
% acceptable forms:
%       plot_nurbs_srf(nurbs, param)
% 
% inputs:
%       nurbs: nurbs curve, surface or volume
%       param: a boolean value, if param=1, plot parameter net of the nurbs
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

% default values
subnum = 100;
hold_flag = ishold;

% polygon = build_visual_mesh_suf( subnum, subnum);       % build visualized mesh
% p = nrbeval (nurbs,polygon.tripts);                     % calculate the coordinates of points on the surface
% patch('Faces',polygon.trimesh, 'Vertices', p','FaceColor',[0,1,1],'LineStyle','none');

p = nrbeval (nurbs, {linspace(nurbs.knots{1}(nurbs.order(1)),nurbs.knots{1}(end-nurbs.order(1)+1),subnum) ...
                       linspace(nurbs.knots{2}(nurbs.order(2)),nurbs.knots{2}(end-nurbs.order(2)+1),subnum)});
h = surfl (squeeze(p(1,:,:)), squeeze(p(2,:,:)), squeeze(p(3,:,:)));   
set(h,'FaceColor',[0,1,1],'LineStyle','none');
% shading interp;      
% plot parameter net of the nurbs surface
if param
    hold on;
    kntcrv = build_visual_knotcurve_suf( nurbs.knots{1}, nurbs.knots{2}, subnum );        
    p = nrbeval (nurbs,kntcrv.linpts)';       % calculate the coordinates of points on the surface
    for j = 1:size(kntcrv.linmesh,1)
        vv = p(kntcrv.linmesh(j,:),:);
        plot3(vv(:,1), vv(:,2),vv(:,3), 'k-');
    end 
end

if (~hold_flag)
    hold off;
end

% axis equal;

end