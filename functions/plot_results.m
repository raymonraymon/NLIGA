function plot_results( geo, vmesh, mesh, flag )
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Plot the results obtained in 2d elastoplasticity
% Input:
%   geo   - nurbs structure
%   vmesh - visualized mesh file
%   mesh  - mesh structure
%   flag  - for color render, 1-ux, 2-uy
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

clf;
nframe = length(vmesh.displacement);
subnum = 50;
tri = build_visual_mesh_suf( subnum, subnum);       % build visualized mesh
trivertex = nrbeval (geo,tri.tripts);               % calculate the coordinates of points on the surface
trivertex = trivertex';
face = tri.trimesh;
p = patch('Faces',face, 'Vertices', trivertex, 'FaceColor','g');         % plot the surface
set(p,'EdgeColor','none');
axis equal;
view(2);
kntcrv = build_visual_knotcurve_suf( geo.knots{1}, geo.knots{2}, subnum );  % build line mesh
linevertices = nrbeval (geo, kntcrv.linpts);  
linevertices = linevertices';
hold on; % plot knot curves
aa = zeros(size(kntcrv.linmesh,1),1);
for j = 1:size(kntcrv.linmesh,1)
    vv = linevertices(kntcrv.linmesh(j,:),:);
    aa(j) = plot(vv(:,1), vv(:,2), 'k-');
end  

coefs = geo.coefs;   % convert homogeneous coordinates to Euclidean space
for j = 1:geo.number(2)
    coefs(1,:,j) = coefs(1,:,j)./coefs(4,:,j);
    coefs(2,:,j) = coefs(2,:,j)./coefs(4,:,j);
end
for i = 1:nframe   % plot each frame  
    gstrs = vmesh.stress{i};
    disp =  vmesh.displacement{i};
    cstrs = patch_recovery( gstrs, mesh, geo);
    [ tri_ns, tri_nd, line_nd ] = compute_nodal_variables( geo, mesh, cstrs, disp, tri, kntcrv );
    if flag == 1       % render ux
        cdata = tri_nd(:,1);
    elseif flag == 2   % render uy
        cdata = tri_nd(:,2);
    elseif flag == 3   % render u magnitude
        cdata = zeros(size(tri_nd,1),1);
        for k=1:size(tri_nd,1)
            cdata(k) = norm(tri_nd(k,1:2));
        end
    elseif flag == 4   %  S11
        cdata = tri_ns(:,1);
    elseif flag == 5   %  S22
        cdata = tri_ns(:,2);
    elseif flag == 6   %  S12
        cdata = tri_ns(:,3);
    elseif flag == 7   % Mises
        cdata = von_mises( tri_ns );
    end
    
    if size(tri_nd,2) == 2, tri_nd(:,3) = 0; end
    tri_nd = tri_nd + trivertex;
    if size(line_nd,2) == 2, line_nd(:,3) = 0; end
    line_nd = line_nd + linevertices;   
  
    set(p, 'Faces',face, 'Vertices', tri_nd, 'FaceVertexCData',cdata);
    set(p,'FaceColor','interp','FaceVertexCData',cdata);
    for k = 1:size(kntcrv.linmesh,1)
        vv = line_nd(kntcrv.linmesh(k,:),:);
        set(aa(k),'xdata', vv(:,1), 'ydata', vv(:,2));
    end       
    title(['Step ',num2str(i)]);
    hcb = colorbar;
    if flag == 1,  title(hcb,'U_x');
    elseif flag == 2,  title(hcb,'U_y');
    elseif flag == 3,  title(hcb,'U, Magnitude');
    elseif flag == 4,  title(hcb,'S_{11}');
    elseif flag == 5,  title(hcb,'S_{22}');
    elseif flag == 6,  title(hcb,'S_{12}');
    elseif flag == 7,  title(hcb,'S, Mises');
    end
    drawnow;  
    pause(0.5);
    
end


end

