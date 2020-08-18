function [ tri_ns, tri_nd, line_nd ] = compute_nodal_variables( geo, mesh, cstrs, disp, tri, kntcrv )
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% compute the nodal stresses and displacements for visualization
% Input: 
%   geo     - geometry structure
%   mesh    - mesh struture
%   cstrs   - stress at the control points
%   disp    - displacements at the control points
%   tri     - visualized polygon structure
%   kntcrv  - visualized knot curve structure
% Output: 
%   tri_ns  - nodal stress on the vertices of each triangles in tri
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

elem_index = find_point_span( mesh, tri.tripts );
tri_ns = zeros(size(tri.tripts,1), size(cstrs,2));
tri_nd = zeros(size(tri.tripts,1), size(disp,2));
line_nd = zeros(size(kntcrv.linpts,1), size(disp,2));
for i = 1:size(tri.tripts,1)    % look over each nodal points
    xi = tri.tripts(i,1);       % u mesh point
    eta = tri.tripts(i,2);      % v mesh point
    e = elem_index(i);              % element number
    sctr = mesh.elNodeCnt(e,:);     % element control points index
    [R,~] = nurbs_derivatives( [xi, eta],geo, mesh );
    tri_ns(i,:) = R * cstrs(sctr,:);
    tri_nd(i,:) = R * disp(sctr,:);
end

line_index = find_point_span( mesh, kntcrv.linpts );
for i = 1:size(kntcrv.linpts,1)
    xi = kntcrv.linpts(i,1);
    eta = kntcrv.linpts(i,2);
    e = line_index(i);   % element number
    sctr = mesh.elNodeCnt(e,:);     % element control points index
    [R,~] = nurbs_derivatives( [xi, eta], geo, mesh );
    line_nd(i,:) = R * disp(sctr,:);
end



end

