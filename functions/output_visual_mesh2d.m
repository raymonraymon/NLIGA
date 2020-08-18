function vmesh = output_visual_mesh2d( fout, mat, geo, mesh, u, step, currentime )
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Main frame for nonlinear isogeometric analysis
%  Input:
%    filename   - filename of visualized mesh
%    mat        - material properties
%    geo        - geometry
%    mesh       - iga mesh structure
%    u          - current displacement
%    step       - iterative step index
%    currentime - current time
%  Output:
%    vmesh      - visualized mesh structure
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

num1 = 100;        
num2 = 10;
polygon = build_visual_mesh_suf( num1, num2 );        % build visualized mesh
kntcrv = build_visual_knotcurve_suf( mesh.uKnots, mesh.vKnots, num1+1 ); % build visualized knot curves

numpts = (num1+1)*(num2+1) + size(kntcrv.linpts,1);
vmesh.nodalpts = zeros(numpts,2);
vmesh.displacement = zeros(numpts,2);
vmesh.stress = zeros(numpts,3);
vmesh.strain = zeros(numpts,3);
elem_index = find_point_span( mesh, polygon.tripts );
for i=1:(num1+1)*(num2+1)
    xi = polygon.tripts(i,1);  % u mesh point
    eta = polygon.tripts(i,2);  % v mesh point
    e = elem_index(i);   % element number
    sctr = mesh.elNodeCnt(e,:);     % element control points index
    exyz = mesh.coords(sctr,:);  % element control points' coordinates
    nn = length(sctr);   % number of control points in the element
    nn2 = nn*2;          % degree of freedom of control points
    sctrB = zeros(1, nn2); 
    sctrB(1:2:nn2) = 2*sctr - 1;   % displacement in x direction
    sctrB(2:2:nn2) = 2*sctr;       % displacement in y direction
    edsp = u(sctrB);
    edsp = reshape(edsp, 2, nn);
%     R = nrbbasisfun ({xi, eta}, geo);
    [R,dRdparm] = nurbs_derivatives( [xi, eta],geo, mesh );
    jmatrix = dRdparm*exyz(:,1:2); 
    dRdx =  jmatrix \ dRdparm;              
    f = edsp * dRdx' + eye(2);   

    [ pk2, ~ ] = constitutive_relation( mesh.dim, mat, f );
    ccy = pk2cauchy(pk2, f);
    vmesh.displacement(i,:) = edsp * R';
    vmesh.nodalpts(i,:) = R*exyz(:,1:2) + vmesh.displacement(i,:);   
    vmesh.stress(i,:) = ccy;
    strain = (f'*f-eye(2))/2;
    vmesh.strain(i,:) = voigt(strain);
end

count = (num1+1)*(num2+1);
line_index = find_point_span( mesh, kntcrv.linpts );
kntcrv.linmesh = kntcrv.linmesh + (num1+1)*(num2+1);
for i = 1:size(kntcrv.linpts,1)
    count = count+1;
    xi = kntcrv.linpts(i,1);
    eta = kntcrv.linpts(i,2);
    e = line_index(i);   % element number
    sctr = mesh.elNodeCnt(e,:);     % element control points index
    exyz = mesh.coords(sctr,:);  % element control points' coordinates
    nn = length(sctr);   % number of control points in the element
    nn2 = nn*2;          % degree of freedom of control points
    sctrB = zeros(1, nn2); 
    sctrB(1:2:nn2) = 2*sctr - 1;   % displacement in x direction
    sctrB(2:2:nn2) = 2*sctr;       % displacement in y direction
    edsp = u(sctrB);
    edsp = reshape(edsp, 2, nn);
%     R = nrbbasisfun ({xi, eta}, geo);
    [R,dRdparm] = nurbs_derivatives( [xi, eta],geo, mesh );
    jmatrix = dRdparm*exyz(:,1:2); 
    dRdx =  jmatrix \ dRdparm;              
    f = edsp * dRdx' + eye(2);   
    [ pk2, ~ ] = constitutive_relation( mesh.dim, mat, f );
    cauchy = pk2cauchy( pk2, f );
    vmesh.displacement(count,:) = edsp * R';
    vmesh.nodalpts(count,:) = R*exyz(:,1:2) + vmesh.displacement(count,:);   
    vmesh.stress(count,:) = cauchy;
    strain = (f'*f-eye(2))/2;
    vmesh.strain(i,:) = voigt(strain);
end

fprintf(fout,'STEP = %d, TIME = %e\n', step, currentime);
for i = 1:numpts
    fprintf(fout,'v %f %f\n', vmesh.nodalpts(i,:));
    fprintf(fout,'d %f %f\n', vmesh.displacement(i,:));
    fprintf(fout,'s %f %f %f\n', vmesh.stress(i,:));
    fprintf(fout,'t %f %f %f\n', vmesh.strain(i,:));
end
for i = 1:size(polygon.trimesh,1)
    fprintf(fout,'f %d %d %d\n', polygon.trimesh(i,:));
end
for i = 1:size(kntcrv.linmesh,1)
    fprintf(fout, 'l ');
    for j = 1:size(kntcrv.linmesh,2)
        fprintf(fout,'%d ', kntcrv.linmesh(i,j));
    end
    fprintf(fout, '\n');
end


end

