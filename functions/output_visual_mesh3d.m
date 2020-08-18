function vmesh = output_visual_mesh3d( fout, mat, geo, mesh, u, step, currentime )
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% output vmesh for 3d nonlinear isogeometric analysis
%  Input:
%    fout - figure handle
%    mat - material properties
%    geo - geometry
%    mesh - iga mesh structure
%    u - current displacement
%    step - iterative step index
%    currentime - current time
%  Output:
%    vmesh - visualized mesh structure
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

num1 = 15;
num2 = 5;
srfs = nrbextract(geo);   % extract six outer faces of the solid model
vmesh = cell(1,6);
numpts = zeros(1,6);
polygon = cell(1,6);
kntcrv = cell(1,6);

for i = 1:6    % iterate six faces
    srf = srfs(i);
    polygon{1,i} = build_visual_mesh_suf( num1, num2 );  % build visualized mesh
    kntcrv{1,i} = build_visual_knotcurve_suf( srf.knots{1}, srf.knots{2}, num1 ); % build visualized knot curves
    numpts(i) = (num1+1)*(num2+1) + size(kntcrv{1,i}.linpts,1);
    vmesh{1,i}.nodalpts = zeros(numpts(i),3);
    vmesh{1,i}.displacement = zeros(numpts(i),3);
    vmesh{1,i}.stress = zeros(numpts(i),6);
    vmesh{1,i}.strain = zeros(numpts(i),6);
    if i == 1     % face u=0 
        polygon{1,i}.tripts = [zeros(size(polygon{1,i}.tripts,1),1),polygon{1,i}.tripts];
        kntcrv{1,i}.linpts = [zeros(size(kntcrv{1,i}.linpts,1),1), kntcrv{1,i}.linpts];
    elseif i == 2  % face u=1
        polygon{1,i}.tripts = [ones(size(polygon{1,i}.tripts,1),1),polygon{1,i}.tripts];
        kntcrv{1,i}.linpts = [ones(size(kntcrv{1,i}.linpts,1),1), kntcrv{1,i}.linpts];
    elseif i == 3  % face v=0
        polygon{1,i}.tripts = [polygon{1,i}.tripts(:,1),zeros(size(polygon{1,i}.tripts,1),1),polygon{1,i}.tripts(:,2)];
        kntcrv{1,i}.linpts = [kntcrv{1,i}.linpts(:,1), zeros(size(kntcrv{1,i}.linpts,1),1), kntcrv{1,i}.linpts(:,2)];
    elseif i == 4  % face v=1
        polygon{1,i}.tripts = [polygon{1,i}.tripts(:,1),ones(size(polygon{1,i}.tripts,1),1),polygon{1,i}.tripts(:,2)];
        kntcrv{1,i}.linpts = [kntcrv{1,i}.linpts(:,1), ones(size(kntcrv{1,i}.linpts,1),1), kntcrv{1,i}.linpts(:,2)];
    elseif i == 5  % face w=0
        polygon{1,i}.tripts = [polygon{1,i}.tripts, zeros(size(polygon{1,i}.tripts,1),1)];
        kntcrv{1,i}.linpts = [ kntcrv{1,i}.linpts, zeros(size(kntcrv{1,i}.linpts,1),1) ];
    elseif i == 6  % face w=1
        polygon{1,i}.tripts = [polygon{1,i}.tripts, ones(size(polygon{1,i}.tripts,1),1)];
        kntcrv{1,i}.linpts = [ kntcrv{1,i}.linpts, ones(size(kntcrv{1,i}.linpts,1),1) ];
    end
    elem_index = find_point_span( mesh, polygon{1,i}.tripts );
    for j = 1:(num1+1)*(num2+1)
        xi = polygon{1,i}.tripts(j,1);   % u mesh point
        eta = polygon{1,i}.tripts(j,2);  % v mesh point
        zeta = polygon{1,i}.tripts(j,3);  % w mesh point
        e = elem_index(j);   % element number
        sctr = mesh.elNodeCnt(e,:);     % element control points index
        exyz = mesh.coords(sctr,:);  % element control points' coordinates
        nn = length(sctr);   % number of control points in the element
        nn3 = nn*3;          % degree of freedom of control points 
        sctrB = zeros(1, nn3); 
        sctrB(1:3:nn3) = 3*sctr - 2;   % displacement in x direction
        sctrB(2:3:nn3) = 3*sctr - 1;   % displacement in y direction
        sctrB(3:3:nn3) = 3*sctr;       % displacement in z direction

        edsp = u(sctrB);
        edsp = reshape(edsp, 3, nn);
%         R = nrbbasisfun ({xi, eta, zeta}, geo);
        [R,dRdparm] = nurbs_derivatives( [xi, eta, zeta], geo, mesh );      
        jmatrix = dRdparm*exyz(:,1:3); 
        dRdx =  jmatrix \ dRdparm;                     
        f = edsp * dRdx' + eye(3);
        [ pk2, ~ ] = constitutive_relation( mesh.dim, mat, f );
        ccy = pk2cauchy(pk2, f);
        vmesh{1,i}.displacement(j,:) = edsp * R';
        vmesh{1,i}.nodalpts(j,:) = R*exyz(:,1:3) + vmesh{1,i}.displacement(j,:);   
        vmesh{1,i}.stress(j,:) = ccy;
        strain = (f'*f-eye(3))/2;
        vmesh{1,i}.strain(j,:) = voigt(strain);               
    end
    count = (num1+1)*(num2+1);
    line_index = find_point_span( mesh, kntcrv{1,i}.linpts );
    kntcrv{1,i}.linmesh = kntcrv{1,i}.linmesh + (num1+1)*(num2+1);    
    for j = 1:size(kntcrv{1,i}.linpts,1)
        count = count+1;
        xi = kntcrv{1,i}.linpts(j,1);
        eta = kntcrv{1,i}.linpts(j,2);
        zeta = kntcrv{1,i}.linpts(j,3);
        e = line_index(j);   % element number
        sctr = mesh.elNodeCnt(e,:);     % element control points index
        exyz = mesh.coords(sctr,:);  % element control points' coordinates
        nn = length(sctr);   % number of control points in the element
        nn3 = nn*3;          % degree of freedom of control points 
        sctrB = zeros(1, nn3); 
        sctrB(1:3:nn3) = 3*sctr - 2;   % displacement in x direction
        sctrB(2:3:nn3) = 3*sctr - 1;   % displacement in y direction
        sctrB(3:3:nn3) = 3*sctr;       % displacement in z direction
        edsp = u(sctrB);
        edsp = reshape(edsp, 3, nn);
%         R = nrbbasisfun ({xi, eta, zeta}, geo);     
        [R,dRdparm] = nurbs_derivatives( [xi, eta, zeta],geo, mesh );      
        jmatrix = dRdparm*exyz(:,1:3); 
        dRdx =  jmatrix \ dRdparm;              
        f = edsp * dRdx' + eye(3);
        [ pk2, ~ ] = constitutive_relation( mesh.dim, mat, f );
        ccy = pk2cauchy(pk2, f);
        vmesh{1,i}.displacement(count,:) = edsp * R';
        vmesh{1,i}.nodalpts(count,:) = R*exyz(:,1:3) + vmesh{1,i}.displacement(count,:);   
        vmesh{1,i}.stress(count,:) = ccy;
        strain = (f'*f-eye(3))/2;
        vmesh{1,i}.strain(count,:) = voigt(strain);  
    end   
end

fprintf(fout,'STEP = %d, TIME = %e\n', step, currentime);
for i = 1:6
    for j = 1:numpts(i)
        fprintf(fout,'v %f %f %f\n', vmesh{1,i}.nodalpts(j,:));
        fprintf(fout,'d %f %f %f\n', vmesh{1,i}.displacement(j,:));
        fprintf(fout,'s %f %f %f %f %f %f\n', vmesh{1,i}.stress(j,:));
        fprintf(fout,'t %f %f %f %f %f %f\n', vmesh{1,i}.strain(j,:));
    end
end
count = 0;
for i = 1:6
    polygon{1,i}.trimesh = polygon{1,i}.trimesh + count;
    count = count + numpts(i);
    for j = 1:size(polygon{1,i}.trimesh,1)
        fprintf(fout,'f %d %d %d\n', polygon{1,i}.trimesh(j,:));
    end
end

count = 0;
for i = 1:6
    kntcrv{1,i}.linmesh = kntcrv{1,i}.linmesh + count;
    count = count + numpts(i);
    for j = 1:size(kntcrv{1,i}.linmesh,1)
        fprintf(fout, 'l ');
        for k = 1:size(kntcrv{1,i}.linmesh,2)
            fprintf(fout,'%d ', kntcrv{1,i}.linmesh(j,k));
        end
        fprintf(fout, '\n');
    end    
end


end

