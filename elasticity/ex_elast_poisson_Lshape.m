%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%  Benchmark problem: a Poisson problem on a L-shaped model;
%
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

clc; 
clear all;

L = 1;
P = [-1,-1];
geo = geo_L_shape( L, P );

% build iga mesh structure
mesh = build_iga_mesh( geo );

dof = 1;             % degree of freedom
% find the boundary nodes
bcodes1 = 1:mesh.nCptsU;
bcodes2    = mesh.nCptsU*(mesh.nCptsV-1)+1:mesh.nCpts;
bcodes3   = 1:mesh.nCptsU:mesh.nCpts;
bcodes4  = mesh.nCptsU:mesh.nCptsU:mesh.nCpts;
bcNodes = unique([bcodes1,bcodes2,bcodes3,bcodes4]');


dbc = [];    % dbc = [node index, direction, prescribed displacement]
dbc = [dbc; bcNodes,  ones(size(bcNodes)),   zeros(size(bcNodes))];

scatdbc = [];
scattbc = [];
if ~isempty(dbc)
    scatdbc = dof * (dbc(:,1)-1) + dbc(:,2);   % scatter dbc
end

% initialize stiffness and force matrices
nDofs = dof * mesh.nCpts;    % total dofs
K = sparse(nDofs,nDofs);     % stiffness matrix 
F = zeros(nDofs,1);          % external force matrix

% use gaussian integration rule
gp_x = mesh.p+1;           % number of integration points in x-direction
gp_y = mesh.q+1;           % number of integration points in y-direction
[gp, wgt] = gauss_quadrature(gp_x, gp_y);   % calculate integration points and its weights
for e = 1:mesh.nElems                 % loop over elements
    sctr   = mesh.elNodeCnt(e,:);     % element control points index
    elDoma = mesh.elDoma(e,:);        % element parametric domain
    elCpts = mesh.coords(sctr,:);     % coordinates of element control points
    nn     = numel(sctr);             % number of control points for each element
    nnElem = nn*dof;                  % dof for each element
    sctrB  = zeros(1, nnElem); 
    for i = 1:dof
        sctrB(i:dof:nnElem) = dof*(sctr-1) + i;  % displacement in i-th direction
    end   
    for ipt = 1:size(gp,1)        % loop over integration points
        pt      = gp(ipt,:);      % reference parametric coordinates for each integration point
        wt      = wgt(ipt);       % weigths for each integration point
        gauPts  = parameter_gauss_mapping( elDoma, pt );       % gauss integration mapping
        j1      = jacobian_gauss_mapping( elDoma );            % jacobian value for gauss mapping   
        [R,ders] = nurbs_derivatives( gauPts, geo, mesh );
        jmatrix = ders*elCpts(:,1:2); 
        j2      = det(jmatrix);
        dR_dx   =  jmatrix \ ders;     
        B = zeros(2,nnElem); 
        B(1,1:nnElem) = dR_dx(1,:);  
        B(2,1:nnElem) = dR_dx(2,:);    
        fac = j1 *j2 * wt;      
        K(sctrB,sctrB) = K(sctrB,sctrB) + B' * B * fac;
        
        x = R * elCpts(:,1:2);
        y = x(2);
        x = x(1);
        fx = 2*pi^2*sin(pi*x)*sin(pi*y);
        F(sctrB) = F(sctrB) + R'* fx * fac;
    end 
end

% solving equations
ndbc = size(dbc,1);   % number of displacement constrained nodes
K(scatdbc,:) = zeros(ndbc, nDofs);
K(scatdbc,scatdbc) = eye(ndbc);
% F(scatdbc,:) = 0;
F(scatdbc,:) = dbc(:,3);        
u = K\F;


% post-processing
num1 = 100;  % if you want to obtain a more smooth results, please provide large numbers for num1 and num2    
num2 = 100;
polygon = build_visual_mesh_suf( num1, num2 );        % build visualized mesh
kntcrv  = build_visual_knotcurve_suf( mesh.uKnots, mesh.vKnots, num1+1 ); % build visualized knot curves
numpts  = (num1+1)*(num2+1) + size(kntcrv.linpts,1);
vmesh.vertices = zeros(numpts,2);
vmesh.displacement = zeros(numpts,1);
vmesh.exactdisplacement = zeros(numpts,1);
vmesh.stress = zeros(numpts,2);
vmesh.exactstress = zeros(numpts,2);
elem_index = find_point_span( mesh, polygon.tripts );
for i=1:(num1+1)*(num2+1)
    xi = polygon.tripts(i,1);   % u mesh point
    eta = polygon.tripts(i,2);  % v mesh point
    e = elem_index(i);          % element number
    sctr = mesh.elNodeCnt(e,:); % element control points index
    exyz = mesh.coords(sctr,:); % element control points' coordinates
    nn = numel(sctr);           % number of control points for each element
    nnElem = nn*dof;            % dof for each element
    sctrB = zeros(1, nnElem); 
    for k = 1:dof
        sctrB(k:dof:nnElem) = dof*(sctr-1) + k;  % displacement in i-th direction
    end
    [R,dRdparm] = nurbs_derivatives( [xi, eta], geo, mesh );
    jmatrix = dRdparm*exyz(:,1:2); 
    dR_dx =  jmatrix \ dRdparm;   
    B = zeros(2,nnElem); 
    B(1,1:nnElem) = dR_dx(1,:);  
    B(2,1:nnElem) = dR_dx(2,:); 
    strain = B * u(sctrB);
    stress = strain;
    edsp   = u(sctrB);
    edsp   = reshape(edsp, dof, nn);
    vmesh.displacement(i) = R * edsp';
    vmesh.vertices(i,1:2) = R*exyz(:,1:2);
    vmesh.stress(i,:) = stress';
    x = vmesh.vertices(i,1);
    y = vmesh.vertices(i,2);
    ux = sin(x*pi)*sin(y*pi);
    vmesh.exactdisplacement(i) = ux;
end
count = (num1+1)*(num2+1);
line_index = find_point_span( mesh, kntcrv.linpts );
kntcrv.linmesh = kntcrv.linmesh + (num1+1)*(num2+1);
for i = 1:size(kntcrv.linpts,1)
    count = count+1;
    xi = kntcrv.linpts(i,1);
    eta = kntcrv.linpts(i,2);
    e = line_index(i);   % element number
    sctr = mesh.elNodeCnt(e,:);  % element control points index
    exyz = mesh.coords(sctr,:);  % element control points' coordinates
    nn = numel(sctr);            % number of control points for each element
    nnElem = nn*dof;             % dof for each element
    sctrB = zeros(1, nnElem); 
    for k = 1:dof
        sctrB(k:dof:nnElem) = dof*(sctr-1) + k;  % displacement in i-th direction
    end
    [R,dRdparm] = nurbs_derivatives( [xi, eta], geo, mesh );
    jmatrix = dRdparm*exyz(:,1:2); 
    dR_dx =  jmatrix \ dRdparm;   
    B = zeros(2,nnElem); 
    B(1,1:nnElem) = dR_dx(1,:);  
    B(2,1:nnElem) = dR_dx(2,:);  
    strain = B * u(sctrB);
    stress = strain;
    edsp   = u(sctrB);
    edsp   = reshape(edsp, dof, nn);
    vmesh.displacement(count) = R * edsp';
    vmesh.vertices(count,1:2) = R*exyz(:,1:2); 
    vmesh.stress(count,:) = stress';    
    
    x = vmesh.vertices(count,1);
    y = vmesh.vertices(count,2);
    ux = sin(x*pi)*sin(y*pi);
    vmesh.exactdisplacement(count) = ux;
end
linmesh = kntcrv.linmesh;
face = polygon.trimesh;
maxnum = max(max(face));
vertices = vmesh.vertices(1:maxnum,:);
displacement = vmesh.displacement(1:maxnum,:);
stress = vmesh.stress(1:maxnum,:);
figure;
p = patch('Faces',face, 'Vertices', vertices);
cdata = displacement(1:maxnum); title('Approximated Solution');
set(p,'FaceColor','interp','FaceVertexCData',cdata);
set(p,'EdgeColor','none');
hold on;
for j = 1:size(linmesh,1)
    vv = vmesh.vertices(linmesh(j,:),:);
    plot(vv(:,1), vv(:,2), 'k-');
end 
colorbar;
axis equal; 
axis off;
hold off;

figure;
displacement = vmesh.exactdisplacement(1:maxnum);
p = patch('Faces',face, 'Vertices', vertices);
cdata = displacement(1:maxnum,1); title('Exact Solution');
set(p,'FaceColor','interp','FaceVertexCData',cdata);
set(p,'EdgeColor','none');
hold on;
for j = 1:size(linmesh,1)
    vv = vmesh.vertices(linmesh(j,:),:);
    plot(vv(:,1), vv(:,2), 'k-');
end 
colorbar;
axis equal; 
axis off;
hold off;

figure;
displacement = vmesh.displacement(1:maxnum) - vmesh.exactdisplacement(1:maxnum);
p = patch('Faces',face, 'Vertices', vertices);
cdata = displacement(1:maxnum,1); title('Absolute Error');
set(p,'FaceColor','interp','FaceVertexCData',cdata);
set(p,'EdgeColor','none');
hold on;
for j = 1:size(linmesh,1)
    vv = vmesh.vertices(linmesh(j,:),:);
    plot(vv(:,1), vv(:,2), 'k-');
end 
colorbar;
axis equal; 
axis off;
hold off;



% compute the error
err1 = 0;
err2 = 0;
% use gaussian integration rule
gp_x = mesh.p+1;           % number of integration points in x-direction
gp_y = mesh.q+1;           % number of integration points in y-direction
[gp, wgt] = gauss_quadrature(gp_x, gp_y);   % calculate integration points and its weights
for e = 1:mesh.nElems                 % loop over elements
    sctr   = mesh.elNodeCnt(e,:);     % element control points index
    elDoma = mesh.elDoma(e,:);        % element parametric domain
    elCpts = mesh.coords(sctr,:);     % coordinates of element control points 
    for ipt = 1:size(gp,1)        % loop over integration points
        pt      = gp(ipt,:);      % reference parametric coordinates for each integration point
        wt      = wgt(ipt);       % weigths for each integration point
        gauPts  = parameter_gauss_mapping( elDoma, pt );       % gauss integration mapping
        [R,ders] = nurbs_derivatives( gauPts, geo, mesh );
        u_h = R * u(sctr);
        coords = R * elCpts(:,1:2);
        x = coords(1);
        y = coords(2);
        u_ex = sin(x*pi)*sin(y*pi);
        err1 = err1 + (u_h - u_ex)^2;
        err2 = err2 + u_ex^2;
    end 
end

err1 = sqrt(err1);
err2 = sqrt(err2);
err2 = err1/err2;

err = [err1 err2]
disp = [mesh.nElems  mesh.nCpts]


