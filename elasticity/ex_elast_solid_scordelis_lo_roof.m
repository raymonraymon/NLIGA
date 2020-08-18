%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%  Elastic analysis of 3D Scordelis-Lo roof solid problem;
%  Converged solution: vertical displacement at midside of free edge =
%  0.3024
%
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

clc; 
clear all;

% build geometry
R = 25;
L = 50;
t = 0.25;
theta = 40;
geo = geo_scordelis_lo_roof_solid( L, R, t, theta );

% build iga mesh structure
mesh = build_iga_mesh( geo );

% material properties
E            = 4.32e8;  % Youngs modulus
nu           = 0.0;  % Poisson ratio
g            = [0,0,-90]';   % gravity loading

D0 = E/((1+nu)*(1-2*nu))*[  1-nu, nu, nu, 0, 0, 0;
                            nu, 1-nu, nu, 0, 0, 0;
                            nu, nu, 1-nu, 0, 0, 0;
                            0, 0, 0, (1-2*nu)/2, 0, 0;
                            0, 0, 0, 0, (1-2*nu)/2, 0;
                            0, 0, 0, 0, 0, (1-2*nu)/2;];
dof = 3;             % degree of freedom

% find the displacement boundary nodes
epsilon=1e-5; 
xzFaceNodes0 = find(abs(mesh.coords(:,2))<epsilon);       % face nodes on x-z plane with y = 0
xzFaceNodes1 = find(abs(mesh.coords(:,2)-L)<epsilon);     % face nodes on x-z plane with y = L/2

dbc = [];    % dbc = [node index, direction, prescribed displacement]
dbc = [dbc; xzFaceNodes0,     ones(size(xzFaceNodes0)),   zeros(size(xzFaceNodes0))];   
dbc = [dbc; xzFaceNodes0,     3*ones(size(xzFaceNodes0)),   zeros(size(xzFaceNodes0))];  
dbc = [dbc; xzFaceNodes1,     ones(size(xzFaceNodes1)),   zeros(size(xzFaceNodes1))];
dbc = [dbc; xzFaceNodes1,     3*ones(size(xzFaceNodes1)),   zeros(size(xzFaceNodes1))];
dbc = [dbc; 1,     2,   0]; % consrtrain the rigid body motion
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
gp_z = mesh.k+1;           % number of integration points in z-direction
[gp, wgt] = gauss_quadrature(gp_x, gp_y, gp_z);   % calculate integration points and its weights
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
        jmatrix = ders*elCpts(:,1:3); 
        j2      = det(jmatrix);
        dR_dx   =  jmatrix \ ders;     
        B = zeros(6,nnElem); 
        B(1,1:3:nnElem) = dR_dx(1,:);  
        B(2,2:3:nnElem) = dR_dx(2,:);
        B(3,3:3:nnElem) = dR_dx(3,:);   
        B(4,1:3:nnElem) = dR_dx(2,:);   
        B(4,2:3:nnElem) = dR_dx(1,:);   
        B(5,2:3:nnElem) = dR_dx(3,:);
        B(5,3:3:nnElem) = dR_dx(2,:); 
        B(6,1:3:nnElem) = dR_dx(3,:); 
        B(6,3:3:nnElem) = dR_dx(1,:); 
        fac = j1 *j2 * wt;      
        K(sctrB,sctrB) = K(sctrB,sctrB) + B' * D0 * B * fac;
    end 
end

% enforce force conditions
% find force boundary nodes and elements
forceNodes = (mesh.nCptsU*mesh.nCptsV*(mesh.nCptsW-1)+1:mesh.nCpts)';  % face nodes on x-y plane with z = 0
srfs = nrbextract(geo);   % extract six outer faces of the solid model
forceMesh = build_iga_mesh( srfs(6) ); 
forceMesh.elNodeCnt = forceMesh.elNodeCnt + mesh.nCptsU*mesh.nCptsV*(mesh.nCptsW-1);

gp_x = mesh.p+1;           % number of integration points in x-direction
gp_y = mesh.q+1;           % number of integration points in y-direction
[gp, wgt] = gauss_quadrature(gp_x, gp_y);   % calculate integration points and its weights
for e = 1:forceMesh.nElems
    sctr   = forceMesh.elNodeCnt(e,:);     % element control points index
    elDoma = forceMesh.elDoma(e,:);        % element parametric domain
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
        [R,ders] = nurbs_derivatives( gauPts, srfs(6), forceMesh );
        jmatrix = ders*elCpts(:,1:3); 
        j2      = norm(cross(jmatrix(1,:),jmatrix(2,:)));
        fac = j1 *j2 * wt;      
        N = zeros(3,nnElem);
        N(1,1:3:end) = R;
        N(2,2:3:end) = R;
        N(3,3:3:end) = R;
        F(sctrB) = F(sctrB) +  N' * g * fac;
    end 
end

% solving equations
ndbc = size(dbc,1);   % number of displacement constrained nodes
K(scatdbc,:) = zeros(ndbc, nDofs);
K(scatdbc,scatdbc) = eye(ndbc);
F(scatdbc,:) = dbc(:,3);     
u = K\F;
msize = get_nurbs_size( geo );
len = norm([msize(1), msize(3), msize(5)] - [msize(2), msize(4), msize(6)]);
factor = len/max(abs(u))/10;  % scale factor for better visualization

% post-processing
num1 = 30;
num2 = 30;
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
        nnElem = nn*dof;                  % dof for each element
        sctrB  = zeros(1, nnElem); 
        for k = 1:dof
            sctrB(k:dof:nnElem) = dof*(sctr-1) + k;  % displacement in i-th direction
        end   

        edsp = u(sctrB);
        edsp = reshape(edsp, dof, nn);
        [R,dRdparm] = nurbs_derivatives( [xi, eta, zeta], geo, mesh );      
        jmatrix = dRdparm*exyz(:,1:3); 
        dR_dx =  jmatrix \ dRdparm;                     
        B = zeros(6,nnElem); 
        B(1,1:3:nnElem) = dR_dx(1,:);  
        B(2,2:3:nnElem) = dR_dx(2,:);
        B(3,3:3:nnElem) = dR_dx(3,:);   
        B(4,1:3:nnElem) = dR_dx(2,:);   
        B(4,2:3:nnElem) = dR_dx(1,:);   
        B(5,2:3:nnElem) = dR_dx(3,:);
        B(5,3:3:nnElem) = dR_dx(2,:); 
        B(6,1:3:nnElem) = dR_dx(3,:); 
        B(6,3:3:nnElem) = dR_dx(1,:); 
        strain = B*u(sctrB);
        stress = D0 * strain;
        vmesh{1,i}.displacement(j,:) = edsp * R';
        vmesh{1,i}.nodalpts(j,:) = R*exyz(:,1:3) + vmesh{1,i}.displacement(j,:)*factor;   
        vmesh{1,i}.stress(j,:) = stress';         
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
        nnElem = nn*dof;                  % dof for each element
        sctrB  = zeros(1, nnElem); 
        for k = 1:dof
            sctrB(k:dof:nnElem) = dof*(sctr-1) + k;  % displacement in i-th direction
        end   
        edsp = u(sctrB);
        edsp = reshape(edsp, 3, nn);
        [R,dRdparm] = nurbs_derivatives( [xi, eta, zeta],geo, mesh );      
        jmatrix = dRdparm*exyz(:,1:3); 
        dR_dx =  jmatrix \ dRdparm;              
        B = zeros(6,nnElem); 
        B(1,1:3:nnElem) = dR_dx(1,:);  
        B(2,2:3:nnElem) = dR_dx(2,:);
        B(3,3:3:nnElem) = dR_dx(3,:);   
        B(4,1:3:nnElem) = dR_dx(2,:);   
        B(4,2:3:nnElem) = dR_dx(1,:);   
        B(5,2:3:nnElem) = dR_dx(3,:);
        B(5,3:3:nnElem) = dR_dx(2,:); 
        B(6,1:3:nnElem) = dR_dx(3,:); 
        B(6,3:3:nnElem) = dR_dx(1,:); 
        strain = B*u(sctrB);
        stress = D0 * strain;
        vmesh{1,i}.displacement(count,:) = edsp * R';
        vmesh{1,i}.nodalpts(count,:) = R*exyz(:,1:3);
%         vmesh{1,i}.nodalpts(count,:) = R*exyz(:,1:3) + vmesh{1,i}.displacement(count,:)*factor;   
        vmesh{1,i}.stress(count,:) = stress'; 
    end   
end

figure(2);
for i = 1:6
    face = polygon{1,i}.trimesh;
    maxnum = max(max(face));
    vertices = vmesh{1,i}.nodalpts(1:maxnum,:);
    displacement = vmesh{1,i}.displacement(1:maxnum,:);
    stress = vmesh{1,i}.stress(1:maxnum,:);
    linmesh = kntcrv{1,i}.linmesh;
    p = patch('Faces',face, 'Vertices', vertices);
    cdata = displacement(1:maxnum,3); title('Vertical Displacement');
    set(p,'FaceColor','interp','FaceVertexCData',cdata);
    set(p,'EdgeColor','none');
    hold on;
    for j = 1:size(linmesh,1)  % undeformed parametric curve
        vv = vmesh{1,i}.nodalpts(linmesh(j,:),:);
        plot3(vv(:,1), vv(:,2),vv(:,3), 'Color',[0.5,0.5,0.5]);
    end 
    for j = 1:size(linmesh,1)  % deformed parametric curve
        vv = vmesh{1,i}.nodalpts(linmesh(j,:),:) + vmesh{1,i}.displacement(linmesh(j,:),1:3)*factor;
        plot3(vv(:,1), vv(:,2),vv(:,3), 'k-');
    end 
end
axis equal;
view(3);

figure(3);
for i = 1:6
    face = polygon{1,i}.trimesh;
    maxnum = max(max(face));
    vertices = vmesh{1,i}.nodalpts(1:maxnum,:);
    displacement = vmesh{1,i}.displacement(1:maxnum,:);
    stress = vmesh{1,i}.stress(1:maxnum,:);
    linmesh = kntcrv{1,i}.linmesh;
    p = patch('Faces',face, 'Vertices', vertices);
    cdata = von_mises(stress(1:maxnum,:)); title('Mises Stress');
    set(p,'FaceColor','interp','FaceVertexCData',cdata);
    set(p,'EdgeColor','none');
    hold on;
    for j = 1:size(linmesh,1)  % undeformed parametric curve
        vv = vmesh{1,i}.nodalpts(linmesh(j,:),:);
        plot3(vv(:,1), vv(:,2),vv(:,3), 'Color',[0.5,0.5,0.5]);
    end 
    for j = 1:size(linmesh,1)  % deformed parametric curve
        vv = vmesh{1,i}.nodalpts(linmesh(j,:),:) + vmesh{1,i}.displacement(linmesh(j,:),1:3)*factor;
        plot3(vv(:,1), vv(:,2),vv(:,3), 'k-');
    end 
end
axis equal;
view(3);

% output the vertical displacement of the midside point
pt = [1e-15, 0.5, 0.5];
[R,ders] = nurbs_derivatives( pt, geo, mesh );
enum = find_point_span( mesh, pt );
sctr = mesh.elNodeCnt(enum,:); % element control points index
exyz = mesh.coords(sctr,:); % element control points' coordinates
vd = R * u(sctr*3)
nn = [mesh.nElems, mesh.nCpts]




