%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%  Benchmark problem: tension of the plate with hole;
%
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

clc; 
clear all;

E            = 1e5;  % Youngs modulus
nu           = 0.3;  % Poisson ratio
q0           = 10;  % uniform tension
a            = 1;    % radius of the hole
L            = 4;    % length of the plate
D0 = E/(1-nu^2)*[  1      nu        0;
                  nu     1          0;
                  0      0      (1-nu)/2  ];
dof = 2;             % degree of freedom
% define the geometry
plate_hole = geo_plate_with_hole(L,a);

% build iga mesh structure
mesh = build_iga_mesh( plate_hole );

% find the boundary nodes
epsilon=1e-5; 
bottomNodes = find(abs(mesh.coords(:,2))<epsilon);
topNodes    = find(abs(mesh.coords(:,2)-L)<epsilon);
leftNodes   = find(abs(mesh.coords(:,1))<epsilon);
rightNodes  = find(abs(mesh.coords(:,1)-L)<epsilon);
forceNodes  = (mesh.nCptsV-1)*mesh.nCptsU+1 : mesh.nCpts;
forceNodes = forceNodes';

dbc = [];    % dbc = [node index, direction, prescribed displacement]
dbc = [dbc; leftNodes,     ones(size(leftNodes)),   zeros(size(leftNodes))];
dbc = [dbc; bottomNodes,     2*ones(size(bottomNodes)),   zeros(size(bottomNodes))];
scatdbc = [];
scattbc = [];
if ~isempty(dbc)
    scatdbc = dof * (dbc(:,1)-1) + dbc(:,2);   % scatter dbc
end

[forceElems, forceKnots] = build_knot_connectivity( mesh.uKnots );
forceElems = forceElems + (mesh.nCptsV-1)*mesh.nCptsU;

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
        [R,ders] = nurbs_derivatives( gauPts, plate_hole, mesh );
        jmatrix = ders*elCpts(:,1:2); 
        j2      = det(jmatrix);
        dR_dx   =  jmatrix \ ders;     
        B = zeros(3,nnElem); 
        B(1,1:2:nnElem) = dR_dx(1,:);  
        B(2,2:2:nnElem) = dR_dx(2,:);   
        B(3,1:2:nnElem) = dR_dx(2,:);   
        B(3,2:2:nnElem) = dR_dx(1,:);   
        fac = j1 *j2 * wt;      
        K(sctrB,sctrB) = K(sctrB,sctrB) + B' * D0 * B * fac;
    end 
end

% boundary conditions imposition
gp_x = mesh.p+1;
[gp, wgt] = gauss_quadrature(gp_x);
for e = 1:size(forceElems,1)
    sctr   = forceElems(e,:);  
    elDoma = forceKnots(e,:);
    elCpts = mesh.coords(sctr,:);
    nn     = numel(sctr);
    nnElem = nn*dof;                 
    sctrB  = zeros(1, nnElem); 
    for i = 1:dof
        sctrB(i:dof:nnElem) = dof*(sctr-1) + i;  % displacement in i-th direction
    end   
    for ipt = 1:size(gp,1) 
        pt     = gp(ipt);                                     % quadrature point
        wt     = wgt(ipt);                                    % quadrature weight
        gauPts = parameter_gauss_mapping( elDoma, pt );       % gauss integration mapping
        j1     = jacobian_gauss_mapping( elDoma );            % jacobian value for gauss mapping      
        [N, dNdxi] = nrbNurbs1DBasisDerivs(gauPts, mesh.p, mesh.uKnots, mesh.coords(forceNodes,4));
        j2     = dNdxi*elCpts(:,1:2);            
        j2     = norm(j2);           
        
        x      = N*elCpts(:,1:2);  % compute exact tractions
        [str,disp]    = exact_solution_plate_hole(x,a,E,nu);
        
        if e <= mesh.nElemU/2    % top edge
            tx  = q0*str(3);
            ty  = q0*str(2);
        else                     % right edge
            tx  = q0*str(1);
            ty  = q0*str(3);
        end
        R = zeros(2,nnElem);
        R(1,1:2:nnElem) = N;
        R(2,2:2:nnElem) = N;
        fac = j1 * j2 * wt;
        F(sctrB) = F(sctrB) + R'* [tx; ty] * fac;
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
num1 = 50;  % if you want to obtain a more smooth results, please provide large numbers for num1 and num2    
num2 = 50;
polygon = build_visual_mesh_suf( num1, num2 );        % build visualized mesh
kntcrv = build_visual_knotcurve_suf( mesh.uKnots, mesh.vKnots, num1+1 ); % build visualized knot curves
numpts = (num1+1)*(num2+1) + size(kntcrv.linpts,1);
vmesh.vertices = zeros(numpts,2);
vmesh.displacement = zeros(numpts,2);
vmesh.exactdisplacement = zeros(numpts,2);
vmesh.stress = zeros(numpts,3);
vmesh.exactstress = zeros(numpts,3);
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
    [R,dRdparm] = nurbs_derivatives( [xi, eta], plate_hole, mesh );
    jmatrix = dRdparm*exyz(:,1:2); 
    dR_dx =  jmatrix \ dRdparm;   
    B(1,1:2:nnElem) = dR_dx(1,:);  
    B(2,2:2:nnElem) = dR_dx(2,:);   
    B(3,1:2:nnElem) = dR_dx(2,:);   
    B(3,2:2:nnElem) = dR_dx(1,:); 
    strain = B * u(sctrB);
    stress = D0 * strain;
    edsp   = u(sctrB);
    edsp   = reshape(edsp, dof, nn);
    vmesh.displacement(i,:) = R * edsp';
    vmesh.vertices(i,1:2) = R*exyz(:,1:2);
    vmesh.stress(i,:) = stress';
    [str,disp]    = exact_solution_plate_hole(vmesh.vertices(i,1:2),a,E,nu);
    vmesh.exactstress(i,:) = q0*str;
    vmesh.exactdisplacement(i,:) = q0*disp;
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
    nn = numel(sctr);           % number of control points for each element
    nnElem = nn*dof;            % dof for each element
    sctrB = zeros(1, nnElem); 
    for k = 1:dof
        sctrB(k:dof:nnElem) = dof*(sctr-1) + k;  % displacement in i-th direction
    end
    [R,dRdparm] = nurbs_derivatives( [xi, eta], plate_hole, mesh );
    jmatrix = dRdparm*exyz(:,1:2); 
    dR_dx =  jmatrix \ dRdparm;   
    B(1,1:2:nnElem) = dR_dx(1,:);  
    B(2,2:2:nnElem) = dR_dx(2,:);   
    B(3,1:2:nnElem) = dR_dx(2,:);   
    B(3,2:2:nnElem) = dR_dx(1,:);  
    strain = B * u(sctrB);
    stress = D0 * strain;
    edsp   = u(sctrB);
    edsp   = reshape(edsp, dof, nn);
    vmesh.displacement(count,:) = R * edsp';
    vmesh.vertices(count,1:2) = R*exyz(:,1:2); 
    vmesh.stress(count,:) = stress';    
    [str,disp]    = exact_solution_plate_hole(vmesh.vertices(count,1:2),a,E,nu);
    vmesh.exactstress(count,:) = q0*str;
    vmesh.exactdisplacement(count,:) = q0*disp;
end

% plot approximate solutions
figure;
face = polygon.trimesh;
maxnum = max(max(face));
vertices = vmesh.vertices(1:maxnum,:);
displacement = vmesh.displacement(1:maxnum,:);
stress = vmesh.stress(1:maxnum,:);
title( 'Approximate Solutions', 'FontSize', 14', 'FontWeight', 'Bold') ;
axis off;
linmesh = kntcrv.linmesh;
colX = linspace( 0.1, 0.65, 3 );
rowY = [0.55, 0.15];
for K = 1:6
    rowId = ceil(K/3);
    colId = K - (rowId-1)*3;
    axes( 'Position', [colX(colId), rowY(rowId), 0.3, 0.3] ) ;
    p = patch('Faces',face, 'Vertices', vertices);
    if K == 1,  cdata = displacement(1:maxnum,1); title('U_x');
    elseif K == 2,  cdata = displacement(1:maxnum,2);  title('U_y');
    elseif K == 3,  cdata = sqrt( displacement(1:maxnum,1).^2 + displacement(1:maxnum,2).^2  ); title('U Magnitude');
    elseif K == 4,  cdata = stress(1:maxnum,1);  title('\sigma_{xx}');  
    elseif K == 5,  cdata = stress(1:maxnum,2);  title('\sigma_{yy}');  
    elseif K == 6,  cdata = stress(1:maxnum,3);  title('\sigma_{xy}');   
    end
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
end

% plot exact solutions
figure;
face = polygon.trimesh;
maxnum = max(max(face));
vertices = vmesh.vertices(1:maxnum,:);
displacement = vmesh.exactdisplacement(1:maxnum,:);
stress = vmesh.exactstress(1:maxnum,:);
title( 'Exact Solutions', 'FontSize', 14', 'FontWeight', 'Bold') ;
axis off;
linmesh = kntcrv.linmesh;
colX = linspace( 0.1, 0.65, 3 );
rowY = [0.55, 0.15];
for K = 1:6
    rowId = ceil(K/3);
    colId = K - (rowId-1)*3;
    axes( 'Position', [colX(colId), rowY(rowId), 0.3, 0.3] ) ;
    p = patch('Faces',face, 'Vertices', vertices);
    if K == 1,  cdata = displacement(1:maxnum,1); title('U_x');
    elseif K == 2,  cdata = displacement(1:maxnum,2);  title('U_y');
    elseif K == 3,  cdata = sqrt( displacement(1:maxnum,1).^2 + displacement(1:maxnum,2).^2  ); title('U Magnitude');
    elseif K == 4,  cdata = stress(1:maxnum,1);  title('\sigma_{xx}');  
    elseif K == 5,  cdata = stress(1:maxnum,2);  title('\sigma_{yy}');  
    elseif K == 6,  cdata = stress(1:maxnum,3);  title('\sigma_{xy}');   
    end
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
end

% plot the absolute error
figure;
face = polygon.trimesh;
maxnum = max(max(face));
vertices = vmesh.vertices(1:maxnum,:);
displacement = vmesh.displacement(1:maxnum,:) - vmesh.exactdisplacement(1:maxnum,:);
stress = vmesh.stress(1:maxnum,:) - vmesh.exactstress(1:maxnum,:);
title( 'Absolute Error', 'FontSize', 14', 'FontWeight', 'Bold') ;
axis off;
linmesh = kntcrv.linmesh;
colX = linspace( 0.1, 0.65, 3 );
rowY = [0.55, 0.15];
for K = 1:6
    rowId = ceil(K/3);
    colId = K - (rowId-1)*3;
    axes( 'Position', [colX(colId), rowY(rowId), 0.3, 0.3] ) ;
    p = patch('Faces',face, 'Vertices', vertices);
    if K == 1,  cdata = displacement(1:maxnum,1); title('U_x');
    elseif K == 2,  cdata = displacement(1:maxnum,2);  title('U_y');
    elseif K == 3,  cdata = sqrt( displacement(1:maxnum,1).^2 + displacement(1:maxnum,2).^2  ); title('U Magnitude');
    elseif K == 4,  cdata = stress(1:maxnum,1);  title('\sigma_{xx}');  
    elseif K == 5,  cdata = stress(1:maxnum,2);  title('\sigma_{yy}');  
    elseif K == 6,  cdata = stress(1:maxnum,3);  title('\sigma_{xy}');   
    end
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
end

