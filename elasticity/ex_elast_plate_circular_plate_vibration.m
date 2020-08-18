%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%  Free vibration analysis of a suqare plate based on Reissner-Mindlin assumption;
%  'nModes' is used to calculate the first nModes mode shapes
%  'bc' is used for boundary conditions
% 
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

clc; 
clear all;
% MATERIAL PREPERTIES
E      = 200*10^9;       % Youngs modulus
nu     = 0.3;            % Poisson ratio
h      = 0.2;            % thickness of the plate
rho    = 8000;           % density
nModes = 20;             % number of mode shapes

bc     = 2;              % 1-free, 2-simply supported,  3-clamped/fixed imposed
                         % on the outer boundary
                     
center = [0,0];          % center of the circle
r      = 1;              % radius of the circle
circle = geo_circle( center, r );   % define geometrical structure

% Build iga mesh structure
mesh = build_iga_mesh( circle );

% find four boundaries
left_nodes   = 1:mesh.nCptsU:mesh.nCpts;
bottom_nodes = 1:mesh.nCptsU; 
right_nodes  = mesh.nCptsU:mesh.nCptsU:mesh.nCpts;
top_nodes    = mesh.nCptsU*(mesh.nCptsV-1)+1:mesh.nCpts;
bc_nodes     = unique([left_nodes, bottom_nodes, right_nodes, top_nodes]);

fixed_dof = [];  %  constrained dofs
if bc == 2
    fixed_dof = bc_nodes*3-2;
elseif bc == 3
    fixed_dof = [bc_nodes*3-2, bc_nodes*3-1, bc_nodes*3 ];
end

% initialize stiffness and force matrices
dof = 3;
nDofs = dof * mesh.nCpts;    % total dofs
K = sparse(nDofs,nDofs);     % stiffness matrix 
F = zeros(nDofs,1);          % external force matrix
M = sparse(nDofs,nDofs);     % mass matrix

% Set the D matrix
lamda = 6/5;
D0 = E*h^3/(12*(1-nu*nu));
Db =D0*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];
Ds =E*h/(2*(1+nu)*lamda)*[1 0;0 1];
D  = {Db, zeros(3,2);zeros(2,3), Ds};
D  = cell2mat(D);

% use gaussian integration rule
gp_x = mesh.p+1;           % number of integration points in x-direction
gp_y = mesh.q+1;           % number of integration points in y-direction
[gp, wgt] = gauss_quadrature(gp_x, gp_y);   % calculate integration points and its weights
for e = 1:mesh.nElems                 % loop over elements
    sctr = mesh.elNodeCnt(e,:);       % element control points index
    elDoma = mesh.elDoma(e,:);        % element parametric domain
    elCpts = mesh.coords(sctr,:);     % coordinates of element control points
    nn = numel(sctr);                 % number of control points for each element
    nnElem = nn*dof;                  % dof for each element
    sctrB = zeros(1, nnElem); 
    for i = 1:dof
        sctrB(i:dof:nnElem) = dof*(sctr-1) + i;  % displacement in i-th direction
    end
    sctrw = 3*sctr-2;
    B = zeros(5, nnElem);
    N = zeros(3, nnElem);
  
    for ipt = 1:size(gp,1)   % loop over integration points
        pt = gp(ipt,:);      % reference parametric coordinates for each integration point
        wt = wgt(ipt);       % weigths for each integration point
        gauPts = parameter_gauss_mapping( elDoma, pt );   % gauss integration mapping
        j1 = jacobian_gauss_mapping( elDoma );            % jacobian value for gauss mapping   
        [R,ders] = nurbs_derivatives( gauPts, circle, mesh );
        jmatrix = ders*elCpts(:,1:2); 
        j2 = det(jmatrix);
        ders =  jmatrix \ ders;     
        B(1,3:3:nnElem) = ders(1,:);    
        B(2,2:3:nnElem) = -ders(2,:);   
        B(3,2:3:nnElem) = -ders(1,:);   
        B(3,3:3:nnElem) = ders(2,:);
        B(4,1:3:nnElem) = -ders(2,:);  
        B(4,2:3:nnElem) = R;
        B(5,1:3:nnElem) = ders(1,:);   
        B(5,3:3:nnElem) = R;
        
        N(1,1:3:nnElem) = R;
        N(2,2:3:nnElem) = R;
        N(3,3:3:nnElem) = R;
        m = rho*[h  0  0;  0  h^3/12  0;  0  0  h^3/12];
        fac = j1 *j2 * wt;      
        K(sctrB,sctrB) = K(sctrB,sctrB) + B' * D * B * fac;   
        M(sctrB,sctrB) = M(sctrB,sctrB) + N' * m * N * fac;
    end   
end

activedof = setdiff((1:nDofs), fixed_dof);
[modeShapes, freq] = eigs(K(activedof,activedof),M(activedof,activedof),nModes,0);
beta = diag(freq).^0.5*r^2*(rho*h/D0).^0.5;
[sortbeta, index_beta] = sort(beta);

% Post-processing
num1 = 30;  % if you want to obtain a more smooth results, please provide large numbers for num1 and num2    
num2 = 30;
polygon = build_visual_mesh_suf( num1, num2 );        % build visualized mesh
kntcrv = build_visual_knotcurve_suf( mesh.uKnots, mesh.vKnots, num1+1 ); % build visualized knot curves
numpts = (num1+1)*(num2+1) + size(kntcrv.linpts,1);
elem_index = find_point_span( mesh, polygon.tripts );
line_index = find_point_span( mesh, kntcrv.linpts );
kntcrv.linmesh = kntcrv.linmesh + (num1+1)*(num2+1);
figure(2);clf;
set( gcf, 'Color', 'White', 'Unit', 'Normalized', 'Position', [0.1,0.1,0.6,0.6] );
title( 'Mode Shapes of a Circular Plate', 'FontSize', 14', 'FontWeight', 'Bold') ;
axis off;
nCol = 5;  nRow = ceil( nModes / nCol );
rowH = 0.6 / nRow ;  colW = 0.8 / nCol;
colX = 0.06 + linspace( 0, 0.94, nCol+1);  
colX = colX(1:end-1) ;
rowY = 0.1 + linspace( 0.8, 0, nRow+1 ) ;  
rowY = rowY(2:end) ;
scale = 0.5/max(abs(modeShapes(:,1)));   % scale factor
for i = 1:nModes     % Loop over each mode shapes
    U = zeros(nDofs,1);
    U(activedof) = modeShapes(:,index_beta(i));
    vmesh.vertices = zeros(numpts,3);
    vmesh.displacement = zeros(numpts,3);
    for j = 1:(num1+1)*(num2+1)
        xi = polygon.tripts(j,1);   % u mesh point
        eta = polygon.tripts(j,2);  % v mesh point
        e = elem_index(j);          % element number
        sctr = mesh.elNodeCnt(e,:); % element control points index
        exyz = mesh.coords(sctr,:); % element control points' coordinates
        nn = numel(sctr);           % number of control points for each element
        nnElem = nn*dof;            % dof for each element
        sctrB = zeros(1, nnElem); 
        sctrB(1:3:nnElem) = 3*sctr-2;  % deflection   
        sctrB(2:3:nnElem) = 3*sctr-1;  % rotation x   
        sctrB(3:3:nnElem) = 3*sctr;    % rotation y
        [R,dRdparm] = nurbs_derivatives( [xi, eta], circle, mesh );
        edsp   = U(sctrB);
        edsp   = reshape(edsp, 3, nn);
        vmesh.displacement(j,:) = R * edsp';
        vmesh.vertices(j,:) = R * exyz(:,1:3);
        vmesh.vertices(j,3) = vmesh.displacement(j,1)*scale;
    end
    count = (num1+1)*(num2+1);
    for j = 1:size(kntcrv.linpts,1)
        count = count+1;
        xi = kntcrv.linpts(j,1);
        eta = kntcrv.linpts(j,2);
        e = line_index(j);   % element number
        sctr = mesh.elNodeCnt(e,:);     % element control points index
        exyz = mesh.coords(sctr,:);  % element control points' coordinates
        nn = numel(sctr);           % number of control points for each element
        nnElem = nn*dof;            % dof for each element
        sctrB = zeros(1, nnElem); 
        sctrB(1:3:nnElem) = 3*sctr-2;  % deflection   
        sctrB(2:3:nnElem) = 3*sctr-1;  % rotation x   
        sctrB(3:3:nnElem) = 3*sctr;    % rotation y
        [R,dRdparm] = nurbs_derivatives( [xi, eta], circle, mesh );
        edsp   = U(sctrB);
        edsp   = reshape(edsp, 3, nn);
        vmesh.displacement(count,:) = R * edsp';
        vmesh.vertices(count,1:2) = R * exyz(:,1:2);
        vmesh.vertices(count,3) = vmesh.displacement(count,1)*scale;   
    end
    rowId = ceil( i / nCol ) ;
    colId = i - (rowId - 1) * nCol ;
    axes( 'Position', [colX(colId), rowY(rowId), colW, rowH] ) ;
    face = polygon.trimesh;
    maxnum = max(max(face));
    vertices = vmesh.vertices(1:maxnum,:);
    p = patch('Faces',face, 'Vertices', vertices);
    box on;
    cdata = vmesh.displacement(1:maxnum,1); 
    set(p,'FaceColor','interp','FaceVertexCData',cdata);
    set(p,'EdgeColor','none');
    hold on;
    linmesh = kntcrv.linmesh;
    for j = 1:size(linmesh,1)
        vv = vmesh.vertices(linmesh(j,:),:);
        plot3(vv(:,1), vv(:,2), vv(:,3), 'k-');
    end 
    view(37,30);
    axis equal; 
    set(gca,'xtick',[],'ytick',[],'ztick',[]);
%     axis off;
    title( [sprintf( 'Modes %d', i), ',  \omega = ', sprintf('%3.3f',beta(index_beta(i)))]) ;  
    hold off;
end






