%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%  Static bending analysis of a circular Mindlin plate subjected to a uniform pressure;
%  'flag' is used to render the plate with obtained results like
%  displacements, bending moments and shear forces
%
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

clc; 
clear all;

E  = 2*10^8;         % Youngs modulus
nu = 0.3;            % Poisson ratio
h  = 0.01;           % thickness of the plate
q0 = -10;            % uniform pressure

flag  = 0;           % render the deformed plate with different fields' variables
                     % 0-all fields, 1-w, 2-rx, 3-ry, 4-Mx, 5-My, 6-Mxy, 7- Qy, 8-Qx;

bc = 1;              % clamped or simply supported Boundary Conditions, 
                     % 1 represents clamped and 2 represents simply supported    
                     
radius = 1;          % side length
center = [0,0];
circle = geo_circle( center, radius );

% Build iga mesh structure
mesh = build_iga_mesh( circle );

% find four boundaries
bottom_nodes = 1:mesh.nCptsU;  
left_nodes = 1:mesh.nCptsU:mesh.nCpts;
top_nodes = mesh.nCptsU*(mesh.nCptsV-1)+1:mesh.nCpts;
right_nodes = mesh.nCptsU:mesh.nCptsU:mesh.nCpts;
bc_nodes = (unique([bottom_nodes,left_nodes,top_nodes,right_nodes]))';
% four boundaries are clamped/fixed 
dof = 3;
dbc = [];    % dbc = [node index, direction, prescribed displacement]
if bc == 1   % clamped
    dbc = [dbc; bc_nodes,     ones(size(bc_nodes)),   zeros(size(bc_nodes))];
    dbc = [dbc; bc_nodes,   2*ones(size(bc_nodes)),   zeros(size(bc_nodes))];
    dbc = [dbc; bc_nodes,   3*ones(size(bc_nodes)),   zeros(size(bc_nodes))];
elseif bc == 2  % simply supported
    dbc = [dbc; bc_nodes,   ones(size(bc_nodes)),   zeros(size(bc_nodes))];
end
scatdbc = [];
scattbc = [];
if ~isempty(dbc)
    scatdbc = dof * (dbc(:,1)-1) + dbc(:,2);   % scatter dbc
end

% Set the D matrix
lamda = 6/5;
D0 = E*h^3/(12*(1-nu*nu));
Db =D0*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];
Ds =E*h/(2*(1+nu)*lamda)*[1 0;0 1];
D  = {Db, zeros(3,2);zeros(2,3), Ds};
D  = cell2mat(D);

% initialize stiffness and force matrices
nDofs = dof * mesh.nCpts;    % total dofs
K = sparse(nDofs,nDofs);   % stiffness matrix 
F = zeros(nDofs,1);        % external force matrix


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
        
        fac = j1 *j2 * wt;      
        K(sctrB,sctrB) = K(sctrB,sctrB) + B' * D * B * fac;  
        F(sctrw)  = F(sctrw) + q0 * R' * fac;
    end 
end

ndbc = size(dbc,1);   % number of displacement constrained nodes
K(scatdbc,:) = zeros(ndbc, nDofs);
K(scatdbc,scatdbc) = eye(ndbc);
F(scatdbc,:) = 0;
F(scatdbc,:) = dbc(:,3);        
u = K\F;

% post-processing
num1 = 50;  % if you want to obtain a more smooth results, please provide large numbers for num1 and num2    
num2 = 50;
polygon = build_visual_mesh_suf( num1, num2 );        % build visualized mesh
kntcrv = build_visual_knotcurve_suf( mesh.uKnots, mesh.vKnots, num1+1 ); % build visualized knot curves
numpts = (num1+1)*(num2+1) + size(kntcrv.linpts,1);
vmesh.vertices = zeros(numpts,3);
vmesh.displacement = zeros(numpts,3);
vmesh.stress = zeros(numpts,5);
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
    sctrB(1:3:nnElem) = 3*sctr-2;  % deflection   
    sctrB(2:3:nnElem) = 3*sctr-1;  % rotation x   
    sctrB(3:3:nnElem) = 3*sctr;    % rotation y
    [R,dRdparm] = nurbs_derivatives( [xi, eta], circle, mesh );
    jmatrix = dRdparm*exyz(:,1:2); 
    dRdx =  jmatrix \ dRdparm;   
    B(1,3:3:nnElem) = dRdx(1,:);    
    B(2,2:3:nnElem) = -dRdx(2,:);   
    B(3,2:3:nnElem) = -dRdx(1,:);   
    B(3,3:3:nnElem) = dRdx(2,:);
    B(4,1:3:nnElem) = -dRdx(2,:);  
    B(4,2:3:nnElem) = R;
    B(5,1:3:nnElem) = dRdx(1,:);   
    B(5,3:3:nnElem) = R;
    strain = B * u(sctrB);
    stress = D * strain;
    edsp   = u(sctrB);
    edsp   = reshape(edsp, 3, nn);
    vmesh.displacement(i,:) = R * edsp';
    vmesh.vertices(i,1:2) = R*exyz(:,1:2);
    vmesh.stress(i,:) = stress';
    if flag == 0        % all fields
        vmesh.vertices(i,3)   = 0;
    elseif flag == 1    % deflection
        vmesh.vertices(i,3)   = vmesh.displacement(i,1);  
    elseif flag == 2    % rotation rx
        vmesh.vertices(i,3)   = vmesh.displacement(i,2);  
    elseif flag == 3    % rotation ry
        vmesh.vertices(i,3)   = vmesh.displacement(i,3);  
    elseif flag == 4    % Mx
        vmesh.vertices(i,3)   = vmesh.stress(i,1);  
    elseif flag == 5    % My
        vmesh.vertices(i,3)   = vmesh.stress(i,2);  
    elseif flag == 6    % Mxy
        vmesh.vertices(i,3)   = vmesh.stress(i,3);  
    elseif flag == 7    % Qy
        vmesh.vertices(i,3)   = vmesh.stress(i,4);  
    elseif flag == 8    % Qx
        vmesh.vertices(i,3)   = vmesh.stress(i,5);  
    end
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
    sctrB(1:3:nnElem) = 3*sctr-2;  % deflection   
    sctrB(2:3:nnElem) = 3*sctr-1;  % rotation x   
    sctrB(3:3:nnElem) = 3*sctr;    % rotation y
    [R,dRdparm] = nurbs_derivatives( [xi, eta], circle, mesh );
    jmatrix = dRdparm*exyz(:,1:2); 
    dRdx =  jmatrix \ dRdparm;   
    B(1,3:3:nnElem) = dRdx(1,:);    
    B(2,2:3:nnElem) = -dRdx(2,:);   
    B(3,2:3:nnElem) = -dRdx(1,:);   
    B(3,3:3:nnElem) = dRdx(2,:);
    B(4,1:3:nnElem) = -dRdx(2,:);  
    B(4,2:3:nnElem) = R;
    B(5,1:3:nnElem) = dRdx(1,:);   
    B(5,3:3:nnElem) = R;
    strain = B * u(sctrB);
    stress = D * strain;
    edsp   = u(sctrB);
    edsp   = reshape(edsp, 3, nn);
    vmesh.displacement(count,:) = R * edsp';
    vmesh.vertices(count,1:2) = R*exyz(:,1:2); 
    vmesh.stress(count,:) = stress';    
    if flag == 0        % all fields
        vmesh.vertices(count,3)   = 0;
    elseif flag == 1    % deflection
        vmesh.vertices(count,3)   = vmesh.displacement(count,1);  
    elseif flag == 2    % rotation rx
        vmesh.vertices(count,3)   = vmesh.displacement(count,2);  
    elseif flag == 3    % rotation ry
        vmesh.vertices(count,3)   = vmesh.displacement(count,3);  
    elseif flag == 4    % Mx
        vmesh.vertices(count,3)   = vmesh.stress(count,1);  
    elseif flag == 5    % My
        vmesh.vertices(count,3)   = vmesh.stress(count,2);  
    elseif flag == 6    % Mxy
        vmesh.vertices(count,3)   = vmesh.stress(count,3);  
    elseif flag == 7    % Qy
        vmesh.vertices(count,3)   = vmesh.stress(count,4);  
    elseif flag == 8    % Qx
        vmesh.vertices(count,3)   = vmesh.stress(count,5);  
    end
end

face = polygon.trimesh;
maxnum = max(max(face));
vertices = vmesh.vertices(1:maxnum,:);
displacement = vmesh.displacement(1:maxnum,:);
stress = vmesh.stress(1:maxnum,:);
figure(2);
if flag == 0     % plot all fields
    title( 'Static Bending of a Circular Plate', 'FontSize', 14', 'FontWeight', 'Bold') ;
    axis off;
    linmesh = kntcrv.linmesh;
    colX = linspace( 0, 0.7, 4 );
    rowY = [0.55, 0.15];
    for K = 1:8
        rowId = ceil(K/4);
        colId = K - (rowId-1)*4;
        axes( 'Position', [colX(colId), rowY(rowId), 0.3, 0.3] ) ;
        p = patch('Faces',face, 'Vertices', vertices);
        if K == 1,  cdata = displacement(1:maxnum,1); title('w');
        elseif K == 2,  cdata = displacement(1:maxnum,2);  title('R_x');
        elseif K == 3,  cdata = displacement(1:maxnum,3);  title('R_y');
        elseif K == 4,  cdata = stress(1:maxnum,1);  title('M_x');  
        elseif K == 5,  cdata = stress(1:maxnum,2);  title('M_y');  
        elseif K == 6,  cdata = stress(1:maxnum,3);  title('M_{xy}');  
        elseif K == 7,  cdata = stress(1:maxnum,4);  title('Q_y');  
        elseif K == 8,  cdata = stress(1:maxnum,5);  title('Q_x');  
        end
        set(p,'FaceColor','interp','FaceVertexCData',cdata);
        set(p,'EdgeColor','none');
        hold on;
        for j = 1:size(linmesh,1)
            vv = vmesh.vertices(linmesh(j,:),:);
            plot3(vv(:,1), vv(:,2), vv(:,3), 'k-');
        end 
        axis equal; 
        axis off;
        hold off;
    end
    return;
else    
    p = patch('Faces',face, 'Vertices', vertices);
end
if flag == 1        % deflection
    cdata = displacement(1:maxnum,1);
elseif flag == 2    % rotation rx
    cdata = displacement(1:maxnum,2);
elseif flag == 3    % rotation ry
    cdata = displacement(1:maxnum,3);
elseif flag == 4    % Mx
    cdata = stress(1:maxnum,1);   
elseif flag == 5    % My
    cdata = stress(1:maxnum,2);
elseif flag == 6    % Mxy
    cdata = stress(1:maxnum,3);
elseif flag == 7    % Qy
    cdata = stress(1:maxnum,4);
elseif flag == 8    % Qx
    cdata = stress(1:maxnum,5);
end
set(p,'FaceColor','interp','FaceVertexCData',cdata);
set(p,'EdgeColor','none');

hold on;
linmesh = kntcrv.linmesh;
for j = 1:size(linmesh,1)
    vv = vmesh.vertices(linmesh(j,:),:);
    plot3(vv(:,1), vv(:,2), vv(:,3), 'k-');
end  
hcb = colorbar;
xsize = [min(vmesh.vertices(:,1)), max(vmesh.vertices(:,1))];
ysize = [min(vmesh.vertices(:,2)), max(vmesh.vertices(:,2))];
zsize = [min(vmesh.vertices(:,3)), max(vmesh.vertices(:,3))];
msize = [xsize, ysize, zsize];
axis(msize);
box on;
if flag == 1,  title(hcb,'w'), view(3);
elseif flag == 2,  title(hcb,'r_x');
elseif flag == 3,  title(hcb,'r_y');
elseif flag == 4,  title(hcb,'M_x');
elseif flag == 5,  title(hcb,'M_y');
elseif flag == 6,  title(hcb,'M_{xy}');
elseif flag == 7,  title(hcb,'Q_y');
elseif flag == 8,  title(hcb,'Q_x');
end


