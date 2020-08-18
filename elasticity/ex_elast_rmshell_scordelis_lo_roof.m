%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%  Elastic analysis of Scordelis-Lo roof shell problem based on 
%  Reissner-Mindlin assumption;
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

% define the geometry
R = 25;
L = 50;
theta = 40;
geo = geo_scordelis_lo_roof_shell(L, R, theta);

% material properties
nu   = 0.0;             % Poisson ratio
E    = 4.32e8;          % Youngs modulus 
kapa = 5/6;             % shear correction factor
h    = 0.25;            % thickness   
Dl = E/(1-nu^2)*[  1, nu, 0, 0, 0, 0;
                   nu, 1, 0, 0, 0, 0;
                   0,  0, 0, 0, 0, 0;
                   0, 0, 0, kapa*(1-nu)/2, 0, 0;
                   0, 0, 0, 0, kapa*(1-nu)/2, 0;
                   0, 0, 0, 0, 0, kapa*(1-nu)/2  ];  % elasticity matrix
dof = 6;             % degree of freedom

% build iga mesh structure
mesh = build_iga_mesh( geo );

% find the boundary nodes
bottomNodes = (1:mesh.nCptsU)';
topNodes    = (mesh.nCptsU*(mesh.nCptsV-1)+1:mesh.nCpts)';

dbc = [];    % dbc = [node index, direction, prescribed displacement]
% rigid diaphragm for bottom edge, ux = uz = theta_y = 0
dbc = [dbc; bottomNodes,     ones(size(bottomNodes)),     zeros(size(bottomNodes))];  
dbc = [dbc; bottomNodes,     3*ones(size(bottomNodes)),   zeros(size(bottomNodes))];
dbc = [dbc; bottomNodes,     5*ones(size(bottomNodes)),   zeros(size(bottomNodes))];
% rigid diaphragm for top edge, ux = uz = theta_y = 0
dbc = [dbc; topNodes,        ones(size(topNodes)),        zeros(size(topNodes))];  
dbc = [dbc; topNodes,        3*ones(size(topNodes)),      zeros(size(topNodes))];
dbc = [dbc; topNodes,        5*ones(size(topNodes)),      zeros(size(topNodes))];
% dbc = [dbc; 1, 2, 0]; % consrtrain the rigid body motion

scatdbc = [];
scattbc = [];
if ~isempty(dbc)
    scatdbc = dof * (dbc(:,1)-1) + dbc(:,2);   % scatter dbc
end

% Compute nodal normal vectors for each control points
% Here we use the normal vectors at greville abscissae for the
% corresponding control points
ga = get_greville_abscissae( geo );
elem_index = find_point_span( mesh, ga );
nodalNormal = zeros(mesh.nCpts,3);
for i = 1:size(ga,1)
    gauPts = ga(i,:);
    e      = elem_index(i);
    sctr   = mesh.elNodeCnt(e,:);
    elCpts = mesh.coords(sctr,:);
    [R,ders] = nurbs_derivatives( gauPts, geo, mesh );
    dxdxi   = ders(1,:)*elCpts(:,1:3);
    dxdeta  = ders(2,:)*elCpts(:,1:3);
    n       = cross(dxdxi, dxdeta);
    nodalNormal(i,:) = n/norm(n);
end

% initialize stiffness and force matrices
nDofs = dof * mesh.nCpts;    % total dofs
K = sparse(nDofs,nDofs);     % stiffness matrix 
F = zeros(nDofs,1);          % external force matrix

% use gaussian integration rule
gp_x = mesh.p+1;           % number of integration points in x-direction
gp_y = mesh.q+1;           % number of integration points in y-direction
gp_z = 2;                  % number of integration points in z-direction
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
        pt      = gp(ipt,1:2);      % reference parametric coordinates for each integration point
        wt      = wgt(ipt);       % weigths for each integration point
        zeta    = gp(ipt,3);
        gauPts  = parameter_gauss_mapping( elDoma, pt );       % gauss integration mapping
        j1      = jacobian_gauss_mapping( elDoma );            % jacobian value for gauss mapping   
        [R,ders] = nurbs_derivatives( gauPts, geo, mesh );
        dxdxi = ders*elCpts(:,1:3); 
        e3 = cross(dxdxi(1,:),dxdxi(2,:));
        e3 = e3/norm(e3);   % normal vector 
        ea = dxdxi(1,:) + dxdxi(2,:);
        ea = ea/norm(ea);
        eb = cross(e3,ea);
        eb = eb/norm(eb);
        e1 = sqrt(2)/2 * (ea - eb);  % local basis vector
        e2 = sqrt(2)/2 * (ea + eb);
        q = [e1; e2; e3]';
        Q = local2global_voigt( q );
        jmatrix(1:2,1:3) = dxdxi + h /2 * zeta * ders * nodalNormal(sctr,:);
        jmatrix(3,1:3)   = R * nodalNormal(sctr,:) * h /2;   % dxdzeta
        j2      = det(jmatrix);
        j2_inv = inv(jmatrix); 
        dzetadx = j2_inv(1:3,3);
        ders(3,:) = 0;
        dRdx   =  jmatrix \ ders;

        Ri  = zeta *h/2 *R;
        Rix =  (dzetadx(1)*R + zeta*dRdx(1,:))*h/2;
        Riy =  (dzetadx(2)*R + zeta*dRdx(2,:))*h/2;
        Riz =  (dzetadx(3)*R + zeta*dRdx(3,:))*h/2;      
        
        B = zeros(6,nnElem);    
        B(1,1:6:nnElem) =    dRdx(1,:);
        B(1,5:6:nnElem) =    Rix .* nodalNormal(sctr,3)';
        B(1,6:6:nnElem) =  - Rix .* nodalNormal(sctr,2)';
        
        B(2,2:6:nnElem) =    dRdx(2,:);
        B(2,4:6:nnElem) =  - Riy .* nodalNormal(sctr,3)';
        B(2,6:6:nnElem) =    Riy .* nodalNormal(sctr,1)';
        
        B(3,3:6:nnElem) =    dRdx(3,:);
        B(3,4:6:nnElem) =    Riz .* nodalNormal(sctr,2)';
        B(3,5:6:nnElem) =  - Riz .* nodalNormal(sctr,1)';
        
        B(4,1:6:nnElem) =    dRdx(2,:);
        B(4,2:6:nnElem) =    dRdx(1,:);
        B(4,4:6:nnElem) =  - Rix .* nodalNormal(sctr,3)';
        B(4,5:6:nnElem) =    Riy .* nodalNormal(sctr,3)';
        B(4,6:6:nnElem) =    Rix .* nodalNormal(sctr,1)' - Riy .* nodalNormal(sctr,2)';        
        
        B(5,2:6:nnElem) =   dRdx(3,:);
        B(5,3:6:nnElem) =   dRdx(2,:);
        B(5,4:6:nnElem) =   Riy .* nodalNormal(sctr,2)' - Riz .* nodalNormal(sctr,3)';
        B(5,5:6:nnElem) = - Riy .* nodalNormal(sctr,1)';
        B(5,6:6:nnElem) =   Riz .* nodalNormal(sctr,1)';       
        
        B(6,1:6:nnElem) =   dRdx(3,:);
        B(6,3:6:nnElem) =   dRdx(1,:);
        B(6,4:6:nnElem) =   Rix .* nodalNormal(sctr,2)';
        B(6,5:6:nnElem) =   Riz .* nodalNormal(sctr,3)' - Rix .* nodalNormal(sctr,1)';
        B(6,6:6:nnElem) = - Riz .* nodalNormal(sctr,2)';

        Dg = Q * Dl * Q';
        fac = j1 *j2 * wt;      
        K(sctrB,sctrB) = K(sctrB,sctrB) + B' * Dg * B * fac;
    end 
end

% enforce force conditions
% use gaussian integration rule
gp_x = mesh.p+1;           % number of integration points in x-direction
gp_y = mesh.q+1;           % number of integration points in y-direction
[gp, wgt] = gauss_quadrature(gp_x, gp_y);   % calculate integration points and its weights
for e = 1:mesh.nElems                 % loop over elements
    sctr   = mesh.elNodeCnt(e,:);     % element control points index
    elDoma = mesh.elDoma(e,:);        % element parametric domain
    elCpts = mesh.coords(sctr,:);     % coordinates of element control points
    strf   = sctr*6-3;
    for ipt = 1:size(gp,1)        % loop over integration points
        pt       = gp(ipt,:);      % reference parametric coordinates for each integration point
        wt       = wgt(ipt);       % weigths for each integration point
        gauPts   = parameter_gauss_mapping( elDoma, pt );       % gauss integration mapping
        j1       = jacobian_gauss_mapping( elDoma );            % jacobian value for gauss mapping   
        [R,ders] = nurbs_derivatives( gauPts, geo, mesh );
        jmatrix  = ders*elCpts(:,1:3); 
        j2       = norm(cross(jmatrix(1,:),jmatrix(2,:)));
        fac      = j1 *j2 * wt;        
        g        = -90;
        F(strf)  = F(strf) +  R' * g * fac;
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
factor = len/max(abs(u(3:6:end)))/10;  % scale factor for better visualization

% post-processing
num1 = 50;  % if you want to obtain a more smooth results, please provide large numbers for num1 and num2    
num2 = 50;
polygon = build_visual_mesh_suf( num1, num2 );        % build visualized mesh
kntcrv = build_visual_knotcurve_suf( mesh.uKnots, mesh.vKnots, num1+1 ); % build visualized knot curves
numpts = (num1+1)*(num2+1) + size(kntcrv.linpts,1);
vmesh.vertices = zeros(numpts,3);
vmesh.displacement = zeros(numpts,6);
vmesh.stress = zeros(numpts,6);
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
    [R,ders] = nurbs_derivatives( [xi, eta], geo, mesh );
    dxdxi = ders*elCpts(:,1:3); 
    e3 = cross(dxdxi(1,:),dxdxi(2,:));
    e3 = e3/norm(e3);   % normal vector 
    ea = dxdxi(1,:) + dxdxi(2,:);
    ea = ea/norm(ea);
    eb = cross(e3,ea);
    eb = eb/norm(eb);
    e1 = sqrt(2)/2 * (ea - eb);  % local basis vector
    e2 = sqrt(2)/2 * (ea + eb);
    q = [e1; e2; e3]';
    Q = local2global_voigt( q );
    jmatrix(1:2,1:3) = dxdxi + h /2 * zeta * ders * nodalNormal(sctr,:);
    jmatrix(3,1:3)   = R * nodalNormal(sctr,:) * h /2;   % dxdzeta
    j2      = det(jmatrix);
    j2_inv = inv(jmatrix); 
    dzetadx = j2_inv(1:3,3);
    ders(3,:) = 0;
    dRdx   =  jmatrix \ ders;

    Ri  = zeta *h/2 *R;
    Rix =  (dzetadx(1)*R + zeta*dRdx(1,:))*h/2;
    Riy =  (dzetadx(2)*R + zeta*dRdx(2,:))*h/2;
    Riz =  (dzetadx(3)*R + zeta*dRdx(3,:))*h/2;      

    B = zeros(6,nnElem);    
    B(1,1:6:nnElem) =    dRdx(1,:);
    B(1,5:6:nnElem) =    Rix .* nodalNormal(sctr,3)';
    B(1,6:6:nnElem) =  - Rix .* nodalNormal(sctr,2)';

    B(2,2:6:nnElem) =    dRdx(2,:);
    B(2,4:6:nnElem) =  - Riy .* nodalNormal(sctr,3)';
    B(2,6:6:nnElem) =    Riy .* nodalNormal(sctr,1)';

    B(3,3:6:nnElem) =    dRdx(3,:);
    B(3,4:6:nnElem) =    Riz .* nodalNormal(sctr,2)';
    B(3,5:6:nnElem) =  - Riz .* nodalNormal(sctr,1)';

    B(4,1:6:nnElem) =    dRdx(2,:);
    B(4,2:6:nnElem) =    dRdx(1,:);
    B(4,4:6:nnElem) =  - Rix .* nodalNormal(sctr,3)';
    B(4,5:6:nnElem) =    Riy .* nodalNormal(sctr,3)';
    B(4,6:6:nnElem) =    Rix .* nodalNormal(sctr,1)' - Riy .* nodalNormal(sctr,2)';        

    B(5,2:6:nnElem) =   dRdx(3,:);
    B(5,3:6:nnElem) =   dRdx(2,:);
    B(5,4:6:nnElem) =   Riy .* nodalNormal(sctr,2)' - Riz .* nodalNormal(sctr,3)';
    B(5,5:6:nnElem) = - Riy .* nodalNormal(sctr,1)';
    B(5,6:6:nnElem) =   Riz .* nodalNormal(sctr,1)';       

    B(6,1:6:nnElem) =   dRdx(3,:);
    B(6,3:6:nnElem) =   dRdx(1,:);
    B(6,4:6:nnElem) =   Rix .* nodalNormal(sctr,2)';
    B(6,5:6:nnElem) =   Riz .* nodalNormal(sctr,3)' - Rix .* nodalNormal(sctr,1)';
    B(6,6:6:nnElem) = - Riz .* nodalNormal(sctr,2)';

    Dg = Q * Dl * Q';
    strain = B * u(sctrB);
    stress = Dg * strain;
    edsp   = u(sctrB);
    edsp   = reshape(edsp, dof, nn);
    vmesh.displacement(i,:) = R * edsp';
    vmesh.vertices(i,1:3) = R*exyz(:,1:3) + vmesh.displacement(i,1:3)*factor;
    vmesh.stress(i,:) = stress';
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
    [R,ders] = nurbs_derivatives( [xi, eta], geo, mesh );
    jmatrix = ders*elCpts(:,1:3); 
    dxdxi = ders*elCpts(:,1:3); 
    e3 = cross(dxdxi(1,:),dxdxi(2,:));
    e3 = e3/norm(e3);   % normal vector 
    ea = dxdxi(1,:) + dxdxi(2,:);
    ea = ea/norm(ea);
    eb = cross(e3,ea);
    eb = eb/norm(eb);
    e1 = sqrt(2)/2 * (ea - eb);  % local basis vector
    e2 = sqrt(2)/2 * (ea + eb);
    q = [e1; e2; e3]';
    Q = local2global_voigt( q );
    jmatrix(1:2,1:3) = dxdxi + h /2 * zeta * ders * nodalNormal(sctr,:);
    jmatrix(3,1:3)   = R * nodalNormal(sctr,:) * h /2;   % dxdzeta
    j2      = det(jmatrix);
    j2_inv = inv(jmatrix); 
    dzetadx = j2_inv(1:3,3);
    ders(3,:) = 0;
    dRdx   =  jmatrix \ ders;

    Ri  = zeta *h/2 *R;
    Rix =  (dzetadx(1)*R + zeta*dRdx(1,:))*h/2;
    Riy =  (dzetadx(2)*R + zeta*dRdx(2,:))*h/2;
    Riz =  (dzetadx(3)*R + zeta*dRdx(3,:))*h/2;      

    B = zeros(6,nnElem);    
    B(1,1:6:nnElem) =    dRdx(1,:);
    B(1,5:6:nnElem) =    Rix .* nodalNormal(sctr,3)';
    B(1,6:6:nnElem) =  - Rix .* nodalNormal(sctr,2)';

    B(2,2:6:nnElem) =    dRdx(2,:);
    B(2,4:6:nnElem) =  - Riy .* nodalNormal(sctr,3)';
    B(2,6:6:nnElem) =    Riy .* nodalNormal(sctr,1)';

    B(3,3:6:nnElem) =    dRdx(3,:);
    B(3,4:6:nnElem) =    Riz .* nodalNormal(sctr,2)';
    B(3,5:6:nnElem) =  - Riz .* nodalNormal(sctr,1)';

    B(4,1:6:nnElem) =    dRdx(2,:);
    B(4,2:6:nnElem) =    dRdx(1,:);
    B(4,4:6:nnElem) =  - Rix .* nodalNormal(sctr,3)';
    B(4,5:6:nnElem) =    Riy .* nodalNormal(sctr,3)';
    B(4,6:6:nnElem) =    Rix .* nodalNormal(sctr,1)' - Riy .* nodalNormal(sctr,2)';        

    B(5,2:6:nnElem) =   dRdx(3,:);
    B(5,3:6:nnElem) =   dRdx(2,:);
    B(5,4:6:nnElem) =   Riy .* nodalNormal(sctr,2)' - Riz .* nodalNormal(sctr,3)';
    B(5,5:6:nnElem) = - Riy .* nodalNormal(sctr,1)';
    B(5,6:6:nnElem) =   Riz .* nodalNormal(sctr,1)';       

    B(6,1:6:nnElem) =   dRdx(3,:);
    B(6,3:6:nnElem) =   dRdx(1,:);
    B(6,4:6:nnElem) =   Rix .* nodalNormal(sctr,2)';
    B(6,5:6:nnElem) =   Riz .* nodalNormal(sctr,3)' - Rix .* nodalNormal(sctr,1)';
    B(6,6:6:nnElem) = - Riz .* nodalNormal(sctr,2)';

    Dg = Q * Dl * Q';
    strain = B * u(sctrB);
    stress = Dg * strain;
    edsp   = u(sctrB);
    edsp   = reshape(edsp, dof, nn);
    vmesh.displacement(count,:) = R * edsp';
    vmesh.vertices(count,1:3) = R*exyz(:,1:3);
%     vmesh.vertices(count,1:3) = R*exyz(:,1:3) + vmesh.displacement(count,1:3)*factor; 
    vmesh.stress(count,:) = stress';    
end

% plot approximate solutions
figure;
face = polygon.trimesh;
maxnum = max(max(face));
vertices = vmesh.vertices(1:maxnum,:);
displacement = vmesh.displacement(1:maxnum,:);
stress = vmesh.stress(1:maxnum,:);
linmesh = kntcrv.linmesh;
title( 'Approximate Solutions', 'FontSize', 14', 'FontWeight', 'Bold') ;
p = patch('Faces',face, 'Vertices', vertices);
cdata = displacement(1:maxnum,3); title('U_z');
set(p,'FaceColor','interp','FaceVertexCData',cdata);
set(p,'EdgeColor','none');
hold on;
for j = 1:size(linmesh,1) % undeformed parametric curve
    vv = vmesh.vertices(linmesh(j,:),:);
    plot3(vv(:,1), vv(:,2),vv(:,3), 'Color',[0.5,0.5,0.5]);
end 
for j = 1:size(linmesh,1) % deformed parametric curve
    vv = vmesh.vertices(linmesh(j,:),:) + vmesh.displacement(linmesh(j,:),1:3)*factor ;
    plot3(vv(:,1), vv(:,2),vv(:,3), 'k-');
end
colorbar;
axis equal; 
view(3);


% output the vertical displacement of the midside point
pt = [1e-15, 0.5];
[R,ders] = nurbs_derivatives( pt, geo, mesh );
enum = find_point_span( mesh, pt );
sctr = mesh.elNodeCnt(enum,:); % element control points index
exyz = mesh.coords(sctr,:); % element control points' coordinates
vd = R * u(sctr*6-3)
nn = [mesh.nElems, mesh.nCpts]


