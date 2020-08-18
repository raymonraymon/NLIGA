%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%  Elastic analysis of the micro-motor rotor;
%  The model is built with eight NURBS patches.
%  Note: here for multipatch problem, we only condider the multiple NURBS patches
%  with same degrees, knot vector and control points on the interfaces.
%
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%


clc; 
clear all;

E            = 169e3;  % Youngs modulus
nu           = 0.262;  % Poisson ratio
q0           = 10e-3;  % uniform tension
D0 = E/(1-nu^2)*[  1      nu        0;
                  nu     1          0;
                  0      0      (1-nu)/2  ];
dof = 2;             % degree of freedom

% define the geometry
r1 = 8*10^(-3);
r2 = 25*10^-3;
r3 = 50*10^-3;
% r1 = 8;
% r2 = 25;
% r3 = 50;
alpha2 = 20;
alpha3 = 25;
geo = geo_micro_motor_rotor_8patch( r1, r2, r3, alpha2, alpha3 );

% build iga mesh structure
mesh = cell(1,length(geo));
for i = 1:length(geo)
    mesh{1,i} = build_iga_mesh( geo{1,i} );   
    mesh{1,i}.gloElNodeCnt = mesh{1,i}.elNodeCnt;
    mesh{1,i}.nodeNum = zeros(mesh{1,i}.nCptsV, mesh{1,i}.nCptsU);
end

% node numbering
for j = 1:mesh{1,1}.nCptsV   % node numbering for patch 1
    mesh{1,1}.nodeNum(j,:) = (j-1)*(mesh{1,1}.nCptsU) + (1:mesh{1,1}.nCptsU);
    mesh{1,1}.globCoords = mesh{1,1}.coords;
end
for kk = 2:5  %  node numbering for patch 2 to 5
    mesh{1,kk}.nodeNum(:,1) = mesh{1,kk-1}.nodeNum(:,end);
    mesh{1,kk}.globCoords = zeros(mesh{1,kk}.nCptsV*(mesh{1,kk}.nCptsU-1),4);
    for j = 1:mesh{1,kk}.nCptsV
        for i = 2:mesh{1,kk}.nCptsU
            aa = (j-1)*(mesh{1,kk}.nCptsU)+i;
            bb = (j-1)*(mesh{1,kk}.nCptsU-1)+i-1;
            mesh{1,kk}.nodeNum(j,i) = bb + mesh{1,kk-1}.nodeNum(end,end);
            mesh{1,kk}.globCoords(bb,:) = mesh{1,kk}.coords(aa,:); 
        end
    end
end
mesh{1,6}.nodeNum(1,:) = mesh{1,1}.nodeNum(end,:);  %  node numbering for patch 6
mesh{1,6}.globCoords = zeros((mesh{1,6}.nCptsV-1)*mesh{1,6}.nCptsU,4);
for j = 2:mesh{1,6}.nCptsV
    aa = (j-1)*(mesh{1,6}.nCptsU) + (1:mesh{1,6}.nCptsU);
    bb = (j-2)*(mesh{1,6}.nCptsU) + (1:mesh{1,6}.nCptsU);
    mesh{1,6}.nodeNum(j,:) = bb + mesh{1,5}.nodeNum(end,end);
    mesh{1,6}.globCoords(bb,:) = mesh{1,6}.coords(aa,:); 
end
mesh{1,7}.nodeNum(1,:) = mesh{1,3}.nodeNum(end,:);  %  node numbering for patch 7
mesh{1,7}.globCoords = zeros((mesh{1,7}.nCptsV-1)*mesh{1,7}.nCptsU,4);
for j = 2:mesh{1,7}.nCptsV
    aa = (j-1)*(mesh{1,7}.nCptsU) + (1:mesh{1,7}.nCptsU);
    bb = (j-2)*(mesh{1,7}.nCptsU) + (1:mesh{1,7}.nCptsU);
    mesh{1,7}.nodeNum(j,:) = bb + mesh{1,6}.nodeNum(end,end);
    mesh{1,7}.globCoords(bb,:) = mesh{1,7}.coords(aa,:); 
end
mesh{1,8}.nodeNum(1,:) = mesh{1,5}.nodeNum(end,:);  %  node numbering for patch 7
mesh{1,8}.globCoords = zeros((mesh{1,8}.nCptsV-1)*mesh{1,8}.nCptsU,4);
for j = 2:mesh{1,8}.nCptsV
    aa = (j-1)*(mesh{1,8}.nCptsU) + (1:mesh{1,8}.nCptsU);
    bb = (j-2)*(mesh{1,8}.nCptsU) + (1:mesh{1,8}.nCptsU);
    mesh{1,8}.nodeNum(j,:) = (j-2)*(mesh{1,8}.nCptsU) + (1:mesh{1,8}.nCptsU) + mesh{1,7}.nodeNum(end,end);
    mesh{1,8}.globCoords(bb,:) = mesh{1,8}.coords(aa,:); 
end

% reconstruct global nodes coordinates
globNodeNum = mesh{1,end}.nodeNum(end,end);
globCoords = [];
for kk = 1:length(mesh)
    globCoords = [globCoords; mesh{1,kk}.globCoords];
end

% reconstruct global connectivity
for kk = 1:length(mesh)
    for j = 1:mesh{1,kk}.nElemV
        for i = 1:mesh{1,kk}.nElemU
            aa = (j-1)*mesh{1,kk}.nElemU + i;
            mesh{1,kk}.gloElNodeCnt(aa,:) = reshape(mesh{1,kk}.nodeNum(j:j+mesh{1,kk}.q, i:i+mesh{1,kk}.p)', 1, mesh{1,kk}.nElemCpts );
        end
    end
end

% find the boundary nodes
epsilon=1e-10; 
bottomNodes = find(abs(globCoords(:,2))<epsilon);
leftNodes   = find(abs(globCoords(:,1))<epsilon);
innerNodes  = [];
for kk = 1:5
    innerNodes = [innerNodes; mesh{1,kk}.nodeNum(1,:)'];
end
innerNodes = unique(innerNodes);

outerNodes = [];
for kk = 6:8
    outerNodes = [outerNodes; mesh{1,kk}.nodeNum(end,:)'];
end
outerNodes = unique(outerNodes);


% hold on;
% plot3(globCoords(innerNodes,1),globCoords(innerNodes,2),globCoords(innerNodes,3),'ro')

dbc = [];    % dbc = [node index, direction, prescribed displacement]
dbc = [dbc; leftNodes,     ones(size(leftNodes)),   zeros(size(leftNodes))];
dbc = [dbc; bottomNodes,     2*ones(size(bottomNodes)),   zeros(size(bottomNodes))];
dbc = [dbc; innerNodes,     ones(size(innerNodes)),   zeros(size(innerNodes))];
dbc = [dbc; innerNodes,     2*ones(size(innerNodes)),   zeros(size(innerNodes))];
scatdbc = [];
scattbc = [];
if ~isempty(dbc)
    scatdbc = dof * (dbc(:,1)-1) + dbc(:,2);   % scatter dbc
end

% initialize stiffness and force matrices
nDofs = dof * globNodeNum;    % total dofs
K = sparse(nDofs,nDofs);     % stiffness matrix 
F = zeros(nDofs,1);          % external force matrix
% use gaussian integration rule
gp_x = mesh{1,1}.p+1;           % number of integration points in x-direction
gp_y = mesh{1,1}.q+1;           % number of integration points in y-direction
[gp, wgt] = gauss_quadrature(gp_x, gp_y);   % calculate integration points and its weights
for jj = 1:length(mesh)
    for e = 1:mesh{1,jj}.nElems                 % loop over elements
        sctr   = mesh{1,jj}.gloElNodeCnt(e,:);  % element control points index
        elDoma = mesh{1,jj}.elDoma(e,:);        % element parametric domain
        elCpts = globCoords(sctr,:);     % coordinates of element control points
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
            [R,ders] = nurbs_derivatives( gauPts, geo{1,jj}, mesh{1,jj} );
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
end

% boundary conditions imposition on outer boundaries
forceElems = cell(1,3);
forceKnots = cell(1,3);
forceNodes = cell(1,3);
fpNum = 6:8;
for ii = 1:3
    [forceElems{1,ii}, forceKnots{1,ii}] = build_knot_connectivity( mesh{1,fpNum(ii)}.uKnots );
    forceNodes{1,ii}  = mesh{1,fpNum(ii)}.nodeNum(end,:);
    for jj = 1:size(forceElems{1,ii},1)
        forceElems{1,ii}(jj,:) = forceNodes{1,ii}(forceElems{1,ii}(jj,:));
    end
end


for kk = 1:3
    gp_x = mesh{1,kk}.p+1;
    [gp, wgt] = gauss_quadrature(gp_x);
    for e = 1:size(forceElems{1,kk},1)
        sctr   = forceElems{1,kk}(e,:);  
        elDoma = forceKnots{1,kk}(e,:);
        elCpts = globCoords(sctr,:);
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
            [N, dNdxi] = nrbNurbs1DBasisDerivs(gauPts, mesh{1,fpNum(kk)}.p, mesh{1,fpNum(kk)}.uKnots, globCoords(forceNodes{1,kk}',4));
            j2     = dNdxi*elCpts(:,1:2);            
            j2     = norm(j2);                
            x      = N*elCpts(:,1:2);  % compute coordinates
            theta = atan(x(2)/x(1));
            tx     = -q0*x(1)/norm(x);
            ty     = -q0*x(2)/norm(x);            
            R = zeros(2,nnElem);
            R(1,1:2:nnElem) = N;
            R(2,2:2:nnElem) = N;
            fac = j1 * j2 * wt;
            F(sctrB) = F(sctrB) + R'* [tx; ty] * fac;
        end
    end
end

% solving equations
ndbc = size(dbc,1);   % number of displacement constrained nodes
K(scatdbc,:) = zeros(ndbc, nDofs);
K(scatdbc,scatdbc) = eye(ndbc);
F(scatdbc,:) = dbc(:,3);        
u = K\F;
min_x = min(globCoords(:,1));
max_x = max(globCoords(:,1));
min_y = min(globCoords(:,2));
max_y = max(globCoords(:,2));
len = sqrt((max_x-min_x)^2 + (max_y-min_y)^2);
% factor = len/max(abs(u))/15;  % scale factor for better visualization
factor = 1.5e6; 
% factor = 1; 
% post-processing
vmesh = cell(1,length(mesh));
polygon = cell(1,length(mesh));
kntcrv = cell(1,length(mesh)); 
num1 = 50;  % if you want to obtain a more smooth results, please provide large numbers for num1 and num2    
num2 = 50;
for kk = 1:length(mesh)
    polygon{1,kk} = build_visual_mesh_suf( num1, num2 );        % build visualized mesh
    kntcrv{1,kk} = build_visual_knotcurve_suf( mesh{1,kk}.uKnots, mesh{1,kk}.vKnots, num1+1 ); % build visualized knot curves
    numpts = (num1+1)*(num2+1) + size(kntcrv{1,kk}.linpts,1);
    vmesh{1,kk}.vertices = zeros(numpts,2);
    vmesh{1,kk}.displacement = zeros(numpts,2);
    vmesh{1,kk}.exactdisplacement = zeros(numpts,2);
    vmesh{1,kk}.stress = zeros(numpts,3);
    vmesh{1,kk}.exactstress = zeros(numpts,3);
    elem_index = find_point_span( mesh{1,kk}, polygon{1,kk}.tripts );
    for i=1:(num1+1)*(num2+1)
        xi = polygon{1,kk}.tripts(i,1);   % u mesh point
        eta = polygon{1,kk}.tripts(i,2);  % v mesh point
        e = elem_index(i);          % element number
        sctr = mesh{1,kk}.gloElNodeCnt(e,:); % element control points index
        exyz = globCoords(sctr,:); % element control points' coordinates
        nn = numel(sctr);           % number of control points for each element
        nnElem = nn*dof;            % dof for each element
        sctrB = zeros(1, nnElem); 
        for k = 1:dof
            sctrB(k:dof:nnElem) = dof*(sctr-1) + k;  % displacement in i-th direction
        end
        [R,dRdparm] = nurbs_derivatives( [xi, eta], geo{1,kk}, mesh{1,kk} );
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
        vmesh{1,kk}.displacement(i,:) = R * edsp';
        vmesh{1,kk}.vertices(i,1:2) = R*exyz(:,1:2)  + vmesh{1,kk}.displacement(i,1:2)*factor; 
        vmesh{1,kk}.stress(i,:) = stress';
    end
    count = (num1+1)*(num2+1);
    line_index = find_point_span( mesh{1,kk}, kntcrv{1,kk}.linpts );
    kntcrv{1,kk}.linmesh = kntcrv{1,kk}.linmesh + (num1+1)*(num2+1);
    for i = 1:size(kntcrv{1,kk}.linpts,1)
        count = count+1;
        xi = kntcrv{1,kk}.linpts(i,1);
        eta = kntcrv{1,kk}.linpts(i,2);
        e = line_index(i);   % element number
        sctr = mesh{1,kk}.gloElNodeCnt(e,:);     % element control points index
        exyz = globCoords(sctr,:);  % element control points' coordinates
        nn = numel(sctr);           % number of control points for each element
        nnElem = nn*dof;            % dof for each element
        sctrB = zeros(1, nnElem); 
        for k = 1:dof
            sctrB(k:dof:nnElem) = dof*(sctr-1) + k;  % displacement in i-th direction
        end
        [R,dRdparm] = nurbs_derivatives( [xi, eta], geo{1,kk}, mesh{1,kk} );
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
        vmesh{1,kk}.displacement(count,:) = R * edsp';
        vmesh{1,kk}.vertices(count,1:2) = R*exyz(:,1:2)  + vmesh{1,kk}.displacement(count,1:2)*factor;  
        vmesh{1,kk}.stress(count,:) = stress';    
    end
end

% plot the solutions
figure;
title( 'Static Analysis of the Micro-Motor Rotor', 'FontSize', 14', 'FontWeight', 'Bold') ;
axis off;
colX = linspace( 0.1, 0.65, 3 );
rowY = [0.55, 0.15];
for K = 1:6
    rowId = ceil(K/3);
    colId = K - (rowId-1)*3;
    axes( 'Position', [colX(colId), rowY(rowId), 0.3, 0.3] ) ;
    for jj = 1:length(mesh)
        face = polygon{1,jj}.trimesh;
        maxnum = max(max(face));
        vertices = vmesh{1,jj}.vertices(1:maxnum,:);
        displacement = vmesh{1,jj}.displacement(1:maxnum,:);
        stress = vmesh{1,jj}.stress(1:maxnum,:);
        linmesh = kntcrv{1,jj}.linmesh;
        p = patch('Faces',face, 'Vertices', vertices);
        if K == 1,  cdata = displacement(1:maxnum,1); title('U_x');
        elseif K == 2,  cdata = displacement(1:maxnum,2);  title('U_y');
        elseif K == 3,  cdata = sqrt( displacement(1:maxnum,1).^2 + displacement(1:maxnum,2).^2  ); title('U Magnitude');
        elseif K == 4,  cdata = stress(1:maxnum,1);  title('\sigma_{xx}');  
        elseif K == 5,  cdata = stress(1:maxnum,2);  title('\sigma_{yy}');  
        elseif K == 6,  cdata = stress(1:maxnum,3);  title('\sigma_{xy}');   
        end
        hold on;
        set(p,'FaceColor','interp','FaceVertexCData',cdata);
        set(p,'EdgeColor','none');
        for j = 1:size(linmesh,1)
            hold on;
            vv = vmesh{1,jj}.vertices(linmesh(j,:),:);
            plot(vv(:,1), vv(:,2), 'k-');
        end 
        axis equal;
        axis off;
    end
    hold off;
    colorbar;
end

% plot von-Mises stress
figure;
for jj = 1:length(mesh)
    face = polygon{1,jj}.trimesh;
    maxnum = max(max(face));
    vertices = vmesh{1,jj}.vertices(1:maxnum,:);
    displacement = vmesh{1,jj}.displacement(1:maxnum,:);
    stress = vmesh{1,jj}.stress(1:maxnum,:);
    linmesh = kntcrv{1,jj}.linmesh;
    p = patch('Faces',face, 'Vertices', vertices);
%     cdata = sqrt( displacement(1:maxnum,1).^2 + displacement(1:maxnum,2).^2  ); 
%     title('Displacement Magnitude');
    cdata = von_mises( stress(1:maxnum,:) );
    title('Von-Mises Stress (No Average)');
    hold on;
    set(p,'FaceColor','interp','FaceVertexCData',cdata);
    set(p,'EdgeColor','none');
    for j = 1:size(linmesh,1)
        hold on;
        vv = vmesh{1,jj}.vertices(linmesh(j,:),:);
        plot(vv(:,1), vv(:,2), 'k-');
    end 
    axis equal;
    axis off;
end
hold off;
colorbar;

% output
nElems = 0;
for jj = 1:length(mesh)
    nElems = nElems + mesh{1,jj}.nElems;
end
nNodes = globNodeNum;
p2n = (num1+1)*(num2+1);
p3n = (num1+1)*num2+1;
p7n = 1;
p2s = von_mises(vmesh{1,2}.stress(p2n,:));
p3s = von_mises(vmesh{1,3}.stress(p3n,:));
p7s = von_mises(vmesh{1,7}.stress(p7n,:));

nn = [nElems, nNodes]
format long
stress = [p2s, p3s, p7s]

maxd = max(sqrt( displacement(1:maxnum,1).^2 + displacement(1:maxnum,2).^2  ))


