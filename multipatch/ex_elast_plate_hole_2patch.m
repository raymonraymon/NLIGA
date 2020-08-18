%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%  Benchmark problem: tension of the plate with hole;
%  The model is built with two NURBS patches.
%  Here for multipatch problem, we only condider the multiple NURBS patches
%  with same degrees, knot vector and control points on the interfaces.
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
plate_hole = geo_plate_with_hole_2patch(L,a);

% build iga mesh structure
mesh = cell(1,length(plate_hole));
for i = 1:length(plate_hole)
    mesh{1,i} = build_iga_mesh( plate_hole{1,i} );   
    mesh{1,i}.gloElNodeCnt = mesh{1,i}.elNodeCnt;
    mesh{1,i}.nodeNum = zeros(mesh{1,i}.nCptsV, mesh{1,i}.nCptsU);
end

% node numbering
for j = 1:mesh{1,1}.nCptsV
    for i = 1:mesh{1,1}.nCptsU
        mesh{1,1}.nodeNum(j,i) = (j-1)*(mesh{1,1}.nCptsU)+i;
    end
end

mesh{1,2}.nodeNum(:,1) = mesh{1,1}.nodeNum(:,end);
for j = 1:mesh{1,2}.nCptsV
    for i = 2:mesh{1,2}.nCptsU
        mesh{1,2}.nodeNum(j,i) = (j-1)*(mesh{1,2}.nCptsU-1)+i-1 + mesh{1,1}.nCpts;
    end
end

% reconstruct global nodes coordinates
globNodeNum = mesh{1,1}.nCpts + mesh{1,2}.nCpts - mesh{1,2}.nCptsV;
globCoords = zeros(globNodeNum,4);
globCoords(1:mesh{1,1}.nCpts,:) = mesh{1,1}.coords;
count = mesh{1,1}.nCpts;
for j = 1:mesh{1,2}.nCptsV
    for i = 2:mesh{1,2}.nCptsU
        count = count+1;
        k = (j-1)*(mesh{1,2}.nCptsU) +i;
        globCoords(count,:) = mesh{1,2}.coords(k,:);
    end
end

% build global node element connectivity
for j = 1:mesh{1,2}.nElemV
    for i = 1:mesh{1,2}.nElemU
        k = (j-1)*mesh{1,2}.nElemU + i;
        mesh{1,2}.gloElNodeCnt(k,:) = reshape(mesh{1,2}.nodeNum(j:j+mesh{1,2}.q, i:i+mesh{1,2}.p)', 1, mesh{1,2}.nElemCpts );
    end
end

% find the boundary nodes
bottomNodes = mesh{1,2}.nodeNum(:,end);
topNodes    = mesh{1,1}.nodeNum(end,:);
leftNodes   = mesh{1,1}.nodeNum(:,1);
rightNodes  = mesh{1,2}.nodeNum(end,:);

dbc = [];    % dbc = [node index, direction, prescribed displacement]
dbc = [dbc; leftNodes,     ones(size(leftNodes)),   zeros(size(leftNodes))];
dbc = [dbc; bottomNodes,     2*ones(size(bottomNodes)),   zeros(size(bottomNodes))];
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
            [R,ders] = nurbs_derivatives( gauPts, plate_hole{1,jj}, mesh{1,jj} );
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

% boundary conditions imposition on top boundary
[forceElems1, forceKnots1] = build_knot_connectivity( mesh{1,1}.uKnots );
forceElems1 = forceElems1 + (mesh{1,1}.nCptsV-1)*mesh{1,1}.nCptsU;
forceNodes1  = (mesh{1,1}.nCptsV-1)*mesh{1,1}.nCptsU+1 : mesh{1,1}.nCpts;
gp_x = mesh{1,1}.p+1;
[gp, wgt] = gauss_quadrature(gp_x);
for e = 1:size(forceElems1,1)
    sctr   = forceElems1(e,:);  
    elDoma = forceKnots1(e,:);
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
        [N, dNdxi] = nrbNurbs1DBasisDerivs(gauPts, mesh{1,1}.p, mesh{1,1}.uKnots, globCoords(forceNodes1',4));
        j2     = dNdxi*elCpts(:,1:2);            
        j2     = norm(j2);                
        x      = N*elCpts(:,1:2);  % compute exact tractions
        [str,disp]    = exact_solution_plate_hole(x,a,E,nu);
        tx  = q0*str(3);
        ty  = q0*str(2);        
        R = zeros(2,nnElem);
        R(1,1:2:nnElem) = N;
        R(2,2:2:nnElem) = N;
        fac = j1 * j2 * wt;
        F(sctrB) = F(sctrB) + R'* [tx; ty] * fac;
    end   
end

% boundary conditions imposition on right boundary
[forceElems2, forceKnots2] = build_knot_connectivity( mesh{1,2}.uKnots );
forceNodes2  = mesh{1,2}.nodeNum(end,:);
for i = 1:size(forceElems2,1)
    forceElems2(i,:) = forceNodes2(forceElems2(i,:))';
end
for e = 1:size(forceElems2,1)
    sctr   = forceElems2(e,:);  
    elDoma = forceKnots2(e,:);
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
        [N, dNdxi] = nrbNurbs1DBasisDerivs(gauPts, mesh{1,2}.p, mesh{1,2}.uKnots, globCoords(forceNodes2',4));
        j2     = dNdxi*elCpts(:,1:2);            
        j2     = norm(j2);               
        x      = N*elCpts(:,1:2);  % compute exact tractions
        [str,disp]    = exact_solution_plate_hole(x,a,E,nu);
        tx  = q0*str(1);
        ty  = q0*str(3);        
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
vmesh = cell(1,length(mesh));
polygon = cell(1,length(mesh));
kntcrv = cell(1,length(mesh)); 
for kk = 1:length(mesh)
    num1 = 50;  % if you want to obtain a more smooth results, please provide large numbers for num1 and num2    
    num2 = 50;
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
        [R,dRdparm] = nurbs_derivatives( [xi, eta], plate_hole{1,kk}, mesh{1,kk} );
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
        vmesh{1,kk}.vertices(i,1:2) = R*exyz(:,1:2);
        vmesh{1,kk}.stress(i,:) = stress';
        [str,disp]    = exact_solution_plate_hole(vmesh{1,kk}.vertices(i,1:2),a,E,nu);
        vmesh{1,kk}.exactstress(i,:) = q0*str;
        vmesh{1,kk}.exactdisplacement(i,:) = q0*disp;
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
        [R,dRdparm] = nurbs_derivatives( [xi, eta], plate_hole{1,kk}, mesh{1,kk} );
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
        vmesh{1,kk}.vertices(count,1:2) = R*exyz(:,1:2); 
        vmesh{1,kk}.stress(count,:) = stress';    
        [str,disp]    = exact_solution_plate_hole(vmesh{1,kk}.vertices(count,1:2),a,E,nu);
        vmesh{1,kk}.exactstress(count,:) = q0*str;
        vmesh{1,kk}.exactdisplacement(count,:) = q0*disp;
    end
end

% plot the solutions
figure;
title( 'Approximate Solutions', 'FontSize', 14', 'FontWeight', 'Bold') ;
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

% % plot Mises stress
% figure;
% title( 'Von Mises Stress (No Average)', 'FontSize', 14', 'FontWeight', 'Bold') ;
% for jj = 1:length(mesh)
%     face = polygon{1,jj}.trimesh;
%     maxnum = max(max(face));
%     vertices = vmesh{1,jj}.vertices(1:maxnum,:);
%     displacement = vmesh{1,jj}.displacement(1:maxnum,:);
%     stress = vmesh{1,jj}.stress(1:maxnum,:);
%     linmesh = kntcrv{1,jj}.linmesh;
%     p = patch('Faces',face, 'Vertices', vertices);
%     cdata = von_mises( stress(1:maxnum,:) );
%     hold on;
%     set(p,'FaceColor','interp','FaceVertexCData',cdata);
%     set(p,'EdgeColor','none');
%     for j = 1:size(linmesh,1)
%         hold on;
%         vv = vmesh{1,jj}.vertices(linmesh(j,:),:);
%         plot(vv(:,1), vv(:,2), 'k-');
%     end 
%     axis equal;
%     axis off;
% end
% hold off;
% colorbar;


% plot the exact solutions
figure;
title( 'Exact Solutions', 'FontSize', 14', 'FontWeight', 'Bold') ;
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
        displacement = vmesh{1,jj}.exactdisplacement(1:maxnum,:);
        stress = vmesh{1,jj}.exactstress(1:maxnum,:);
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

% plot the absolute error
figure;
title( 'Absolute Error', 'FontSize', 14', 'FontWeight', 'Bold') ;
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
        displacement = vmesh{1,jj}.displacement(1:maxnum,:) - vmesh{1,jj}.exactdisplacement(1:maxnum,:);
        stress = vmesh{1,jj}.stress(1:maxnum,:) - vmesh{1,jj}.exactstress(1:maxnum,:);
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