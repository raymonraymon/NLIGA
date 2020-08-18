%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%  Poisson problem on a L-shaped plane model
%  The model is built with two conforming NURBS patches.
%
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

clc; 
clear all;

dof = 1;             % degree of freedom
% define the geometry
geo = geo_Lshape_2patch(1,[-1,-1]);

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
mesh{1,2}.nodeNum(:,1) = mesh{1,1}.nodeNum(:,end);
mesh{1,2}.globCoords = zeros(mesh{1,2}.nCptsV*(mesh{1,2}.nCptsU-1),4);
for j = 1:mesh{1,2}.nCptsV
    for i = 2:mesh{1,2}.nCptsU
        aa = (j-1)*(mesh{1,2}.nCptsU)+i;
        bb = (j-1)*(mesh{1,2}.nCptsU-1)+i-1;
        mesh{1,2}.nodeNum(j,i) = bb + mesh{1,1}.nodeNum(end,end);
        mesh{1,2}.globCoords(bb,:) = mesh{1,2}.coords(aa,:); 
    end
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
bcNodes1 = mesh{1,1}.nodeNum(1,:)';
bcNodes2 = mesh{1,2}.nodeNum(1,:)';
bcNodes3 = mesh{1,2}.nodeNum(:,end);
bcNodes4 = mesh{1,2}.nodeNum(end,:)';
bcNodes5 = mesh{1,1}.nodeNum(end,:)';
bcNodes6 = mesh{1,1}.nodeNum(:,1);
bcNodes = unique([bcNodes1; bcNodes2; bcNodes3; bcNodes4; bcNodes5; bcNodes6]);

dbc = [];    % dbc = [node index, direction, prescribed displacement]
dbc = [dbc; bcNodes,  ones(size(bcNodes)),   zeros(size(bcNodes))];

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
    num1 = 100;  % if you want to obtain a more smooth results, please provide large numbers for num1 and num2    
    num2 = 100;
    polygon{1,kk} = build_visual_mesh_suf( num1, num2 );        % build visualized mesh
    kntcrv{1,kk} = build_visual_knotcurve_suf( mesh{1,kk}.uKnots, mesh{1,kk}.vKnots, num1+1 ); % build visualized knot curves
    numpts = (num1+1)*(num2+1) + size(kntcrv{1,kk}.linpts,1);
    vmesh{1,kk}.vertices = zeros(numpts,2);
    vmesh{1,kk}.displacement = zeros(numpts,1);
    vmesh{1,kk}.exactdisplacement = zeros(numpts,1);
    vmesh{1,kk}.stress = zeros(numpts,2);
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
        B = zeros(2,nnElem); 
        B(1,1:nnElem) = dR_dx(1,:);  
        B(2,1:nnElem) = dR_dx(2,:); 
        strain = B * u(sctrB);
        stress = strain;
        edsp   = u(sctrB);
        edsp   = reshape(edsp, dof, nn);
        vmesh{1,kk}.displacement(i,1) = R * edsp';
        vmesh{1,kk}.vertices(i,1:2) = R*exyz(:,1:2);
        vmesh{1,kk}.stress(i,:) = stress';
        x = vmesh{1,kk}.vertices(i,1);
        y = vmesh{1,kk}.vertices(i,2);
        ux = sin(x*pi)*sin(y*pi);
        vmesh{1,kk}.exactdisplacement(i) = ux;
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
        B = zeros(2,nnElem); 
        B(1,1:nnElem) = dR_dx(1,:);  
        B(2,1:nnElem) = dR_dx(2,:);  
        strain = B * u(sctrB);
        stress = strain;
        edsp   = u(sctrB);
        edsp   = reshape(edsp, dof, nn);
        vmesh{1,kk}.displacement(count) = R * edsp';
        vmesh{1,kk}.vertices(count,1:2) = R*exyz(:,1:2); 
        vmesh{1,kk}.stress(count,:) = stress';    
        x = vmesh{1,kk}.vertices(i,1);
        y = vmesh{1,kk}.vertices(i,2);
        ux = sin(x*pi)*sin(y*pi);
        vmesh{1,kk}.exactdisplacement(i) = ux;
    end
end

% plot the approximate solutions
figure;
for jj = 1:length(mesh)
    face = polygon{1,jj}.trimesh;
    maxnum = max(max(face));
    vertices = vmesh{1,jj}.vertices(1:maxnum,:);
    displacement = vmesh{1,jj}.displacement(1:maxnum,:);
    stress = vmesh{1,jj}.stress(1:maxnum,:);
    linmesh = kntcrv{1,jj}.linmesh;
    p = patch('Faces',face, 'Vertices', vertices);
    cdata = displacement(1:maxnum);  title('Approximated Solution');
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

% plot the exact solutions
figure;
for jj = 1:length(mesh)
    face = polygon{1,jj}.trimesh;
    maxnum = max(max(face));
    vertices = vmesh{1,jj}.vertices(1:maxnum,:);
    displacement = vmesh{1,jj}.exactdisplacement(1:maxnum,:);
    stress = vmesh{1,jj}.stress(1:maxnum,:);
    linmesh = kntcrv{1,jj}.linmesh;
    p = patch('Faces',face, 'Vertices', vertices);
    cdata = displacement(1:maxnum);  title('Exact Solution');
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

% plot the absolute error
figure;
for jj = 1:length(mesh)
    face = polygon{1,jj}.trimesh;
    maxnum = max(max(face));
    vertices = vmesh{1,jj}.vertices(1:maxnum,:);
    displacement = vmesh{1,jj}.displacement(1:maxnum,:) - vmesh{1,jj}.exactdisplacement(1:maxnum,:);
    stress = vmesh{1,jj}.stress(1:maxnum,:);
    linmesh = kntcrv{1,jj}.linmesh;
    p = patch('Faces',face, 'Vertices', vertices);
    cdata = displacement(1:maxnum);  title('Absolute Error');
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

% compute the error
err1 = 0;
err2 = 0;
% use gaussian integration rule
gp_x = mesh{1,1}.p+1;           % number of integration points in x-direction
gp_y = mesh{1,1}.q+1;           % number of integration points in y-direction
[gp, wgt] = gauss_quadrature(gp_x, gp_y);   % calculate integration points and its weights
for jj = 1:length(mesh)
    for e = 1:mesh{1,jj}.nElems                 % loop over elements
        sctr   = mesh{1,jj}.gloElNodeCnt(e,:);  % element control points index
        elDoma = mesh{1,jj}.elDoma(e,:);        % element parametric domain
        elCpts = globCoords(sctr,:);     % coordinates of element control points
        for ipt = 1:size(gp,1)        % loop over integration points
            pt      = gp(ipt,:);      % reference parametric coordinates for each integration point
            wt      = wgt(ipt);       % weigths for each integration point
            gauPts  = parameter_gauss_mapping( elDoma, pt );       % gauss integration mapping
            [R,ders] = nurbs_derivatives( gauPts, geo{1,jj}, mesh{1,jj} );
            u_h = R * u(sctr);
            coords = R * elCpts(:,1:2);
            y = coords(2);
            x = coords(1);
            u_ex = sin(x*pi)*sin(y*pi);
            err1 = err1 + (u_h - u_ex)^2;
            err2 = err2 + u_ex^2;
        end 
    end
end

err1 = sqrt(err1);
err2 = sqrt(err2);
err2 = err1/err2;

nelems = 0;
ncpts = 0;
for jj = 1:length(mesh)
    nelems = nelems + mesh{1,jj}.nElems;
    ncpts = ncpts + mesh{1,jj}.nCpts;
end
format long
err = [err1 err2]
dof = [nelems globNodeNum]