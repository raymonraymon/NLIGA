function output_plastic( D, eltype, fout, mat, geo, mesh, iu, u, step, curtime)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%   Output the displacement of the control points and the stresses on the
%    gauss points.
%  Input:
%    fout - output file handle
%    mat - material properties
%    geo - geometry
%    mesh - iga mesh structure
%    iu - displacement in current step
%    u - current displacement
%    step - iterative step index
%    currentime - current time
%  Output:
%    the file will be saved into fout
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

global SIGMA XQ;

if eltype == 10             % plane element
    dof = 2;                % degree of freedom
elseif eltype == 20
    dof = 3; 
end
if mesh.dim == 2            % two dimensional
    gp_x = mesh.p+1;        % number of integration points in x-direction
    gp_y = mesh.q+1;        % number of integration points in y-direction
    [gp, ~] = gauss_quadrature(gp_x, gp_y);   % calculate integration points and its weights
elseif mesh.dim == 3   % three dimensional
    gp_x = mesh.p+1;        % number of integration points in x-direction
    gp_y = mesh.q+1;        % number of integration points in y-direction
    gp_z = mesh.k+1;        % number of integration points in y-direction
    [gp, ~] = gauss_quadrature(gp_x, gp_y, gp_z);   % calculate integration points and its weights
end
                              
count = 0;                            % count for each integration point
for e = 1:mesh.nElems                % loop over elements
    sctr = mesh.elNodeCnt(e,:);       % element control points index
    elDoma = mesh.elDoma(e,:);        % element parametric domain
    elCpts = mesh.coords(sctr,:);     % coordinates of element control points
    nn = length(sctr);                % number of control points for each element
    nnElem = nn*dof;                  % dof for each element
    sctrB = zeros(1, nnElem); 
    for i = 1:dof
        sctrB(i:dof:nnElem) = dof*(sctr-1) + i;  % displacement in i-th direction
    end
    
    elDisp = iu(sctrB);
    elDisp = reshape(elDisp, dof, nn);
    
    for ipt = 1:size(gp,1)            % loop over integration points
        count = count + 1;
        pt = gp(ipt,:);      % reference parametric coordinates for each integration point
        gauPts = parameter_gauss_mapping( elDoma, pt );   % gauss integration mapping  
        [~,ders] = nurbs_derivatives( gauPts, geo, mesh );
        jmatrix = ders*elCpts(:,1:dof); 
        ders =  jmatrix \ ders;              
        deps  = elDisp*ders';
        
        [ stress, alpha, ep, ~ ] = material_plasticity( D, dof, mat, deps, count );
        
        if dof == 2
            SIGMA(1:4, count) = stress;
            XQ(:,count) = [alpha;ep];
        elseif dof == 3
            SIGMA(1:6, count) = stress;
            XQ(:,count) = [alpha;ep];
        end
    end 
end

npts = mesh.nCpts;
nelms = mesh.nElems;
if eltype == 10      % plane element
    dof = 2;    % degree of freedom
    egp = (mesh.p+1)*(mesh.q+1); % gauss point in each element
elseif eltype == 20  % solid element
    dof = 3;    % degree of freedom
    egp = (mesh.p+1)*(mesh.q+1)*(mesh.k+1); % gauss point in each element
end

fprintf(fout,'TIME = %11.3e STEP = %d\r\nNodal Displacements\r\n',curtime,step);
if eltype == 10,    fprintf(fout,'Node           U1          U2\r\n');
elseif eltype == 20, fprintf(fout,'Node           U1          U2          U3\r\n');
end

for i = 1:npts
    fprintf(fout, '%5d ',i );
    for j = 1:dof
        fprintf(fout,'%11.3e', u(dof*(i-1)+j));
    end
    fprintf(fout, '\r\n');
end

fprintf(fout,'Element Stress\r\n');
if dof == 2,    fprintf(fout,'        S11         S22         S12\r\n');
elseif dof == 3, fprintf(fout,'        S11         S22         S33         S12         S23         S13\r\n');
end
for i=1:nelms
    fprintf(fout,'Element %5d\r\n',i);
    j = (i-1)*egp;
    if eltype == 10      % plane element
        fprintf(fout,'%11.3e %11.3e %11.3e\r\n',SIGMA(1:3,j+1:j+egp));
    elseif eltype == 20  % solid element
        fprintf(fout,'%11.3e %11.3e %11.3e %11.3e %11.3e %11.3e\r\n',SIGMA(1:6,j+1:j+egp));
    end

end

fprintf(fout,'\r\n\r\n');

end