function [ K, R ] = globalstiffness_hyper( eltype, geo, mesh, mat, u )
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%  Calculate global stiffness K and residual R in hyperelasticity
%  Input:
%    eltype - element types,  
%           10 - plane strain element
%           20 - solid element
%    geo - nurbs geometry
%    mesh - iga mesh structure
%    mat - material definition
%    u - current displacements
%  Output:
%    K - stiffness matric
%    R - residual vector
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

if eltype == 10             % plane element
    dof = 2;                % degree of freedom
elseif eltype == 20
    dof = 3; 
end
if mesh.dim == 2            % two dimensional
    gp_x = mesh.p+1;        % number of integration points in x-direction
    gp_y = mesh.q+1;        % number of integration points in y-direction
    [gp, wgt] = gauss_quadrature(gp_x, gp_y);   % calculate integration points and its weights
elseif mesh.dim == 3   % three dimensional
    gp_x = mesh.p+1;        % number of integration points in x-direction
    gp_y = mesh.q+1;        % number of integration points in y-direction
    gp_z = mesh.k+1;        % number of integration points in y-direction
    [gp, wgt] = gauss_quadrature(gp_x, gp_y, gp_z);   % calculate integration points and its weights
end
ndofs = dof * mesh.nCpts;      % total dofs

K = sparse(ndofs,ndofs);             % reserve stiffness matrix
R = zeros(ndofs,1);               % reserve residual matrix
count = 0;
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
    
    elDisp = u(sctrB);
    elDisp = reshape(elDisp, dof, nn);
    
    for ipt = 1:size(gp,1)            % loop over integration points
        count = count + 1;
        pt = gp(ipt,:);      % reference parametric coordinates for each integration point
        wt = wgt(ipt);       % weigths for each integration point
        gauPts = parameter_gauss_mapping( elDoma, pt );   % gauss integration mapping
        j1 = jacobian_gauss_mapping( elDoma );     % jacobian value for gauss mapping   
        [~,ders] = nurbs_derivatives( gauPts,geo, mesh );
        jmatrix = ders*elCpts(:,1:dof); 
        j2 = det(jmatrix);
        ders =  jmatrix \ ders;              
        fac = j1 *j2 * wt;      
        F = elDisp * ders' + eye(dof);      
%         [ stress, dtan ] = constitutive_relation( mesh.dim, mat, F );
        [ stress, dtan ] = material_hyperelasticity( mesh.dim, mat, F );
        if dof == 2
            BN = zeros(3,nnElem);
            BG = zeros(4,nnElem);
            for i = 1:nn
                BN(:,i*2-1:i*2) = [ F(1,1)*ders(1,i)     F(2,1)*ders(1,i);
                                    F(1,2)*ders(2,i)     F(2,2)*ders(2,i);
                                    F(1,1)*ders(2,i)+ F(1,2)*ders(1,i)  F(2,1)*ders(2,i)+F(2,2)*ders(1,i)];

                BG(:,i*2-1:i*2) = [ders(1,i)          0;
                                   ders(2,i)          0;
                                          0   ders(1,i);
                                          0   ders(2,i);];
            end
        elseif dof == 3
            BN = zeros(6,nn*3);
            BG = zeros(9,nn*3);
            for i = 1:nn
                BN(:,i*3-2:i*3) = [ F(1,1)*ders(1,i)     F(2,1)*ders(1,i)      F(3,1)*ders(1,i);
                                    F(1,2)*ders(2,i)     F(2,2)*ders(2,i)      F(3,2)*ders(2,i);
                                    F(1,3)*ders(3,i)     F(2,3)*ders(3,i)      F(3,3)*ders(3,i);
                                    F(1,1)*ders(2,i)+ F(1,2)*ders(1,i)  F(2,1)*ders(2,i)+F(2,2)*ders(1,i)   F(3,1)*ders(2,i)+F(3,2)*ders(1,i);
                                    F(1,2)*ders(3,i)+ F(1,3)*ders(2,i)  F(2,2)*ders(3,i)+F(2,3)*ders(2,i)   F(3,2)*ders(3,i)+F(3,3)*ders(2,i);
                                    F(1,3)*ders(1,i)+ F(1,1)*ders(3,i)  F(2,3)*ders(1,i)+F(2,1)*ders(3,i)   F(3,3)*ders(1,i)+F(3,1)*ders(3,i);];
                
                BG(:,i*3-2:i*3) = [ders(1,i)          0          0;
                                   ders(2,i)          0          0;
                                   ders(3,i)          0          0;
                                          0   ders(1,i)          0;
                                          0   ders(2,i)          0;
                                          0   ders(3,i)          0;
                                          0          0   ders(1,i);
                                          0          0   ders(2,i);
                                          0          0   ders(3,i)];
            end
        end
        R(sctrB) = R(sctrB) + fac*BN'*stress;  

        if dof == 2
            sigma = [stress(1) stress(3);
                     stress(3) stress(2)];
            stan = zeros(4);
            stan(1:2,1:2) = sigma;
            stan(3:4,3:4) = sigma;
        elseif dof == 3
            sigma = [stress(1) stress(4) stress(6);
                     stress(4) stress(2) stress(5);
                     stress(6) stress(5) stress(3)];
            stan=zeros(9);
            stan(1:3,1:3) = sigma;
            stan(4:6,4:6) = sigma;
            stan(7:9,7:9) = sigma;
        end
        EK = BN'*dtan*BN + BG'*stan*BG;              % element stiffness matrix
        K(sctrB,sctrB) = K(sctrB,sctrB) + fac*EK;    % assemble global stiffness                    
    end
    
end



end

