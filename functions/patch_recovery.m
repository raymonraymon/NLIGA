function cstrs = patch_recovery( gstrs, mesh, geo)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% recover the stress at Gauss points to the control points
% Input: 
%    gstrs   - stress at Gauss points
%    mesh    - mesh struture
%    nelems  - number of elements
%    geo     - geometry structure
% Output:
%    cstrs   - stress at the control points
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

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
N = zeros(size(gstrs,1), mesh.nCpts ); 
count = 0;
for e = 1:mesh.nElems                % loop over elements
    sctr = mesh.elNodeCnt(e,:);       % element control points index
    elDoma = mesh.elDoma(e,:);        % element parametric domain
    for ipt = 1:size(gp,1)            % loop over integration points
        count = count + 1;
        pt = gp(ipt,:);      % reference parametric coordinates for each integration point
        gauPts = parameter_gauss_mapping( elDoma, pt );    % gauss integration mapping
        [R,~] = nurbs_derivatives( gauPts, geo, mesh ); %
        N(count,sctr) = R;
    end
end

cstrs = (N'*N)\(N'*gstrs);


end

