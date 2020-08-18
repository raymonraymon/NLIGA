function tbc = force_boundary_condition ( ndofs, fb, t, d )
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%  impose force boundary condition on the given boundary of a two dimensional mesh
%  Input:
%    ndofs - dimensions of the stiffness matrix
%    fb    - the boundary to be imposed force boundary condition
%    t     - magnitude of the force
%    b     - direction of the force
%  Output:
%    tbc   - force boundary condition
%       tbc(:,1) - global number of control points
%       tbc(:,2) - force direction, 1 - x direction, 2 - y direction
%       tbc(:,3) - force magnitude on the corresponding control points
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

[Q, W] = gauss_quadrature(fb.p+1);
pts = zeros(fb.p+1,2);

f = zeros(ndofs,1);
tbc = zeros(ndofs,3);

for e = 1:fb.nElems
    elDoma = fb.elDoma(e,:);    
    sctr = fb.elNodeCnt_global(e,:); 
    sctrl = fb.elNodeCnt_local(e,:);
    coords_pts = fb.coords(sctrl,:);
        
    pts(:,1) = coords_pts(:,1)./coords_pts(:,4);
    pts(:,2) = coords_pts(:,2)./coords_pts(:,4);
    
    if(d == 1)
        sctrw =  2*sctr-1;
    else
        sctrw =  2*sctr;
    end
    
    % loop over Gauss points
    for gp = 1:size(W,1)
        pt =  Q(gp,:);
        wt = W(gp);
        
        xi = parameter_gauss_mapping(elDoma,pt);        
        j1 = jacobian_gauss_mapping( elDoma );
        [R, dRdxi] = nrbNurbs1DBasisDerivs (xi,fb.p,fb.knot,fb.coords(:,4));
        j2 = vecmag(dRdxi*pts);
        
        f(sctrw) = f(sctrw) + R'*t*j1*j2*wt;
    end

end

for i = 1:ndofs
    tbc(i,1) = ceil(i/2);
    tbc(i,2) = mod(i+1,2)+1;
    tbc(i,3) = f(i);
end

del = (tbc(:, 3)==0);
tbc(del,:) = [];

end