function nliga( eltype, geo, mesh, mat, dbc, tbc, fout )
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Main frame for nonlinear isogeometric analysis
%  Input:
%    eltype - element types,  
%           10 - plane strain element
%           20 - solid element
%    geo - nurbs geometry
%    mesh - iga mesh structure
%    mat - material definition
%    dbc - displacements boundary conditions
%    tbc - tractions boundary conditions 
%    fout - output figure handle for visualization
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

% initialization
if eltype == 10      % plane element
    dof = 2;    % degree of freedom
    egp = (mesh.p+1)*(mesh.q+1); % gauss point in each element
elseif eltype == 20  % solid element
    dof = 3;    % degree of freedom
    egp = (mesh.p+1)*(mesh.q+1)*(mesh.k+1); % gauss point in each element
end
ndofs = dof * mesh.nCpts;   % total dofs
ngp   = egp*mesh.nElems;    % total integration points
scatdbc = [];
scattbc = [];
if ~isempty(dbc)
    scatdbc = dof * (dbc(:,1)-1) + dbc(:,2);   % scatter dbc
end
if ~isempty(tbc) 
    scattbc = dof * (tbc(:,1)-1) + tbc(:,2);   % scatter tbc
end
tol = 1e-6;           % tolerance of convergence
maxit = 20;           % maximum iterative steps 
reit = 0;             % reduction index
maxreit = 6;          % maximum load/displacement step reduction 
ndbc = size(dbc,1);   % number of displacement constrained nodes
ntbc = size(tbc,1);   % number of displacement constrained nodes
u = zeros(ndofs,1);   % nodal displacements

cu = zeros(ndofs,1);  % converged nodal displacements
step = 0;             % load/displacement step index
curtime = 0;          % current time
timeInterval = 0.25;  % initial time interval
cnit = [];            % record the iterative steps

if ( mat(1) >= 20 && mat(1) < 40 )        % belongs to plastic materials
    D = plastic_init_setting(mat, ngp, eltype);
end
if ( mat(1) >= 10 && mat(1) < 20 )            % belongs to hyperelastic materials
    % % output initial undeformed geometries
    if mesh.dim == 2
        output_visual_mesh2d( fout, mat, geo, mesh, u, step, curtime );
    elseif mesh.dim == 3
        output_visual_mesh3d( fout, mat, geo, mesh, u, step, curtime );
    end
end

while curtime ~= 1    % get to the end
    curtime = curtime + timeInterval;
    if curtime > 1
        timeInterval = 1 - curtime + timeInterval;
        curtime = 1; 
    end
    err = 1e6;        % record iterative error
    perr = err;
    nit = 0;          % record iterative times
    fprintf(1,'\n \t time    time step    iter \t  residual \n');
    iu = zeros(ndofs,1);   % record iterative displacements
    while (err > tol) && ( nit <= maxit)
        nit = nit+1;
        if ( mat(1) >= 10 && mat(1) < 20 )            % belongs to hyperelastic materials
            [ k, r ] = globalstiffness_hyper( eltype, geo, mesh, mat, u );
        elseif ( mat(1) >= 20 && mat(1) < 40 )        % belongs to plastic materials
            [ k, r ] = globalstiffness_plastic( D, eltype, geo, mesh, mat, iu );
        end
        f = zeros(ndofs,1);         % define external force
        if ntbc~=0                  % enforce traction conditions  
            f(scattbc) = tbc(:,3);
        end
           
        if ndbc~=0                  % enforce displacement conditions
            k(scatdbc,:) = zeros(ndbc, ndofs);
            k(scatdbc,scatdbc) = eye(ndbc);
            f(scatdbc,:) = 0;
            if nit == 1
                f(scatdbc,:) = dbc(:,3); 
            end          
        end
        b = curtime*f - r;          % define right side of the governing equation
        if ndbc~=0  
            b(scatdbc) = curtime*dbc(:,3) - u(scatdbc);   
        end
        du = k\b;                   % solve equation
        
        alldof = 1:ndofs;
        freedof = setdiff(alldof, scatdbc);    % nodes without displacement constraint
        u = u + du;                 % update displacement 
        iu = iu + du;               % update increment displacement
        if nit > 1                  % compute iterative error
            num = b(freedof)' * b(freedof);
            denom = 1+f(freedof)' * f(freedof);
            err = num/denom; 
        end
        % output current time step and iterative error
        fprintf(1,'%10.5f %10.3e %5d %14.5e \n',curtime,timeInterval,nit,err); 
        if err/perr > 1E3 && nit > 2
            nit = maxit+1;   % lf solution diverge extremely, go to next iteration
        else
            perr = err;
        end
    end
    
    if  nit <= maxit               % converged 
        reit = 0;                  % reset reduction index
        step = step + 1;           % increase converged steps by 1
        cu = u;                    % update converged displacement
        cnit = [cnit, nit];
        if length(cnit) >=2 && all(cnit(end-1:end) <= 5)
            timeInterval = timeInterval*1.5;  % increase the increment by times 1.5 
        end
        if ( mat(1) >= 10 && mat(1) < 20 )            % belongs to hyperelastic materials
            % output visualized mesh file with 'filename'
            if mesh.dim == 2
                output_visual_mesh2d( fout, mat, geo, mesh, u, step, curtime );
            elseif mesh.dim == 3
                output_visual_mesh3d( fout, mat, geo, mesh, u, step, curtime );
            end
        elseif ( mat(1) >= 20 && mat(1) < 40 )        % belongs to plastic materials           
            output_plastic( D, eltype, fout, mat, geo, mesh, iu, u, step, curtime);          
        end
        
    else                           % not converged
        if reit <= maxreit         % refine time interval and continue iterating
            curtime = curtime - timeInterval;   % recover current time step
            timeInterval = timeInterval/4;      % refine time interval
            reit = reit+1;         % increase reduction index
            u = cu;                % recover current displacement from last converged displacement
        else
            return;                % stop analysis
        end
    end
end
end

