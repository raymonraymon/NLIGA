%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%  Hyperelastic isogeometric analysis of 2D curved beam under a prescribed 
%  large vertical displcaement and horizontal displacement.
%
%  The difference between the script 'ex_hyper_curvedbeam' and 
%  'ex_hyper_curvedbeam2' is the boundary condition.
%
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

clc;
clear all;

% Build geometrical model
plate = geo_curved_beam;

% Build iga mesh structure
mesh = build_iga_mesh( plate );

% Enforce displacement boundary conditions 
dbc = [];        % dbc = [node index, node dof, prescribed displacement]
tol = 1e-8;
side_node_1 = find( abs(mesh.coords(:,1)) < tol );   % left boundary
side_node_2 = find( abs(mesh.coords(:,2)) < tol );   % bottom boundary
% bottom boundary is clamped/fixed 
dbc = [dbc; side_node_2,   ones(length(side_node_2),1),   zeros(length(side_node_2),1)];
dbc = [dbc; side_node_2,   2*ones(length(side_node_2),1),   zeros(length(side_node_2),1)];
% left boundary is enforced by a prescribed vertical displacement and
% horizontal displacement with magnitude 18
dbc = [dbc; side_node_1,   ones(length(side_node_1),1),   18*ones(length(side_node_1),1)];
dbc = [dbc; side_node_1,   2*ones(length(side_node_1),1),   -18*ones(length(side_node_1),1)];

% Enforce traction boundary conditions 
tbc = [];        % dbc = [node index, node dof, prescribed displacement]

% Determine material properties
% Note that definition 'mat' is different for nonlinear materials, the first number
% is always the index to define material categories
% mat - [index, parameter1, parameter2,...]
% index - [10-20) - hyperelastic material
% for incompressible:
% mat - [ 10, K, A10]                  % Neo-Hookean
%     - [ 11, K, A10, A01]             % Mooney-Rivlin
%     - [ 12, K, A10, A20, A30]        % Yeoh
%     - [ 13, K, A10, A01, A20, A30]   % Bidtanerman
% K = 0-incompressible, >0 nearly-incompressible
mat=[10, 0, 80, 20 ];  % mat = [index K, A10 A01 A20 A30]

% Nonlinear analysis
eltype = 10;    % element type: 10 - plane strain element, 20- solid element
filename = 'curvedBeam2';
fname = get_output_file_name(filename);
fout = fopen(fname,'w'); 
nliga( eltype, plate, mesh, mat, dbc, tbc, fout );

vmesh = read_visual_mesh( fname );               % read the saved .msh file
outstr = von_mises(vmesh.stress{1,end});
midstr = outstr(2);
maxstr = max(outstr);
ux = vmesh.displacement{1,end}(2,1);
uy = vmesh.displacement{1,end}(2,2);
aa = [mesh.nElems, mesh.nCpts, midstr,maxstr,ux,uy];
flag = 7;
if mesh.dim == 2
    vmesh = plot_visual_mesh2d(flag, filename);
elseif mesh.dim == 3
    vmesh = plot_visual_mesh3d(flag, filename);
end

