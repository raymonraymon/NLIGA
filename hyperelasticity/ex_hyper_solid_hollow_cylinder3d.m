%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%  Hyperelastic isogeometric analysis of a hollow cylinder
%
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

clc;
clear all;

% Build geometrical model
cube = geo_quarter_hollow_cylinder3d;

% Build iga mesh structure
mesh = build_iga_mesh( cube );

% Enforce displacement boundary conditions 
dbc = [];        % dbc = [node index, node dof, prescribed displacement]
tol = 1e-8;
L = 15;           % the side length of the cube
side_node_1 = find( abs(mesh.coords(:,2)) < tol );     % left face
side_node_2 = find( abs(mesh.coords(:,2)-L) < tol );   % right face
side_node_3 =  find( abs(mesh.coords(:,1)) < tol );    % x = 0
side_node_4 =  find( abs(mesh.coords(:,3)) < tol );
% hold on;
% sctr  = side_node_4;
% plot3(mesh.coords(sctr,1),mesh.coords(sctr,2),mesh.coords(sctr,3),'o')
% the left face is simply supported 
dbc = [dbc; side_node_1,   2*ones(length(side_node_1),1),   zeros(length(side_node_1),1)];
% the right face is enforced by a prescribed displacement 
dbc = [dbc; side_node_2,   2*ones(length(side_node_2),1),   75*ones(length(side_node_2),1)];
dbc = [dbc; side_node_3,   ones(length(side_node_3),1),   zeros(length(side_node_3),1)];
dbc = [dbc; side_node_4,   3*ones(length(side_node_4),1),   zeros(length(side_node_4),1)];
% Enforce traction boundary conditions 
tbc = [];        % tbc = [node index, node dof, prescribed force]

% mat - [index, parameter1, parameter2,...]
% index - [10-20) - hyperelastic material
% for incompressible:
% mat - [ 10, K, A10]                  % Neo-Hookean
%     - [ 11, K, A10, A01]             % Mooney-Rivlin
%     - [ 12, K, A10, A20, A30]        % Yeoh
%     - [ 13, K, A10, A01, A20, A30]   % Bidtanerman
% K = 0-incompressible, >0 nearly-incompressible
% mat=[11, 1E5, 80, 20];  % mat = [index A10 A01 K]
mat = [11,1e5,0.5516,0.1379];
% mat=[11, 1e5, 0.1863 0.00979];
% Nonlinear analysis
eltype = 20;    % element type: 10 - plane strain element, 20- solid element
filename = 'quarder_hollow_cylinder3d';
fname = get_output_file_name(filename);
fout = fopen(fname,'w');
nliga( eltype, cube, mesh, mat, dbc, tbc, fout );

flag = 2;
if mesh.dim == 2
    vmesh = plot_visual_mesh2d(flag, filename);
elseif mesh.dim == 3
    vmesh = plot_visual_mesh3d(flag, filename);
end

