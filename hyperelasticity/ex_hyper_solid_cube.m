%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%  Hyperelastic isogeometric analysis of a cube
%
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

clc;
clear all;

% Build geometrical model
cube = geo_cube3d([0,0,0],1);

% Build iga mesh structure
mesh = build_iga_mesh( cube );

% Enforce displacement boundary conditions 
dbc = [];        % dbc = [node index, node dof, prescribed displacement]
tol = 1e-8;
len = 1;           % the side length of the cube

bottom_face_nodes = find( abs(mesh.coords(:,3)) < tol );        % nodes on the bottom face
top_face_nodes = find( abs(mesh.coords(:,3)-len) < tol );       % nodes on the top face
left_face_nodes = find( abs(mesh.coords(:,1)) < tol );          % nodes on the left face
right_face_nodes = find( abs(mesh.coords(:,1)-len) < tol );     % nodes on the right face
front_face_nodes = find( abs(mesh.coords(:,2)) < tol );         % nodes on the front face
back_face_nodes = find( abs(mesh.coords(:,2)-len) < tol );      % nodes on the rear face
% the bottom face is simply supported 
dbc = [dbc; bottom_face_nodes,   3*ones(length(bottom_face_nodes),1),   zeros(length(bottom_face_nodes),1)];
% the left face is simply supported 
dbc = [dbc; left_face_nodes,   ones(length(left_face_nodes),1),   zeros(length(left_face_nodes),1)];
% the front face is simply supported 
dbc = [dbc; back_face_nodes,   2*ones(length(back_face_nodes),1),   zeros(length(back_face_nodes),1)];
% the right face is enforced by a prescribed displacement 
dbc = [dbc; right_face_nodes,   ones(length(right_face_nodes),1),   5*ones(length(right_face_nodes),1)];

% Enforce traction boundary conditions 
tbc = [];        % dbc = [node index, node dof, prescribed displacement]

% mat - [index, parameter1, parameter2,...]
% index - [10-20) - hyperelastic material
% for incompressible:
% mat - [ 10, K, A10]                  % Neo-Hookean
%     - [ 11, K, A10, A01]             % Mooney-Rivlin
%     - [ 12, K, A10, A20, A30]        % Yeoh
%     - [ 13, K, A10, A01, A20, A30]   % Bidtanerman
% K = 0-incompressible, >0 nearly-incompressible
% mat=[11, 1E5, 80, 20];  % mat = [index A10 A01 K]
mat=[11, 1E5, 0.1863, 0.00979]; 
% Nonlinear analysis
eltype = 20;    % element type: 10 - plane strain element, 20- solid element
filename = 'cube';
fname = get_output_file_name(filename);
fout = fopen(fname,'w'); 
nliga( eltype, cube, mesh, mat, dbc, tbc, fout );

flag = 1;
if mesh.dim == 2
    vmesh = plot_visual_mesh2d(flag, filename);
elseif mesh.dim == 3
    vmesh = plot_visual_mesh3d(flag, filename);
end

