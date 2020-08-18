%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%  Hyperelastic isogeometric analysis of 3d square with a hole
%
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

clc;
clear all;

% Build geometrical model
cube = geo_plate_with_hole3d;

% Build iga mesh structure
mesh = build_iga_mesh( cube );

% Enforce displacement boundary conditions 
dbc = [];        % dbc = [node index, node dof, prescribed displacement]
L = 82.5;  % side length
R = 6.35;  % radius
h = 10; % height
tol = 1e-8;
bottom_face_nodes = find( abs(mesh.coords(:,3)) < tol );        % nodes on the bottom face
top_face_nodes = find( abs(mesh.coords(:,3)-h) < tol );       % nodes on the top face
left_face_nodes = find( abs(mesh.coords(:,1)) < tol );          % nodes on the left face
right_face_nodes = find( abs(mesh.coords(:,1)-L) < tol );     % nodes on the right face
front_face_nodes = find( abs(mesh.coords(:,2)) < tol );         % nodes on the front face
rear_face_nodes = find( abs(mesh.coords(:,2)-L) < tol );      % nodes on the rear face
% the bottom face is simply supported 
dbc = [dbc; bottom_face_nodes,   3*ones(length(bottom_face_nodes),1),   zeros(length(bottom_face_nodes),1)];
% the bottom face is simply supported 
% dbc = [dbc; top_face_nodes,   3*ones(length(top_face_nodes),1),   zeros(length(top_face_nodes),1)];
% the left face is simply supported 
dbc = [dbc; left_face_nodes,   ones(length(left_face_nodes),1),   zeros(length(left_face_nodes),1)];
% the front face is simply supported 
dbc = [dbc; front_face_nodes,   2*ones(length(front_face_nodes),1),   zeros(length(front_face_nodes),1)];
% the rear face is simply supported 
dbc = [dbc; rear_face_nodes,   2*ones(length(rear_face_nodes),1),   zeros(length(rear_face_nodes),1)];
% the right face is enforced by a prescribed displacement 
dbc = [dbc; right_face_nodes,   ones(length(right_face_nodes),1),   10*ones(length(right_face_nodes),1)];

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
mat=[11, 1E5, 80, 20];  % mat = [index A10 A01 K]

% Nonlinear analysis
eltype = 20;    % element type: 10 - plane strain element, 20- solid element
filename = 'plate_with_hole3d';
fname = get_output_file_name(filename);
fout = fopen(fname,'w');
nliga( eltype, cube, mesh, mat, dbc, tbc, fout );

flag = 1;
if mesh.dim == 2
    vmesh = plot_visual_mesh2d(flag, filename);
elseif mesh.dim == 3
    vmesh = plot_visual_mesh3d(flag, filename);
end

