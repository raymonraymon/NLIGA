%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%  Elastoplastic isogeometric analysis of a cube
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
side_node_1 = find( abs(mesh.coords(:,1)) < tol );   % left boundary
side_node_2 = find( abs(mesh.coords(:,2)) < tol );   % bottom boundary
bottom_face_nodes = find( abs(mesh.coords(:,3)) < tol );        % nodes on the bottom face
top_face_nodes = find( abs(mesh.coords(:,3)-len) < tol );       % nodes on the bottom face
left_face_nodes = find( abs(mesh.coords(:,1)) < tol );          % nodes on the bottom face
right_face_nodes = find( abs(mesh.coords(:,1)-len) < tol );     % nodes on the bottom face
front_face_nodes = find( abs(mesh.coords(:,2)) < tol );         % nodes on the bottom face
back_face_nodes = find( abs(mesh.coords(:,2)-len) < tol );      % nodes on the bottom face

% the bottom face is simply supported 
dbc = [dbc; bottom_face_nodes,   3*ones(length(bottom_face_nodes),1),   zeros(length(bottom_face_nodes),1)];
dbc = [dbc; bottom_face_nodes,   2*ones(length(bottom_face_nodes),1),   zeros(length(bottom_face_nodes),1)];
dbc = [dbc; bottom_face_nodes,   1*ones(length(bottom_face_nodes),1),   zeros(length(bottom_face_nodes),1)];

% the front face is enforced by a prescribed displacement 
dbc = [dbc; top_face_nodes,   3*ones(length(top_face_nodes),1),   2.6*10^-3*ones(length(top_face_nodes),1)];


% Enforce traction boundary conditions 
tbc = [];        % dbc = [node index, node dof, prescribed displacement]

% index - [20-40) - elastoplastic material
index = 21;              % material flag
E = 200E9;               % Young's modulus
nu = 0.29;               % Poisso n ratio
sigmay0 = 200E6;         % initial yield stress
H = (20*2/(20-2))*1e10;  % plastic modulus
beta = 0;                % Bauschinger parameter
mat= [index, E, nu, sigmay0, H, beta];

% Nonlinear analysis
eltype = 20;    % element type: 10 - plane strain element, 20- solid element
filename = 'cube_pla';
fname = get_output_file_name(filename);
fout = fopen(fname,'w'); 
nliga( eltype, cube, mesh, mat, dbc, tbc, fout );
fclose(fout);
