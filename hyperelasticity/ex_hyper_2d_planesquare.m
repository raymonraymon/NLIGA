%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Hyperelastic tension example, Plane square
%
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

clc;
clear all;

% Build geometrical model
plate = geo_square([0,0],10);

% Build iga mesh structure
mesh = build_iga_mesh( plate );

% Enforce displacement boundary conditions 
dbc = [];        % dbc = [node index, node dof, prescribed displacement]
L = 10;  % side length
tol = 1e-8;
side_node_1 = find( abs(mesh.coords(:,1)) < tol );
side_node_2 = find( abs(mesh.coords(:,1) - L) < tol );
side_node_3 = find(abs(mesh.coords(:,2)) < tol);
dbc = [dbc; side_node_1,   ones(length(side_node_1),1),   zeros(length(side_node_1),1)];
dbc = [dbc; side_node_3,   2*ones(length(side_node_3),1),   zeros(length(side_node_3),1)];
dbc = [dbc; side_node_2,   1*ones(length(side_node_2),1),   10*ones(length(side_node_2),1)];

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
% mat=[11, 0, 0.1863 0.00979, -0.00186, 0.0000451 ];  % mat = [13 0 A10 A01 A20 A30]
mat=[11, 0, 0.1863 0.00979];  % mat = [10 0 A10 A01]
% Nonlinear analysis
eltype = 10;    % element type: 10 - plane strain element, 20- solid element

filename = 'planeSquare';
fname = get_output_file_name(filename);
fout = fopen(fname,'w'); 
nliga( eltype, plate, mesh, mat, dbc, tbc, fout );
fclose(fout);

flag = 1;
if mesh.dim == 2
    vmesh = plot_visual_mesh2d(flag, filename);
elseif mesh.dim == 3
    vmesh = plot_visual_mesh3d(flag, filename);
end
