%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Hyperelastic example, two dimensional plate with hole 
%
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

clc;
clear all;

% Build geometrical model
L = 82.5;  % side length
R = 6.35;  % radius
platehole = geo_plate_with_hole(L, R);
% Build iga mesh structure
mesh = build_iga_mesh( platehole );

% Enforce displacement boundary conditions 
dbc = [];        % dbc = [node index, node dof, prescribed displacement]

tol = 1e-8;
side_node_1 = find( abs(mesh.coords(:,2)) < tol );      % bottom edge
side_node_2 = find( abs(mesh.coords(:,1)) < tol );      % left edge
side_node_3 = find( abs(mesh.coords(:,2) - L) < tol );  % top edge
side_node_4 = find( abs(mesh.coords(:,1) - L) < tol );  % right edge
dbc = [dbc; side_node_1, 2*ones(length(side_node_1),1),   zeros(length(side_node_1),1)];   % bottome edge is simply supported
dbc = [dbc; side_node_2,   ones(length(side_node_2),1),   zeros(length(side_node_2),1)];   % left edge is simply supported
dbc = [dbc; side_node_3, 2*ones(length(side_node_3),1),   zeros(length(side_node_3),1)];   % top edge is simply supported
dbc = [dbc; side_node_4,   ones(length(side_node_4),1), 50.8*ones(length(side_node_4), 1)];  % right edge is enforced by a prescribed displacement

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
mat=[11, 0, 0.1863 0.00979, -0.00186, 0.0000451 ];  % mat = [13 0 A10 A01 A20 A30]
% mat=[12, 0, 0.1863 -0.00186, 0.0000451];  % mat = [10 0 A10 A01]
% Nonlinear analysis
eltype = 10;    % element type: 10 - plane strain element, 20- solid element
filename = 'plate_with_hole';
fname = get_output_file_name(filename);
fout = fopen(fname,'w'); 
nliga( eltype, platehole, mesh, mat, dbc, tbc, fout );

% flag - color map: 1-U1, 2-U2, 3-U magnitude, 4-S11, 5-S22, 6-S12, 7-mises
flag = 2;
if mesh.dim == 2
    vmesh = plot_visual_mesh2d(flag, filename);
elseif mesh.dim == 3
    vmesh = plot_visual_mesh3d(flag, filename);
end


totalElem = mesh.nElems
totalNode = mesh.nCpts

s11 = vmesh.stress{1,end}(1,1)
uy = vmesh.displacement{1,end}(1,2)
