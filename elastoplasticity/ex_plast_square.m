%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%  Elastoplastic isogeometric analysis of 2D square a prescribed 
%  vertical displcaement.
%
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

clc;
clear all;

% Build geometrical model
plate =  geo_square( [0,0], 1 );

% Build iga mesh structure
mesh = build_iga_mesh( plate );

% Build four edges
curve = extract_iga_boundary(mesh);

% % Enforce traction boundary conditions
tbc = [];
% noDofs = mesh.nCpts*2;
% T = 320e6;          %Magnitude of force
% d = 1;            %Direction of force
% forcedcurve = curve.RightCurve;
% tbc = force_boundary_condition ( noDofs, forcedcurve, T, d );


% Enforce displacement boundary conditions 
dbc = [];        % dbc = [node index, node dof, prescribed displacement]
tol = 1e-8;
side_node_1 = find( abs(mesh.coords(:,2)-1) < tol );   % top boundary
side_node_2 = find( abs(mesh.coords(:,2)) < tol );     % bottom boundary

side_node_3 = find( abs(mesh.coords(:,1)-1) < tol );   % right boundary
side_node_4 = find( abs(mesh.coords(:,1)) < tol );     % left boundary
% 

% left boundary is clamped/fixed 
dbc = [dbc; side_node_4,   ones(length(side_node_4),1),   zeros(length(side_node_4),1)];
dbc = [dbc; side_node_4,   2*ones(length(side_node_4),1),   zeros(length(side_node_4),1)];
% right boundary is subjected to a precribed displacement
dbc = [dbc; side_node_3,   1*ones(length(side_node_3),1),   0.0056*ones(length(side_node_3),1)];

% Determine material properties
% Note that definition 'mat' is different for nonlinear materials, the first number
% is always the index to define material categories
% index - [20-40) - elastoplastic material

% E = 200E9;
% nu = 0.29;
% LAMDA = E*nu/(1+nu)/(1-2*nu);
% MU = E/2/(1+nu);
% Y0 = 200E6;
% H = (20*2/(20-2))*1e10;
% BETA = 0;
% index = 21;
% mat=[index LAMDA MU BETA H Y0 E nu];
index = 21;
E = 200E9;
nu = 0.29;
sigmay0 = 200E6;
H = (20*2/(20-2))*1e10;
beta = 0;
mat= [index, E, nu, sigmay0, H, beta];

% Nonlinear analysis
eltype = 10;    % element type: 10 - plane strain element, 20- solid element
filename = 'square_pla';
fname = fileparts(mfilename('fullpath')); % get current path
index_dir = strfind(fname,'\');
str_temp = fname(1:index_dir(end));
fname = [str_temp,'output\', filename, '.msh'];

fout = fopen(fname,'w'); 
nliga( eltype, plate, mesh, mat, dbc, tbc, fout );
fclose(fout);
vmesh = read_results( fname, eltype, mesh.nCpts, mesh.nElems,  mesh.nElemCpts );
fclose(fout);
flag = 3;
plot_results( plate, vmesh, mesh, flag );
% vmesh = plot_visual_mesh(flag, filename);

