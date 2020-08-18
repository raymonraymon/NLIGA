%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%  Elastoplastic isogeometric analysis of 2D circular a prescribed 
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
plate = geo_quarter_annulus(1,2);

% Build iga mesh structure
mesh = build_iga_mesh( plate );
% Build four edges
curve = extract_iga_boundary(mesh);
% Enforce traction boundary conditions
tbc = [];
noDofs = mesh.nCpts*2;
T = -4.5e7;          % Magnitude of force
d = 2;               % Direction of force
forcedcurve = curve.RightCurve;
tbc = force_boundary_condition ( noDofs, forcedcurve, T, d );

% Enforce displacement boundary conditions 
dbc = [];        % dbc = [node index, node dof, prescribed displacement]
tol = 1e-8;
side_node_1 = find( abs(mesh.coords(:,2)) < tol );   % right boundary
side_node_2 = find( abs(mesh.coords(:,1)) < tol );   % left boundary

% left boundary is clamped/fixed 
dbc = [dbc; side_node_2,   ones(length(side_node_2),1),   zeros(length(side_node_2),1)];
dbc = [dbc; side_node_2,   2*ones(length(side_node_2),1),   zeros(length(side_node_2),1)];

% right boundary a prescribed displacement is imposed in y-direction
% dbc = [dbc; side_node_1,   2*ones(length(side_node_1),1),   -5.6*10^-3*ones(length(side_node_1),1)];

% Determine material properties
% Note that definition 'mat' is different for nonlinear materials, the first number
% is always the index to define material categories
% index - [20-40) - elastoplastic material
index = 21;              % material flag
E = 200E9;               % Young's modulus
nu = 0.29;               % Poisso n ratio
sigmay0 = 200E6;         % initial yield stress
Et = 20E9;
H = (E*Et)/(E-Et);
% H = (20*2/(20-2))*1e10;  % plastic modulus
beta = 0.0;                % Bauschinger parameter
mat= [index, E, nu, sigmay0, H, beta];

% Nonlinear analysis
eltype = 10;    % element type: 10 - plane strain element, 20- solid element
filename = 'circular_pla';
fname = get_output_file_name(filename);
fout = fopen(fname,'w'); 
nliga( eltype, plate, mesh, mat, dbc, tbc, fout );
fclose(fout);

vmesh = read_results( fname, eltype, mesh.nCpts, mesh.nElems,  mesh.nElemCpts );
flag = 2;
plot_results( plate, vmesh, mesh, flag );

