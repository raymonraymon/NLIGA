%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% This file is to visualize the results of different examples which you can
% choose after implementation.
% For three dimensional solid
%    flag - color map: 1-U1, 2-U2, 3-U3, 4, U magnitude, 5-S11, 6-S22
%                    7-S33, 8-S12, 9-S23, 10-S31, 11-mises
% For two dimensional surface
%    flag - color map: 1-U1, 2-U2, 3-U magnitude, 4-S11, 5-S22, 6-S12, 7-mises
% Please comment out the example you want and run the script.
% ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

clc;
clear;
dbstop if error


%% the curved beam example
filename = 'curvedBeam2';
% vmesh = plot_visual_mesh2d(7, filename);
steps = 1:5:21;
hold on;
plot_overlay2d( 7, filename, steps );
%-------------------------
% filename = 'curvedBeam';
% vmesh = plot_visual_mesh2d(7, filename);

%% the cube example
% filename = 'cube';
% vmesh = plot_visual_mesh3d(1, filename);

%% the square plane example
% filename = 'planeSquare';
% vmesh = plot_visual_mesh2d(1, filename);

%% the square plane example
% filename = 'plate_with_hole';
% vmesh = plot_visual_mesh2d(1, filename);

%% the overlay plot of the cube example
% filename = 'cube';
% steps = 1:7;
% plot_overlay3d( 4, filename, steps );
% view(3)

%% the overlay plot of the cube example




