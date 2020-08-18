function pinchedCylinder = geo_pinched_cylinder_solid(R, L, t)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Build a three-dimensional pinched cylinder SOLID model
% Note that only one-eight of the classical pinched cylinder model is built
% 'geo_pinched_cylinder_shell' is the middle surface of this solid model.
% Input:
%   R - radius 
%   L - length
%   t - thickness
% 
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

if nargin ~= 3    % default parameters
    R = 300;
    L = 600;
    t = 3;
end

theta = 45;
w = cos(theta/180*pi);
r1 = R-t/2;
r2 = R+t/2;
% define control points and knots vector
coefs = zeros(4,3,2,2);
coefs(:,:,1,1) = [0, 0, r1, 1; r1*w, 0, r1*w, w; r1, 0, 0, 1; ]';
coefs(:,:,2,1) = [0, L/2, r1, 1; r1*w, L/2*w, r1*w, w; r1, L/2, 0, 1;]';
coefs(:,:,1,2) = [0, 0, r2, 1; r2*w, 0, r2*w, w; r2, 0, 0, 1; ]';
coefs(:,:,2,2) = [0, L/2, r2, 1; r2*w, L/2*w, r2*w, w; r2, L/2, 0, 1;]';
knots{1} = [0 0 0 1 1 1];
knots{2} = [0 0 1 1];
knots{3} = [0 0 1 1];

% build nurbs solid by using control points and knots vector
pinchedCylinder = nrbmak(coefs,knots);

% degeree elevate
pinchedCylinder = nrbdegelev(pinchedCylinder,[1,2,2]);

% insert knots
RefinementX = 5;    % the number of knots inseted in u direction 
RefinementY = RefinementX;    % the number of knots inseted in v direction
RefinementZ = 1;    % the number of knots inseted in w direction
iuknots = 1/(RefinementX+1):1/(RefinementX+1):RefinementX/(RefinementX+1);
ivknots = 1/(RefinementY+1):1/(RefinementY+1):RefinementY/(RefinementY+1);
iwknots = 1/(RefinementZ+1):1/(RefinementZ+1):RefinementZ/(RefinementZ+1);
pinchedCylinder = nrbkntins(pinchedCylinder, {iuknots,ivknots,iwknots});

% plot nurbs solid

plot_nurbs(pinchedCylinder, 0,1);

view(3);
axis equal
end

