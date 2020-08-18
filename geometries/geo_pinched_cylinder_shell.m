function pinchedCylinder = geo_pinched_cylinder_shell(R, L)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Build a pinched cylinder SHELL model
% Note that we only build one-eight of the classical pinched cylinder model
% Input:
%   R - radius 
%   L - length
% 
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

if nargin ~= 2    % default parameters
    R = 300;
    L = 600;
end

theta = 45;
w = cos(theta/180*pi);
% define control points and knots vector
coefs = zeros(4,3,2);
coefs(:,:,1) = [0, 0, R, 1; R*w, 0, R*w, w; R, 0, 0, 1; ]';
coefs(:,:,2) = [0, L/2, R, 1; R*w, L/2*w, R*w, w; R, L/2, 0, 1;]';
knots{1} = [0 0 0 1 1 1];
knots{2} = [0 0 1 1];

% build nurbs solid by using control points and knots vector
pinchedCylinder = nrbmak(coefs,knots);

% degeree elevate
pinchedCylinder = nrbdegelev(pinchedCylinder,[1,2]);

% insert knots
RefinementX = 15;    % the number of knots inseted in u direction 
RefinementY = RefinementX;    % the number of knots inseted in v direction
iuknots = 1/(RefinementX+1):1/(RefinementX+1):RefinementX/(RefinementX+1);
ivknots = 1/(RefinementY+1):1/(RefinementY+1):RefinementY/(RefinementY+1);
pinchedCylinder = nrbkntins(pinchedCylinder, {iuknots,ivknots});

% plot nurbs solid
plot_nurbs(pinchedCylinder, 0,1);
view(3);
axis equal
end

