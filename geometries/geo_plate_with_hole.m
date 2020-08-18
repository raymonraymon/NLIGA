function plate_hole = geo_plate_with_hole(length, radius)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% build one-quarter of a square plate with hole
% Input:
%   length - length of the side = 1/2 side length of the initial square
%   radius - radius of the hole
%
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

if nargin ~= 2    % default parameters
    length = 4;  % side length
    radius = 1;  % radius
end
% define control points and knots vector
L = length;
R = radius;

coefs = zeros(4,4,3);
alpha = 22.5/180*pi;
w = cos(alpha)^2;
coefs(:,:,1) = [0, R, 0, 1; R*tan(alpha)*w, R*w, 0, w; R*w, R*tan(alpha)*w, 0, w; R, 0, 0, 1]';
coefs(:,:,2) = [0, (R+L)/2, 0, 1; (R+L)/2*tan(alpha)*w, (R+L)/2*w, 0, w; (R+L)/2*w, (R+L)/2*tan(alpha)*w, 0, w; (R+L)/2, 0, 0, 1]';
coefs(:,:,3) = [0, L, 0, 1; L*w, L*w, 0, w; L*w, L*w, 0, w; L, 0, 0, 1]';
knots{1} = [0 0 0 0.5 1 1 1];
knots{2} = [0 0 0 1 1 1];

% build nurbs solid by using control points and knots vector
plate_hole = nrbmak(coefs, knots);

% degeree elevate
plate_hole = nrbdegelev(plate_hole,[1,1]);

% insert knots
RefinementX = 1;    % the number of knots inseted in u direction 
RefinementY = RefinementX;    % the number of knots inseted in v direction
iuknots = 1/(RefinementX+1):1/(RefinementX+1):RefinementX/(RefinementX+1);
ivknots = 1/(RefinementY+1):1/(RefinementY+1):RefinementY/(RefinementY+1);
plate_hole = nrbkntins(plate_hole, {iuknots, ivknots});

% plot nurbs solid
plot_nurbs(plate_hole, 0,1);
view(2);
axis equal
end

