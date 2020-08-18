function plate_hole = geo_plate_with_hole_2patch(length, radius)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Using two patches to build one-quarter of a square plate with hole
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
L = length;
R = radius;
alpha = 22.5/180*pi;
w = cos(alpha);

% define control points and knots vector of patch 1
coefs1 = zeros(4,3,3);
coefs1(:,:,1) = [0, R, 0, 1; R*tan(alpha)*w, R*w, 0, w; R*cos(2*alpha),R*sin(2*alpha), 0, 1]';
coefs1(:,:,2) = [0, (R+L)/2, 0, 1; (R+L)/2*tan(alpha)*w, (R+L)/2*w, 0, w; (R*cos(2*alpha)+L)/2, (R*sin(2*alpha)+L)/2, 0, 1]';
coefs1(:,:,3) = [0, L, 0, 1; L/2*w, L*w, 0, w; L, L, 0, 1]';
knots1{1} = [0 0 0 1 1 1];
knots1{2} = [0 0 0 1 1 1];
% build patch 1 by using control points and knots vector
plate_hole1 = nrbmak(coefs1, knots1);

% define control points and knots vector of patch 2
coefs2 = zeros(4,3,3);
coefs2(:,:,1) = [R*cos(2*alpha),R*sin(2*alpha), 0, 1; R*w, R*tan(alpha)*w,  0, w; R, 0,0, 1]';
coefs2(:,:,2) = [ (R*cos(2*alpha)+L)/2, (R*sin(2*alpha)+L)/2, 0, 1; (R+L)/2*w, (R+L)/2*tan(alpha)*w, 0, w; (R+L)/2,0, 0, 1]';
coefs2(:,:,3) = [ L, L, 0, 1; L*w, L/2*w,  0, w; L, 0, 0, 1;]';
knots2{1} = [0 0 0 1 1 1];
knots2{2} = [0 0 0 1 1 1];
% build patch 2 by using control points and knots vector
plate_hole2 = nrbmak(coefs2, knots2);


% degeree elevate  of patch 1 and 2
plate_hole1 = nrbdegelev(plate_hole1,[1,1]);
plate_hole2 = nrbdegelev(plate_hole2,[1,1]);

% insert knots into patch 1 and patch 2
RefinementX = 2;    % the number of knots inseted in u direction 
RefinementY = 3;    % the number of knots inseted in v direction
iuknots = 1/(RefinementX+1):1/(RefinementX+1):RefinementX/(RefinementX+1);
ivknots = 1/(RefinementY+1):1/(RefinementY+1):RefinementY/(RefinementY+1);
plate_hole1 = nrbkntins(plate_hole1, {iuknots ivknots});
plate_hole2 = nrbkntins(plate_hole2, {iuknots ivknots});


% plot patch 1 and 2
figure
plot_nurbs(plate_hole1, 0,1);
hold on;
plot_nurbs(plate_hole2, 0,1);
axis([-inf,inf,-inf,inf]);
view(2);

plate_hole = cell(1,2);
plate_hole{1} = plate_hole1;
plate_hole{2} = plate_hole2;


end

