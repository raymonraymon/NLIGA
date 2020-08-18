function circle = geo_circle( center, radius )
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Build a circular plate
% Input:
%   center - coordinates of the center of the circular plate [x y]
%   radius - radius of the circle
%
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

if nargin ~= 2    % default parameters
    center = [0,0];
    radius = 1.0;
end

x = center(1);
y = center(2);
r = radius;
rad = pi/180;
w = cos(45*rad);

% define control points and knots vector
coefs = zeros(4,3,3);
coefs(:,:,1) = [x-r*cos(rad*45)    y-r*sin(45*rad)    0    1;
                            x*w    (y-sqrt(2)*r)*w    0    w;
                x+r*cos(rad*45)    y-r*sin(45*rad)    0    1]';
coefs(:,:,2) = [(x-sqrt(2)*r)*w                y*w    0    w;
                              x                  y    0    1;
                (x+sqrt(2)*r)*w                y*w    0    w]';
coefs(:,:,3) = [x-r*cos(rad*45)    y+r*sin(45*rad)    0    1;
                            x*w    (y+sqrt(2)*r)*w    0    w;
                x+r*cos(rad*45)    y+r*sin(45*rad)    0    1]';
knots{1}=[0 0 0 1 1 1];
knots{2}=[0 0 0 1 1 1];

% build nurbs solid by using control points and knots vector
circle = nrbmak(coefs, knots);

% degeree elevate
circle = nrbdegelev(circle,[1,1]);

% insert knots
RefinementX = 7;    % the number of knots inseted in u direction 
RefinementY = 7;    % the number of knots inseted in v direction
iuknots = 1/(RefinementX+1):1/(RefinementX+1):RefinementX/(RefinementX+1);
ivknots = 1/(RefinementY+1):1/(RefinementY+1):RefinementY/(RefinementY+1);
circle = nrbkntins(circle, {iuknots ivknots});

% plot nurbs solid

plot_nurbs(circle, 0,1);

view(2);
axis equal
end

