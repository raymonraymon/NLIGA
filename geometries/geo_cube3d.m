function solid = geo_cube3d( pts, length )
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% build a unit solid model
% Input:
%   pts - coordinates the the bottom-left corner vertex on the front face, [x y z]
%   length - length of the square
%
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

if nargin ~= 2    % default parameters
    pts = [0,0,0];
    length = 1.0;
end

% define control points and knots vector
x = pts(1);
y = pts(2);
z = pts(3);
L = length;
coefs = zeros(4,2,2,2);
coefs(:,:,1,1) = [x, y,   z,   1; x+L, y,   z,   1]';
coefs(:,:,2,1) = [x, y+L, z,   1; x+L, y+L, z,   1]';
coefs(:,:,1,2) = [x, y,   z+L, 1; x+L, y,   z+L, 1]';
coefs(:,:,2,2) = [x, y+L, z+L, 1; x+L, y+L, z+L, 1]';

knots{1} = [0 0 1 1];
knots{2} = [0 0 1 1];
knots{3} = [0 0 1 1];

% build nurbs solid by using control points and knots vector
solid = nrbmak(coefs, knots);

% degeree elevate
solid = nrbdegelev(solid,[2,2,2]);

% insert knots
iuknots = [0.5];
ivknots = [0.5];
iwknots = [0.5];
solid = nrbkntins(solid, {iuknots ivknots iwknots});

% plot nurbs solid
plot_nurbs(solid, 1,1);

axis equal

end

