function square = geo_square( pts, length )
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Build a plane square
% Input:
%   pts - coordinates of the bottom-left corner vertex, [x y]
%   length - length of the square
%
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

if nargin ~= 2    % default parameters
    pts = [0,0];  
    length = 1.0; 
end
x = pts(1);
y = pts(2);
L = length;
% define control points and knots vector
coefs = zeros(4,2,2);
coefs(:,:,1) = [x, y, 0, 1; x+L, y, 0, 1;]';
coefs(:,:,2) = [x, y+L, 0, 1; x+L, y+L, 0, 1;]';
knots{1} = [0 0 1 1];
knots{2} = [0 0 1 1];

% build nurbs solid by using control points and knots vector
square = nrbmak(coefs, knots);

% degeree elevate
square = nrbdegelev(square,[2,2]);

% insert knots
RefinementX = 4;    % the number of knots inseted in u direction 
RefinementY = 4;    % the number of knots inseted in v direction
iuknots = 1/(RefinementX+1):1/(RefinementX+1):RefinementX/(RefinementX+1);
ivknots = 1/(RefinementY+1):1/(RefinementY+1):RefinementY/(RefinementY+1);
square = nrbkntins(square, {iuknots ivknots});

% plot nurbs solid

plot_nurbs(square, 0,1);
axis equal;
view(2);

end

