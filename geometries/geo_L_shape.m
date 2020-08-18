function Lshape = geo_L_shape(L,P)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Build a L-shaped plane
%   L - side length
%   P - bottom-left corner POINT coordinates with form of [x, y]
%
%       |----------|
%       |          |
%       |          |
%       |          |
%       |          --------------|
%       |                        |
%       |                        | L
% P(x,y)|________________________|
%
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

if nargin ~= 2    % default parameters
    L = 1.0;  
    P = [-1,-1]; 
end
x = P(1);
y = P(2);
% define control points and knots vector
coefs = zeros(4,4,3);
coefs(:,:,1) = [x, y+2*L, 0, 1; x, y,0, 1; x, y, 0, 1; x+2*L, y, 0, 1;]';
coefs(:,:,2) = [x+L/2, y+2*L, 0, 1; x+L/2, y+L, 0, 1; x+L, y+L/2, 0, 1; x+2*L, y+L/2, 0, 1;]';
coefs(:,:,3) = [x+L, y+2*L, 0, 1; x+L, y+L, 0, 1; x+L, y+L, 0, 1; x+2*L, y+L, 0, 1;]';
knots{1} = [0 0 0 0.5 1 1 1];
knots{2} = [0 0 0 1 1 1];

% build nurbs solid by using control points and knots vector
Lshape = nrbmak(coefs,knots);

% degeree elevate
Lshape = nrbdegelev(Lshape,[1 2]);

% insert knots
RefinementY = 12;    % the number of knots inseted in u direction 
RefinementX = RefinementY*2;    % the number of knots inseted in v direction
iuknots = 1/(RefinementX+1):1/(RefinementX+1):RefinementX/(RefinementX+1);
ivknots = 1/(RefinementY+1):1/(RefinementY+1):RefinementY/(RefinementY+1);
Lshape = nrbkntins(Lshape,{iuknots ivknots});

% plot nurbs patch
plot_nurbs(Lshape, 0, 1);
view(2);
axis equal

end

