function Lshape = geo_Lshape_3patch(L,P)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Build a conforming two-patch L-shaped plane model
%   L - side length
%   P - bottom-left corner POINT coordinates with form of [x, y]
%   see the same model given in the function 'geo_L_shape' and 'geo_Lshape_2patch'
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

Lshape = cell(1,3);
% define control points and knots vector of patch 1
coefs1 = zeros(4,2,2);
coefs1(:,:,1) = [x, y, 0, 1; x+L, y, 0, 1]';
coefs1(:,:,2) = [x, y+L, 0, 1; x+L, y+L, 0, 1]';
knots1{1} = [0 0 1 1];
knots1{2} = [0 0 1 1];
% build patch 1 by using control points and knots vector
Lshape{1,1} = nrbmak(coefs1, knots1);

% define control points and knots vector of patch 2
coefs2 = zeros(4,2,2);
coefs2(:,:,1) = [x+L, y, 0, 1; x+2*L, y, 0, 1]';
coefs2(:,:,2) = [x+L, y+L, 0, 1; x+2*L, y+L, 0, 1]';
knots2{1} = [0 0 1 1];
knots2{2} = [0 0 1 1];
% build patch 2 by using control points and knots vector
Lshape{1,2} = nrbmak(coefs2, knots2);

% define control points and knots vector of patch 3
coefs3 = zeros(4,2,2);
coefs3(:,:,1) = [x, y+L, 0, 1; x+L, y+L, 0, 1]';
coefs3(:,:,2) = [x, y+2*L, 0, 1; x+L, y+2*L, 0, 1]';
knots3{1} = [0 0 1 1];
knots3{2} = [0 0 1 1];
% build patch 2 by using control points and knots vector
Lshape{1,3} = nrbmak(coefs3, knots3);


% degeree elevate  of patch 1 and 2
Lshape{1,1} = nrbdegelev(Lshape{1,1},[2,2]);
Lshape{1,2} = nrbdegelev(Lshape{1,2},[2,2]);
Lshape{1,3} = nrbdegelev(Lshape{1,3},[2,2]);

% insert knots into patch 1 and patch 2
RefinementX = 9;    % the number of knots inseted in u direction 
RefinementY = RefinementX;    % the number of knots inseted in v direction
iuknots = 1/(RefinementX+1):1/(RefinementX+1):RefinementX/(RefinementX+1);
ivknots = 1/(RefinementY+1):1/(RefinementY+1):RefinementY/(RefinementY+1);
Lshape{1,1} = nrbkntins(Lshape{1,1}, {iuknots ivknots});
Lshape{1,2} = nrbkntins(Lshape{1,2}, {iuknots ivknots});
Lshape{1,3} = nrbkntins(Lshape{1,3}, {iuknots ivknots});

% plot patch 1 and 2
figure
plot_nurbs(Lshape{1,1}, 0,1);
hold on;
plot_nurbs(Lshape{1,2}, 0,1);
plot_nurbs(Lshape{1,3}, 0,1);
axis([-inf,inf,-inf,inf]);
view(2);




end

