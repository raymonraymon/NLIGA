function crv = geo_closed_curve
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% build a cubic closed curve
%
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

% define the control points
P0 = [1.5, 0, 0, 1];
P1 = [0.5, 0.5, 0, 1];
P2 = [-0.5, 0.5, 0, 1];
P3 = [-1.5, 0, 0, 1];
P4 = [-0.5, -0.5, 0, 1];
P5 = [0.5, -0.5, 0, 1];
P6 = P0;
P7 = P1;
P8 = P2;
coefs = [P0; P1; P2; P3; P4; P5; P6; P7; P8]';

% define knots vector
knots = [-1, -2/6, -1/6, 0, 1/6, 2/6, 3/6, 4/6, 5/6, 1, 1+1/6, 1+2/6, 1+4/6];

% construct nurbs curve`1
crv = nrbmak(coefs,knots);

plot_nurbs(crv, 1,1);
view(2);
axis equal
end

