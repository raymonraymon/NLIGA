function solid = geo_quarter_hollow_cylinder3d
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% build a quarter of hollow cylinder with inder radius r and outer radius R
% the bottom face is on yz plane
%
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%


% define control points and knots vector
r = 7;   % inner radius
R = 10;  % outer radius
L = 15;  % length
w = sqrt(2)/2;   % weight
coefs = zeros(4,3,2,2);
coefs(:,:,1,1) = [r,0,0,1; r*w,0,r*w,w; 0,0,r,1;]';
coefs(:,:,2,1) = [R,0,0,1; R*w,0,R*w,w; 0,0,R,1; ]';
coefs(:,:,1,2) = [r,L,0,1; r*w,L*w,r*w,w; 0,L,r,1;]';
coefs(:,:,2,2) = [R,L,0,1; R*w,L*w,R*w,w; 0,L,R,1;]';
knots{1} = [0 0 0 1 1 1];
knots{2} = [0 0 1 1];
knots{3} = [0 0 1 1];

% build nurbs solid by using control points and knots vector
solid = nrbmak(coefs, knots);

% degeree elevate
solid = nrbdegelev(solid,[1,2,2]);

% insert knots
iuknots = [0.5 ];
ivknots = [ ];
iwknots = [0.5 ];
solid = nrbkntins(solid, {iuknots ivknots iwknots});

% plot nurbs solid
plot_nurbs(solid, 1,1);
axis equal
end

