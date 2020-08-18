function srf = geo_test_surface
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% construct a test surface
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

coefs = zeros(4,4,4);
r = 0.75;
coefs(:,:,1) = [0,0,0,1; 1,0,r,1; 3,0,r,1; 4,0,0,1]';
coefs(:,:,2) = [0,1,r,1; 1,1,1,1; 3,1,1,1; 4,1,r,1]';
coefs(:,:,3) = [0,3,r,1; 1,3,1,1; 3,3,1,1; 4,3,r,1]';
coefs(:,:,4) = [0,4,0,1; 1,4,r,1; 3,4,r,1; 4,4,0,1]';

knots{1} = [0 0 0 0.5 1 1 1];
knots{2} = [0 0 0 0.5 1 1 1];

srf = nrbmak(coefs,knots);
plot_nurbs(srf, 1, 1);
axis equal
end

