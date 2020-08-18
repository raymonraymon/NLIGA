function solid = geo_wineglass_volume
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% build a wine galss volume from the surface
% 'geo_wineglass_surface'
%
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

srf = geo_wineglass_surface;
coef1 = srf.coefs;
coefs = zeros([size(coef1),2]);
coefs(:,:,:,1) = coef1;
coef2 = coef1;
scale = 0.9;
for i = 1:3   % convert to Euclidean space
    for j = 1:size(coef2,3)
        coef2(i,:,j) = coef2(i,:,j)./coef2(4,:,j); 
    end
end
for i = 1:2   % scale x and y coordinates
    for j = 1:size(coef2,3)
        coef2(i,:,j) = coef2(i,:,j).*scale;
    end
end

for i = 1:3   % convert back to homogeneous space
    for j = 1:size(coef2,3)
        coef2(i,:,j) = coef2(i,:,j).*coef2(4,:,j); 
    end
end

coefs(:,:,:,2) = coef2;
knots3{1} =  srf.knots{1,1};
knots3{2} =  srf.knots{1,2};
knots3{3} = [0 0 1 1];
solid = nrbmak(coefs,knots3);
clf;
plot_nurbs(solid, 1, 1);
axis equal
end

