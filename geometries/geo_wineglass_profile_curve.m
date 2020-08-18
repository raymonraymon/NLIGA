function crv = geo_wineglass_profile_curve
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% build the profile curve of a red wine glass
%
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

coefs = [0,-36,0,1;   0,-35.5,2,1; 0,-4,2,1; 
         0,-4,30,1;   0,-4,50,1;   0,-4,80,1; 
         0,-44,120,1; 0,-36,160,1; 0,-28,200,1]';
knots = [0 0 0 0.02 0.15 0.3 0.5 0.6 0.9 1 1 1];
crv = nrbmak(coefs, knots);

% plot_nurbs(crv, 1, 1);
nrbctrlplot(crv);
view([1,0,0]);
axis equal
end

