function srf = geo_wineglass_surface
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% build a wine glass surface from the prfile curve
% 'geo_wineglass_profile_curve' by using revolving
%
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

crv = geo_wineglass_profile_curve;
pnt = [0 0 0];
vec = [0 0 1];
srf = nrbrevolve(crv,pnt,vec);
msize = get_nurbs_size( srf );
axis(msize);
clf;

plot_nurbs(srf, 1, 1);
axis equal
end

