function curved_beam = geo_curved_beam
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% build curved beam
%
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

% define control points and knots vector
R = 20;
T = 1;
w = sqrt(2)/2;
coefs = zeros(4,3,2);
coefs(:,:,1) = [0 R 0 1; R*w, R*w, 0, w; R, 0, 0, 1]';
coefs(:,:,2) = [0 (R+T) 0 1; (R+T)*w, (R+T)*w, 0, w; (R+T),0, 0, 1]';
knots{1} = [0 0 0 1 1 1];
knots{2} = [0 0  1 1];

% build nurbs solid by using control points and knots vector
curved_beam = nrbmak(coefs, knots);

% degeree elevate
curved_beam = nrbdegelev(curved_beam,[1,2]);

% insert knots
RefinementX = 64;    % the number of knots inseted in u direction 
RefinementY = 2;    % the number of knots inseted in v direction
iuknots = 1/(RefinementX+1):1/(RefinementX+1):RefinementX/(RefinementX+1);
ivknots = 1/(RefinementY+1):1/(RefinementY+1):RefinementY/(RefinementY+1);
curved_beam = nrbkntins(curved_beam, {iuknots ivknots});

% plot nurbs solid
plot_nurbs(curved_beam, 0,1);
axis equal;
view(2);

end


