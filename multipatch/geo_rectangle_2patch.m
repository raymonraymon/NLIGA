function rectangle = geo_rectangle_2patch
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Using two patches to build a rectangle
%
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

rectangle = cell(1,2);
L = 1.0;
% define control points and knots vector of patch 1
coefs1 = zeros(4,2,3);
coefs1(:,:,1) = [0, 0, 0, 1; 1,0,0,1]';
coefs1(:,:,2) = [0,0.5*L,0,1; 0.75*L,0.5*L,0,1;]';
coefs1(:,:,3) = [0,1*L,0,1; 1*L,1*L,0,1]';
knots1{1} = [0 0 1 1];
knots1{2} = [0 0 0 1 1 1];
% build patch 1 by using control points and knots vector
rectangle{1,1} = nrbmak(coefs1, knots1);

% define control points and knots vector of patch 2
coefs2 = zeros(4,2,3);
coefs2(:,:,1) = [L,0,0,1; 2*L,0,0,1; ]';
coefs2(:,:,2) = [0.75*L, 0.5*L, 0, 1; 2*L,0.5*L,0,1;]';
coefs2(:,:,3) = [L,L,0,1; 2*L,L,0,1]';
coefs2(1,:,:) = coefs2(1,:,:) + 0.08;
knots2{1} = [0 0 1 1];
knots2{2} = [0 0 0 1 1 1];
% build patch 2 by using control points and knots vector
rectangle{1,2} = nrbmak(coefs2, knots2);

% degeree elevate  of patch 1 and 2
rectangle{1,1} = nrbdegelev(rectangle{1,1},[2,1]);
rectangle{1,2} = nrbdegelev(rectangle{1,2},[2,1]);

% insert knots into patch 1 and patch 2
RefinementX = 3;    % the number of knots inseted in u direction 
RefinementY = 3;    % the number of knots inseted in v direction
iuknots = 1/(RefinementX+1):1/(RefinementX+1):RefinementX/(RefinementX+1);
ivknots = 1/(RefinementY+1):1/(RefinementY+1):RefinementY/(RefinementY+1);
rectangle{1,1} = nrbkntins(rectangle{1,1}, {iuknots ivknots});
rectangle{1,2} = nrbkntins(rectangle{1,2}, {iuknots ivknots});

figure
for i = 1:length(rectangle)
    hold on;
    plot_nurbs(rectangle{1,i}, 0,1);
end
axis([-inf,inf,-inf,inf,-inf,inf]);
view(2);


end

