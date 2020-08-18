function cooks = geo_cooks_membrane
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% build cooks model
%                          /|
%                      /    |
%                  /        |  L3
%               /           |
%            |            /
%            |          /
%     L2     |        /
%            |      /
%            |    /
%            |  /
%                    L1
%
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

% define control points and knots vector
L1 = 48;
L2 = 44;
L3 = 16;
coefs = zeros(4,2,2);
coefs(:,:,1) = [0 0 0 1; L1, L2, 0, 1]';
coefs(:,:,2) = [0, L2, 0, 1; L1, L2+L3, 0, 1]';
knots{1} = [0 0  1 1];
knots{2} = [0 0  1 1];

% build nurbs solid by using control points and knots vector
cooks = nrbmak(coefs, knots);

% degeree elevate
cooks = nrbdegelev(cooks,[1,1]);

% insert knots
RefinementX = 2;    % the number of knots inseted in u direction 
RefinementY = 2;    % the number of knots inseted in v direction
iuknots = 1/(RefinementX+1):1/(RefinementX+1):RefinementX/(RefinementX+1);
ivknots = 1/(RefinementY+1):1/(RefinementY+1):RefinementY/(RefinementY+1);
cooks = nrbkntins(cooks, {iuknots ivknots});

% plot nurbs solid
plot_nurbs(cooks, 1,1);
view(2);
axis equal
end

