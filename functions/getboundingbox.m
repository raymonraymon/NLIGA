function bb = getboundingbox(points)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% get the bounding box of some points
% 
% acceptable forms:
%       bb = getboundingbox(points)
% 
% inputs:
%       points: some points in 2-D, 3-D or 4-D
%       where size(points,1) = 3
% 
% outputs:
%       bb: coordinates of the bounding box 
%       where bb(3,2) = [xmin xmax; ymin ymax; zmin zmax];
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
tic;

bb = zeros(3,2);
d = ndims(points);          % dimension of points matrix

% get x, y, z coordinates respectively
if d == 4
    x = points(1,:,:,:);
    y = points(2,:,:,:);
    z = points(3,:,:,:);
elseif d == 3
    x = points(1,:,:);
    y = points(2,:,:);
    z = points(3,:,:);
else
    x = points(1,:);
    y = points(2,:);
    z = points(3,:);
end

% get extremum in x,y,z
bb(1,1) = min(x(:));
bb(1,2) = max(x(:));
bb(2,1) = min(y(:));
bb(2,2) = max(y(:));
bb(3,1) = min(z(:));
bb(3,2) = max(z(:));

% disp(['get bounding box time ', num2str(toc)]);

end