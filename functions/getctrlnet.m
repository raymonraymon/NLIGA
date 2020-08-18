function [center, edge] = getctrlnet(points)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% get control net of nurbs curve, surface or volume
% 
% acceptable forms:
%       [center, edge] = getctrlnet(nurbs)
% 
% inputs:
%       points: control points [x,y,z] of nurbs curve, surface or volume
%       where points in 2-D, 3-D or 4-D, size(points,1) = 3
% 
% outputs:
%       center: coordinates of control points 
%       edge:   edge correlation of control net
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
tic;

d = ndims(points);                          % dimension of points matrix

% get control points and edges of control net  
if d == 4                                   % get control net of nurbs volume  
%     get all control points
    ctrlnum = size(points,2)*size(points,3)*size(points,4);
    center = zeros(ctrlnum,3);
    c = 1;
    for k = 1:size(points,4)
        for j = 1:size(points,3)
            for i = 1:size(points,2)
                center(c,:) = reshape (points(:,i,j,k), 1, 3);
                c = c+1;
            end
        end
    end

%         get all edges of control net
    edgenum = (size(points,2)*(size(points,3)-1)+(size(points,2)-1)*size(points,3))*size(points,4)...
                +size(points,2)*size(points,3)*(size(points,4)-1);
    edge = zeros(edgenum,6);
    c = 1;
    for ii = 1:size(points,2)
        for jj = 1:size(points,3)
            coefs = reshape (points(1:3,ii,jj,:), 3, []);
            for i = 2:size(coefs,2)
                edge(c,:) = reshape([coefs(:,i-1);coefs(:,i)],1,6);
                c = c+1;
            end                
        end
        for kk = 1:size(points, 4)
            coefs = reshape (points(1:3,ii,:,kk), 3, []);
            for i = 2:size(coefs,2)
                edge(c,:) = reshape([coefs(:,i-1);coefs(:,i)],1,6);
                c = c+1;
            end
        end
    end
    for jj = 1:size(points, 3)
        for kk = 1:size(points, 4)
            coefs = reshape (points(1:3,:,jj,kk), 3, []);
            for i = 2:size(coefs,2)
                edge(c,:) = reshape([coefs(:,i-1);coefs(:,i)],1,6);
                c = c+1;
            end
        end
    end

elseif d == 3                               % plot control net of nurbs surface
%         get all control points
    ctrlnum = size(points,2)*size(points,3);
    center = zeros(ctrlnum,3);
    c = 1;
    for j = 1:size(points,3)
        for i = 1:size(points,2)
            center(c,:) = reshape (points(:,i,j), 1, 3);
            c = c+1;
        end
    end

%         get all edges of control net
    edgenum = size(points,2)*(size(points,3)-1)+(size(points,2)-1)*size(points,3);
    edge = zeros(edgenum,6);
    c = 1;        
    for ii = 1:size(points, 2)
        coefs = reshape(points(1:3,ii,:), 3, []);
        for i = 2:size(coefs,2)
            edge(c,:) = reshape([coefs(:,i-1);coefs(:,i)],1,6);
            c = c+1;
        end
    end
    for jj = 1:size(points, 3)
        coefs = reshape(points(1:3,:,jj), 3, []);
        for i = 2:size(coefs,2)
            edge(c,:) = reshape([coefs(:,i-1);coefs(:,i)],1,6);
            c = c+1;
        end                    
    end

else                                        % plot control net of nurbs curve
% 	get all control points  
    ctrlnum = size(points,2);
    center = zeros(ctrlnum,3);
    for i = 1:ctrlnum
        center(i,:) = reshape(points(:,i),1,3);
    end
%     get all edges of control net   
    edgenum = size(points,2)-1;
    edge = zeros(edgenum,6);
    for i = 2:size(points,2)
        edge(i-1,:) = reshape([points(:,i-1);points(:,i)],1,6);
    end
end

% disp(['get control net time ', num2str(toc)]);

end