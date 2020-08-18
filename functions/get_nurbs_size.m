function msize = get_nurbs_size( nrb )
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Get the bounding box of a NURBS structure
% Input-
%    nrb - NURBS structure, curve, surface, or volume
% Output
%    msize - [xmin, xmax, ymin, ymax, zmin, zmax]
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

xmin = 1e10;  ymin = 1e10;   zmin = 1e10;
xmax = 1e-10; ymax = 1e-10;  zmax= 1e-10;

len = length(nrb.number);
coefs = nrb.coefs;
if len == 1       % one-dimensional 
    for i = 1:3    % convert the coordinates to Euclidean space
        coefs(i,:) = coefs(i,:)./coefs(4,:);
    end
    xmax = max( coefs(1,:) );
    xmin = min( coefs(1,:) );
    ymax = max( coefs(2,:) );
    ymin = min( coefs(2,:) ); 
    zmax = max( coefs(3,:) );
    zmin = min( coefs(3,:) );
elseif len == 2   % two-dimensional
    for i = 1:3    % convert the coordinates to Euclidean space
        for j = 1:size(coefs,3)
            coefs(i,:,j) = coefs(i,:,j)./coefs(4,:,j);
        end
    end
    xmax = max(max( coefs(1,:,:) ));
    xmin = min(min( coefs(1,:,:) ));
    ymax = max(max( coefs(2,:,:) ));
    ymin = min(min( coefs(2,:,:) )); 
    zmax = max(max( coefs(3,:,:) ));
    zmin = min(min( coefs(3,:,:) ));
elseif len == 3   % three-dimensional
    for i = 1:3    % convert the coordinates to Euclidean space
        for j = 1:size(coefs,3)
            for k = 1:size(coefs,4)
                coefs(i,:,j,k) = coefs(i,:,j,k)./coefs(4,:,j,k);
            end
        end
    end
    xmax = max(max(max( coefs(1,:,:,:) )));
    xmin = min(min(min( coefs(1,:,:,:) )));
    ymax = max(max(max( coefs(2,:,:,:) )));
    ymin = min(min(min( coefs(2,:,:,:) ))); 
    zmax = max(max(max( coefs(3,:,:,:) )));
    zmin = min(min(min( coefs(3,:,:,:) )));    
end


msize = [xmin, xmax, ymin, ymax, zmin, zmax];


end

