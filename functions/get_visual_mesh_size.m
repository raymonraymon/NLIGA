function msize = get_visual_mesh_size( vmesh )
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Get the bounding box of the visualized mesh
% Input-
%    vmesh - visualized mesh file
% Output
%    msize - [xmin, xmax, ymin, ymax] 
%          or [xmin, xmax, ymin, ymax, zmin, zmax]
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

if size(vmesh.vertices{1,1},2) == 2        % two dimensions
    xmin = 1e10;  ymin = 1e10; 
    xmax = 1e-10; ymax = 1e-10; 
    for i = 1:length(vmesh.vertices)
        xmaxtmp = max(vmesh.vertices{1,i}(:,1));
        xmintmp = min(vmesh.vertices{1,i}(:,1));
        ymaxtmp = max(vmesh.vertices{1,i}(:,2));
        ymintmp = min(vmesh.vertices{1,i}(:,2));
        if xmin > xmintmp   
            xmin = xmintmp;   
        end
        if xmax < xmaxtmp 
            xmax = xmaxtmp;
        end
        if ymin > ymintmp
            ymin = ymintmp;
        end
        if ymax < ymaxtmp
            ymax = ymaxtmp;
        end
    end
    msize = [xmin, xmax, ymin, ymax];
elseif size(vmesh.vertices{1,1},2) == 3    % three dimensions
    xmin = 1e10;  ymin = 1e10;   zmin = 1e10;
    xmax = 1e-10; ymax = 1e-10;  zmax= 1e-10;
    for i = 1:length(vmesh.vertices)
        xmaxtmp = max(vmesh.vertices{1,i}(:,1));
        xmintmp = min(vmesh.vertices{1,i}(:,1));
        ymaxtmp = max(vmesh.vertices{1,i}(:,2));
        ymintmp = min(vmesh.vertices{1,i}(:,2));
        zmaxtmp = max(vmesh.vertices{1,i}(:,3));
        zmintmp = min(vmesh.vertices{1,i}(:,3));
        if xmin > xmintmp   
            xmin = xmintmp;   
        end
        if xmax < xmaxtmp 
            xmax = xmaxtmp;
        end
        if ymin > ymintmp
            ymin = ymintmp;
        end
        if ymax < ymaxtmp
            ymax = ymaxtmp;
        end
        if zmin > zmintmp
            zmin = zmintmp;
        end
        if zmax < zmaxtmp
            zmax = zmaxtmp;
        end
    end
    msize = [xmin, xmax, ymin, ymax, zmin, zmax];
end



end

