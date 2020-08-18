function vmesh = plot_visual_mesh3d(flag, filename)
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Plot the results of visualized 3d mesh model
%  Input:
%    flag - color map: 1-U1, 2-U2, 3-U3, 4, U magnitude, 5-S11, 6-S22
%                      7-S33, 8-S12, 9-S23, 10-S31, 11-mises
%    filename - file name for read
%  Output:
%    vmesh - visualized mesh structure
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

fname = fileparts(mfilename('fullpath'));        % get current path
index_dir = strfind(fname,'\');                  
str_temp = fname(1:index_dir(end));              % get the name of previous folder
fname = [str_temp,'output\', filename, '.msh'];  % get the file saved in the folder ...\Output\.msh
vmesh = read_visual_mesh( fname );               % read the saved .msh file

% % make movie
fig = figure;
% videoName = [str_temp, 'output\', filename, '.mp4'];
% aviobj = VideoWriter(videoName);
% aviobj.FrameRate = 2; % frame rate
% open(aviobj);

face = cell2mat(vmesh.face(1));
maxnum = max(max(face));
vertices = cell2mat(vmesh.vertices(1));
trivertex = vertices(1:maxnum,:);
displacement = cell2mat(vmesh.displacement(1));
stress = cell2mat(vmesh.stress(1));  
linmesh = cell2mat(vmesh.linmesh(1));
p = patch('Faces',face, 'Vertices', trivertex);

% principal stress
if flag == 1
    cdata = displacement(1:maxnum,1);
elseif flag == 2
    cdata = displacement(1:maxnum,2);
elseif flag == 3
    cdata = displacement(1:maxnum,3);  
elseif flag == 4
    cdata = zeros(maxnum,1);
    for k = 1:maxnum
        cdata(k) = norm(displacement(k,:));
    end 
elseif flag == 5
    cdata = stress(1:maxnum,1);
elseif flag == 6
    cdata = stress(1:maxnum,2);
elseif flag == 7
    cdata = stress(1:maxnum,3); 
elseif flag == 8
    cdata = stress(1:maxnum,4);
elseif flag == 9
    cdata = stress(1:maxnum,5);
elseif flag == 10
    cdata = stress(1:maxnum,6); 
elseif flag == 11
    cdata = von_mises( stress(1:maxnum,:) );
end
set(p,'FaceColor','interp','FaceVertexCData',cdata);
set(p,'EdgeColor','none');
hold on; % plot knot curves
for i = 1:size(linmesh,1)
    vv = vertices(linmesh(i,:),:);
    aa(i) = plot3(vv(:,1), vv(:,2), vv(:,3), 'k-');
end
hold off;
axis equal;
msize = get_visual_mesh_size(vmesh);
axis(msize);
box on;
colorbar;
xlabel('x');
ylabel('y');
zlabel('z');
% set(gcf,'color','black');
% MOV = getframe(fig);
% writeVideo(aviobj,MOV);
% fmat(1,1) = getframe;
pause(1.0);

if length(vmesh.vertices) > 1
    for i = 2:length(vmesh.vertices)
        face = cell2mat(vmesh.face(i));
        maxnum = max(max(face));
        vertices = cell2mat(vmesh.vertices(i));
        trivertex = vertices(1:maxnum,:);
        displacement = cell2mat(vmesh.displacement(i));
        stress = cell2mat(vmesh.stress(i));    
        linmesh = cell2mat(vmesh.linmesh(i));
        if flag == 1
            cdata = displacement(1:maxnum,1);
        elseif flag == 2
            cdata = displacement(1:maxnum,2);
        elseif flag == 3
            cdata = displacement(1:maxnum,3);
        elseif flag == 4
            cdata = zeros(maxnum,1);
            for k = 1:maxnum
                cdata(k) = norm(displacement(k,:));
            end 
        elseif flag == 5
            cdata = stress(1:maxnum,1);
        elseif flag == 6
            cdata = stress(1:maxnum,2);
        elseif flag == 7
            cdata = stress(1:maxnum,3); 
        elseif flag == 8
            cdata = stress(1:maxnum,4);
        elseif flag == 9
            cdata = stress(1:maxnum,5);
        elseif flag == 10
            cdata = stress(1:maxnum,6); 
        elseif flag == 11
            cdata = von_mises( stress );
        end
        set(p, 'Faces',face, 'Vertices', trivertex, 'FaceVertexCData',cdata);
        for j = 1:size(linmesh,1)
            vv = vertices(linmesh(j,:),:);
            set(aa(j),'xdata', vv(:,1), 'ydata', vv(:,2),'zdata', vv(:,3));
        end
        title(['Step ',num2str(i)]);
        hcb = colorbar;
        if flag == 1,  title(hcb,'U_x');
        elseif flag == 2,  title(hcb,'U_y');
        elseif flag == 3,  title(hcb,'U_z');
        elseif flag == 4,  title(hcb,'U Magnitude');
        elseif flag == 5,  title(hcb,'S_{11}');
        elseif flag == 6,  title(hcb,'S_{22}');
        elseif flag == 7,  title(hcb,'S_{33}');
        elseif flag == 8,  title(hcb,'S_{12}');
        elseif flag == 9,  title(hcb,'S_{23}');
        elseif flag == 10,  title(hcb,'S_{31}');
        elseif flag == 11,  title(hcb,'Mises');
        end
        drawnow;
%         MOV = getframe(fig);
%         writeVideo(aviobj,MOV);
%         fmat(i+1,1) = getframe;   
        pause(0.5);
    end
end
% pause(2.0);
% close(aviobj);
end

