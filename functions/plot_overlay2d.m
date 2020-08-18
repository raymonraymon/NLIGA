function plot_overlay2d( flag, filename, steps )
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Plot the overlay of the two-dimensional results.
%  Input:
%    flag - color map: 1-U1, 2-U2, 3-U magnitude, 4-S11, 5-S22, 6-S12, 7-mises
%    filename - filename for read
%    steps - the steps to be plot
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

outstr = zeros(length(vmesh.stress),1);
ux = zeros(length(vmesh.stress),1);
uy = zeros(length(vmesh.stress),1);
for i = 1:length(vmesh.stress)
    outstr(i) = von_mises(vmesh.stress{1,i}(2,:));
    ux(i) = vmesh.displacement{1,i}(2,1);
    uy(i) = vmesh.displacement{1,i}(2,2);
end


for i = 1:length(steps)
    j = steps(i);
    face = cell2mat(vmesh.face(j));
    maxnum = max(max(face));
    vertices = cell2mat(vmesh.vertices(j));
    trivertex = vertices(1:maxnum,:);
    displacement = cell2mat(vmesh.displacement(j));
    stress = cell2mat(vmesh.stress(j));    
    linmesh = cell2mat(vmesh.linmesh(j));
    hold on;
    p = patch('Faces',face, 'Vertices', trivertex);   
    hcb = colorbar;
    if flag == 1
        cdata = displacement(1:maxnum,1);
        title(hcb,'U_x');
    elseif flag == 2
        cdata = displacement(1:maxnum,2);
        title(hcb,'U_y');
    elseif flag == 3
        cdata = sqrt(displacement(1:maxnum,1).^2 + displacement(1:maxnum,2).^2);
        title(hcb,'U Magnitude');
    elseif flag == 4
        cdata = stress(1:maxnum,1);
        title(hcb,'S_{11}');
    elseif flag == 5
        cdata = stress(1:maxnum,2);
        title(hcb,'S_{22}');
    elseif flag == 6
        cdata = stress(1:maxnum,3);
        title(hcb,'S_{12}');
    elseif flag == 7
        cdata = von_mises( stress(1:maxnum,:) );
        title(hcb,'Mises');
    end
  
    set(p,'FaceColor','interp','FaceVertexCData',cdata);
    set(p,'EdgeColor','none');
    for k = 1:size(linmesh,1)
        vv = vertices(linmesh(k,:),:);
        aa(k) = plot(vv(:,1), vv(:,2), 'k-');
    end  
end
axis equal;

end

