%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%  This script will show the node numbering of the two-patch rectangle
%
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%


clc;
clear all;

% define the geometry
geo = geo_rectangle_2patch;
% alpha(0)
% build iga mesh structure
mesh = cell(1,length(geo));
for i = 1:length(geo)
    mesh{1,i} = build_iga_mesh( geo{1,i} );   
    mesh{1,i}.gloElNodeCnt = mesh{1,i}.elNodeCnt;
    mesh{1,i}.nodeNum = zeros(mesh{1,i}.nCptsV, mesh{1,i}.nCptsU);
end

% node numbering
for j = 1:mesh{1,1}.nCptsV
    for i = 1:mesh{1,1}.nCptsU
        mesh{1,1}.nodeNum(j,i) = (j-1)*(mesh{1,1}.nCptsU)+i;
    end
end
mesh{1,2}.nodeNum(:,1) = mesh{1,1}.nodeNum(:,end);
for j = 1:mesh{1,2}.nCptsV
    for i = 2:mesh{1,2}.nCptsU
        mesh{1,2}.nodeNum(j,i) = (j-1)*(mesh{1,2}.nCptsU-1)+i-1 + mesh{1,1}.nCpts;
    end
end

% reconstruct global nodes coordinates
globNodeNum = mesh{1,1}.nCpts + mesh{1,2}.nCpts - mesh{1,2}.nCptsV;
globCoords = zeros(globNodeNum,4);
globCoords(1:mesh{1,1}.nCpts,:) = mesh{1,1}.coords;
count = mesh{1,1}.nCpts;
for j = 1:mesh{1,2}.nCptsV
    for i = 2:mesh{1,2}.nCptsU
        count = count+1;
        k = (j-1)*(mesh{1,2}.nCptsU) +i;
        globCoords(count,:) = mesh{1,2}.coords(k,:);
    end
end

min_x = min(globCoords(:,1));
max_x = max(globCoords(:,1));
min_y = min(globCoords(:,2));
max_y = max(globCoords(:,2));
len = sqrt((max_x-min_x)^2 + (max_y-min_y)^2);
hold on;
plot3(globCoords(:,1),globCoords(:,2),globCoords(:,3),'ro');
for kk = 1:size(globCoords,1)
    text(globCoords(kk,1)+len/150,globCoords(kk,2)+len/150,globCoords(kk,3)+0.001,num2str(kk),...
        'color','k','FontSize',20,'FontName','Times New Roman');
end

rightInterPts = zeros(mesh{1,2}.nCptsV,4);
for ii = 1:mesh{1,2}.nCptsV
    aa = ii*mesh{1,1}.nCptsU;
    bb = (ii-1)*mesh{1,2}.nCptsU + 1;
    rightInterPts(ii,:) = mesh{1,2}.coords(bb,:);
    text(rightInterPts(ii,1)+len/150,rightInterPts(ii,2)+len/150,rightInterPts(ii,3)+0.001,num2str(aa),...
        'color','k','FontSize',20,'FontName','Times New Roman');
end
plot3(rightInterPts(:,1),rightInterPts(:,2),rightInterPts(:,3),'ro');
axis off
% axis([-0.2,2.2,-0.2,1.2]);
title('Global Nodes Numbering','color','k','FontSize',12);

geo = geo_rectangle_2patch;
hold on;
for kk = 1:length(mesh)
    plot3(mesh{1,kk}.coords(:,1),mesh{1,kk}.coords(:,2),mesh{1,kk}.coords(:,3)+0.001,'ro');
    for jj = 1:size(mesh{1,kk}.coords,1)
        text(mesh{1,kk}.coords(jj,1)+len/150,mesh{1,kk}.coords(jj,2)+len/150,mesh{1,kk}.coords(jj,3)+0.001,num2str(jj),'color','k','FontSize',20,'FontName','Times New Roman')
    end
end
axis off
% axis([-0.2,2.2,-0.2,1.2]);
title('Local Nodes Numbering','color','k','FontSize',12);


