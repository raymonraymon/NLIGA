function tri_mesh = build_visual_mesh_suf( num1, num2 )
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%  Build parametric triangle structure for visualization
%  Input:
%    num1,  segment number for u-side direction
%    num2,  segment number for v-side direction
%  Output:
%    tri_mesh - visualized triangular mesh
%               tri_mesh.trinum - total triangles
%               tri_mesh.tripts - coordinates of triangular vertices
%               tri_mesh.trimesh - element-node connectivity
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%


trinum = num1*num2*2;        % total elements
trimesh = zeros(trinum,3);   % mesh connectivity
count = 1;
for j = 1:num2
    for i = 1:num1
        p1 = (j-1)*(num1+1)+i;
        p2 = (j-1)*(num1+1)+i+1;
        p3 = (j-1)*(num1+1)+i+1+(num1+1);
        p4 = (j-1)*(num1+1)+i+(num1+1);
        trimesh(count,:) = [p1 p2 p4];
        trimesh(count+1,:) = [p2 p3 p4];
        count = count+2;
    end
end

offset = 0;
uu = linspace(0+offset,1-offset,num1+1);
vv = linspace(0+offset,1-offset,num2+1);
count = 1;
tripts = zeros((num1+1)*(num2+1),2);
for j = 1:num2+1
    for i = 1:num1+1
        tripts(count,:) = [uu(i), vv(j)];   % parametric nodal points
        count = count+1;
    end
end

tri_mesh.trinum = trinum;
tri_mesh.tripts = tripts;
tri_mesh.trimesh = trimesh;

end

