function elem_index = find_point_span( mesh, pts )
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% find the elem number of a point belongs to
% geo: nurbs strucutre
% pts: parametric points with format: [u0, v0; u1, v1; ...]
%  Find the element that a point belongs to 
%  Input:
%    mesh - iga mesh structure
%    pts - parametric points with format: [u0, v0; u1, v1; ...]
%  Output:
%    elem_index - obtained element index
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

if mesh.dim == 2
    uknot = unique( mesh.uKnots );
    vknot = unique( mesh.vKnots );
    elem_u = length(uknot)-1;
    elem_v = length(vknot)-1;
    elem_index = zeros(1,size(pts,1));

    for i=1:size(pts,1)
        for j = 1:elem_u
            if (pts(i,1) >= uknot(j)) && (pts(i,1) < uknot(j+1))
                uindex = j;
                break;
            end
            if pts(i,1) == uknot(elem_u+1);
                uindex = elem_u;
            end
        end
        for j = 1:elem_v
            if (pts(i,2) >= vknot(j)) && (pts(i,2) < vknot(j+1))
                vindex = j;
                break;
            end
            if pts(i,2) == vknot(elem_v+1);
                vindex = elem_v;
            end
        end
        num = (vindex-1)*elem_u+uindex;
        elem_index(i) = num;
    end
elseif mesh.dim == 3
    uknot = unique( mesh.uKnots );
    vknot = unique( mesh.vKnots );
    wknot = unique( mesh.wKnots );
    elem_u = length(uknot)-1;
    elem_v = length(vknot)-1;
    elem_w = length(wknot)-1;
    elem_index = zeros(1,size(pts,1));

    for i=1:size(pts,1)
        for j = 1:elem_u
            if (pts(i,1) >= uknot(j)) && (pts(i,1) < uknot(j+1))
                uindex = j;
                break;
            end
            if pts(i,1) == uknot(elem_u+1);
                uindex = elem_u;
            end
        end
        for j = 1:elem_v
            if (pts(i,2) >= vknot(j)) && (pts(i,2) < vknot(j+1))
                vindex = j;
                break;
            end
            if pts(i,2) == vknot(elem_v+1);
                vindex = elem_v;
            end
        end
        for j = 1:elem_w
            if (pts(i,3) >= wknot(j)) && (pts(i,3) < wknot(j+1))
                windex = j;
                break;
            end
            if pts(i,3) == wknot(elem_w+1);
                windex = elem_w;
            end
        end
        num = (windex-1)*elem_u*elem_v + (vindex-1)*elem_u + uindex;
        elem_index(i) = num;
    end
end

end
