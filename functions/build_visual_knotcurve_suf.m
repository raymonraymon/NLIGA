function [ kntcrv ] = build_visual_knotcurve_suf( uknots, vknots, num )
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%  Build parametric knots lines structure for visualization
%  Input:
%    uknots, surface knots vector in u direction
%    vknots, surface knots vector in v direction
%    num,    segment number for refinement of knot lines
%  Output:
%    kntcrv - all knot lines
%           kntcrv.linpts - segment vertices
%           kntcrv.linmesh - segment-node connectivity
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

uknot = unique( uknots);
vknot = unique( vknots);

unum = length(uknot);
vnum = length(vknot);
linpts = zeros(2,(unum+vnum)*num);
linmesh = zeros((unum+vnum),num);
count = 1;
offset = 0;
for i = 1:length(vknot)
    linpts(:,(count-1)*num+1:count*num) = [linspace(0+offset,1-offset,num); ones(1,num)*vknot(i)];
    linmesh(count,:) = (count-1)*num+(1:num);
    count = count+1;
end

for i = 1:length(uknot)
    linpts(:,(count-1)*num+1:count*num) = [ones(1,num)*uknot(i); linspace(0+offset,1-offset,num)];
    linmesh(count,:) = (count-1)*num+(1:num);
    count = count+1;
end

kntcrv.linpts = linpts';
kntcrv.linmesh = linmesh;

end

