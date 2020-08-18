function [R, ders] = nurbs_derivatives( gauPts, geo, mesh )
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%  Calculate derivatives of a given point
%  Input:
%    gauPts - parameter points
%    geo - geometry
%    mesh - mesh structure of the geo
%  Output:
%    R - basis functions
%    ders - derivatives
% IMPORTANT: If the C codes can't be compiled successfully, please 
%            follow the instructions below
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

num = length(gauPts);
if num == 1
    [R, dRdxi] = nrbNurbs1DBasisDerivs (gauPts,mesh.p,mesh.uKnots,mesh.coords(:,4)');
    ders = dRdxi;
elseif num == 2
    [R, deru, derv] = nrbNurbs2DBasisDerivs(gauPts, mesh.p,mesh.q,mesh.uKnots,mesh.vKnots,mesh.coords(:,4)');
    ders = [deru; derv];
elseif num == 3
    [R, deru, derv, derw] = nrbNurbs3DBasisDerivs(gauPts,mesh.p,mesh.q,mesh.k, mesh.uKnots,mesh.vKnots,mesh.wKnots, mesh.coords(:,4)');    
    ders = [deru; derv; derw];
end

%%  INSTRUCTION: Please comment the above codes If the C codes can't be compiled successfully and
%  uncomment the following codes for calculation of derivatives using
%***************************************************************
% num = length(gauPts);
% if num == 2
%     [deru, derv]  = nrbbasisfunder ({gauPts(1), gauPts(2)}, geo);
%     ders = [deru; derv];
%     R = nrbbasisfun ({gauPts(1), gauPts(2)}, geo);
% elseif num == 3
%     [deru, derv, derw]  = nrbbasisfunder ({gauPts(1), gauPts(2), gauPts(3)}, geo);
%     ders = [deru; derv; derw];
%     R = nrbbasisfun ({gauPts(1), gauPts(2), gauPts(3)}, geo);
% end
%***************************************************************




end

