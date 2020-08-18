function curve = extract_iga_boundary ( mesh )
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%  extract the four boundaries of the two dimensional mesh
%  Input:
%    mesh - two dimensional mesh structure
%  Output:
%    curve - boundary mesh structure
%       curve.BottomCurve - corresponding to u = [0,1], v = 0
%       curve.TopCurve    - corresponding to u = [0,1], v = 1
%       curve.LeftCurve   - corresponding to u = 0, v = [0,1]
%       curve.RightCurve  - corresponding to u = 1, v = [0,1]
%           .p - degrees
%           .knots - knot vectors
%           .nCpts - total number of control points          
%           .coords - coordinates of control points, [x1 y1 z1 w1; x2 y2 z2 w2 ...]                       
%           .nElems - total number of elements
%           .nElemCpts -  total number of control points in each element
%           .elNodeCnt - total element node corresponding to the mesh
%           .elNodeCnt_local - local element node connectivity
%           .elNodeCnt_global - global element node connectivity, corresponding to the mesh
%           .elDoma - elements' parametric domain [u1, u2; ...]
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

%bottom curve of body mesh
curve.BottomCurve.p         = mesh.p;
curve.BottomCurve.knot      = mesh.uKnots;
curve.BottomCurve.nCpts     = mesh.nCptsU;
curve.BottomCurve.elNodeCnt = 1:mesh.nCptsU;
curve.BottomCurve.coords    = mesh.coords(curve.BottomCurve.elNodeCnt,:);
curve.BottomCurve.nElems    = mesh.nElemU;
curve.BottomCurve.nElemCpts = curve.BottomCurve.p+1;

[elNodeCnt, elDoma] = build_knot_connectivity( curve.BottomCurve.knot );
curve.BottomCurve.elNodeCnt_local  = elNodeCnt;
curve.BottomCurve.elNodeCnt_global = elNodeCnt;
curve.BottomCurve.elDoma = elDoma;

%top curve of body mesh
curve.TopCurve.p         = mesh.p;
curve.TopCurve.knot      = mesh.uKnots;
curve.TopCurve.nCpts	 = mesh.nCptsU;
curve.TopCurve.elNodeCnt = mesh.nCptsU*(mesh.nCptsV-1)+1:mesh.nCpts;
curve.TopCurve.coords	 = mesh.coords(curve.TopCurve.elNodeCnt,:);
curve.TopCurve.nElems    = mesh.nElemU;
curve.TopCurve.nElemCpts = curve.TopCurve.p+1;

[elNodeCnt, elDoma] = build_knot_connectivity( curve.TopCurve.knot );
curve.TopCurve.elNodeCnt_local  = elNodeCnt;
curve.TopCurve.elNodeCnt_global = elNodeCnt+(mesh.nCptsV-1)*mesh.nCptsU;
curve.TopCurve.elDoma = elDoma;

%left curve of body mesh
curve.LeftCurve.p         = mesh.q;
curve.LeftCurve.knot      = mesh.vKnots;
curve.LeftCurve.nCpts     = mesh.nCptsV;
curve.LeftCurve.elNodeCnt = 1:mesh.nCptsU:mesh.nCpts;
curve.LeftCurve.coords    = mesh.coords(curve.LeftCurve.elNodeCnt,:);
curve.LeftCurve.nElems    = mesh.nElemV;
curve.LeftCurve.nElemCpts = curve.LeftCurve.p+1;

[elNodeCnt, elDoma] = build_knot_connectivity( curve.LeftCurve.knot );
curve.LeftCurve.elNodeCnt_local  = elNodeCnt;
for j = 1:size(elNodeCnt,2)
    for i = 1:size(elNodeCnt,1)
        elNodeCnt(i,j) = elNodeCnt(i,j)+(elNodeCnt(i,j)-1)*(mesh.nCptsU-1);
    end
end
curve.LeftCurve.elNodeCnt_global = elNodeCnt;
curve.LeftCurve.elDoma = elDoma;

%right curve of body mesh
curve.RightCurve.p         = mesh.q;
curve.RightCurve.knot      = mesh.vKnots;
curve.RightCurve.nCpts     = mesh.nCptsV;
curve.RightCurve.elNodeCnt = mesh.nCptsU:mesh.nCptsU:mesh.nCpts;
curve.RightCurve.coords    = mesh.coords(curve.RightCurve.elNodeCnt,:);
curve.RightCurve.nElems    = mesh.nElemV;
curve.RightCurve.nElemCpts = curve.RightCurve.p+1;

[elNodeCnt, elDoma] = build_knot_connectivity( curve.RightCurve.knot );
curve.RightCurve.elNodeCnt_local  = elNodeCnt;
curve.RightCurve.elNodeCnt_global = elNodeCnt*mesh.nCptsU;
curve.RightCurve.elDoma = elDoma;

end