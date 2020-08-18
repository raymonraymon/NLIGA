function mesh = build_iga_mesh( geo )
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%  Input:
%    geo - NURBS geometry structure
%  Output:
%    mesh - mesh structure
%           mesh.dim  - dimension. 1- one dimensional(curve); 2- two dimensional(surface);
%                                  3- three dimensional(volume)
%           mesh.p/mesh.q/mesh.k - degrees in u/v/w -directions
%           mesh.uKnots/mesh.vKnots/mesh.wKnots - knot vectors in u/v/w -directions
%           mesh.nCptsU/mesh.nCptsV/mesh.nCptsW 
%                       - number of control points in u/v/w directions
%           mesh.nCpts - total number of control points          
%           mesh.coords - coordinates of control points, [x1 y1 z1 w1; x2 y2 z2 w2 ...]                       
%           mesh.nElemU/mesh.nElemV/mesh.nElemW - number ofelements in u/v/w -directions
%           mesh.nElems - total number of elements
%           mesh.nElemCpts -  total number of control points in each element
%           mesh.elNodeCnt - element node connectivity
%           mesh.elDoma - elements' parametric domain [u1, u2, v1, v2, w1, w2; ...]
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

flag = length(geo.number);
if flag == 1    % curve
    mesh.dim = 1;                     % dimension 
    mesh.p = geo.order - 1;               % degree of element
    mesh.uKnots = geo.knots;          % knot vector
    mesh.nCpts = geo.number;          % total number of control points
    mesh.coords = geo.coefs';         % coordinates of control points 
    for i = 1:3
        mesh.coords(:,i) = mesh.coords(:,i)./mesh.coords(:,4);    % transfer [xw yw zw w] to [x y z w] 
    end 
    
    [elNodeCnt, elDoma] = build_knot_connectivity( geo.knots );
    mesh.nElem = size(elNodeCnt,1);  % total number of elements
    mesh.nElemCpts = geo.order;      % total number of control points in each element
    mesh.elNodeCnt = elNodeCnt;       % element node connectivity
    mesh.elDoma = elDoma;           % elements' parametric domain [u1, u2; u1 u2...]
    
elseif flag == 2  % surface
    mesh.dim = 2;                     % dimension 
    mesh.p = geo.order(1)-1;         % degree in u-direction
    mesh.q = geo.order(2)-1;         % degree in v-direction
    mesh.uKnots = geo.knots{1};      % knot vector in u-direction
    mesh.vKnots = geo.knots{2};      % knot vector in v-direction
    mesh.nCptsU = geo.number(1);    % number of control points in u-direction
    mesh.nCptsV = geo.number(2);    % number of control points in v-direction
    mesh.nCpts  = mesh.nCptsU * mesh.nCptsV;    % total of control points
    mesh.coords = reshape(geo.coefs, 4, geo.number(1) * geo.number(2))';  % coordinates of control points
    for i = 1:3
        mesh.coords(:,i) = mesh.coords(:,i)./mesh.coords(:,4);    % transfer [xw yw zw w] to [x y z w] 
    end
    
    % separately build one dimensional connectivity along the two directions
    [elNodeCntU, elDomaU] = build_knot_connectivity( geo.knots{1} );
    [elNodeCntV, elDomaV] = build_knot_connectivity( geo.knots{2} );     
    mesh.nElemU = size(elNodeCntU,1);               % number of elemens in u direction
    mesh.nElemV = size(elNodeCntV,1);               % number of elemens in v direction
    mesh.nElems = mesh.nElemU * mesh.nElemV;      % total number of elements
    mesh.nElemCpts = geo.order(1) * geo.order(2);   % total number of nodes in one element
    mesh.elNodeCnt = zeros( mesh.nElems, mesh.nElemCpts );  % element node connectivity
    mesh.elDoma = zeros(mesh.nElems, 4);      % elements' parametric domain [u1 u2 v1 v2; u1 u2 v1 v2;...]
    count = 0;
    for j = 1:mesh.nElemV
        for i = 1:mesh.nElemU
            count = count + 1;
            for hh = 1:geo.order(2)
                for gg = 1:geo.order(1)
                    qq = (hh-1)*geo.order(1) + gg;
                    % build element-node connectivity
                    mesh.elNodeCnt(count,qq) = (elNodeCntV(j,hh)-1)*geo.number(1) + elNodeCntU(i,gg);   
                end
            end
            mesh.elDoma(count,:) = [elDomaU(i,:) elDomaV(j,:)];   % record elements' parametric domain
        end
    end    
 
elseif flag == 3  % volume
    mesh.dim = 3;                     % dimension 
    mesh.p = geo.order(1)-1;         % degree in u-direction
    mesh.q = geo.order(2)-1;         % degree in v-direction
    mesh.k = geo.order(3)-1;         % degree in w-direction
    mesh.uKnots = geo.knots{1};      % knot vector in u-direction
    mesh.vKnots = geo.knots{2};      % knot vector in v-direction
    mesh.wKnots = geo.knots{3};      % knot vector in w-direction
    mesh.nCptsU = geo.number(1);    % number of control points in u-direction
    mesh.nCptsV = geo.number(2);    % number of control points in v-direction
    mesh.nCptsW = geo.number(3);    % number of control points in w-direction
    mesh.nCpts  = mesh.nCptsU * mesh.nCptsV * mesh.nCptsW;    % total of control points
    mesh.coords = reshape(geo.coefs, 4, geo.number(1) * geo.number(2) * geo.number(3))';  % coordinates of control points
    for i = 1:3
        mesh.coords(:,i) = mesh.coords(:,i)./mesh.coords(:,4);    % transfer [xw yw zw w] to [x y z w] 
    end 
    
    % separately build one dimensional connectivity along the three directions
    [elNodeCntU, elDomaU] = build_knot_connectivity( geo.knots{1} );
    [elNodeCntV, elDomaV] = build_knot_connectivity( geo.knots{2} );   
    [elNodeCntW, elDomaW] = build_knot_connectivity( geo.knots{3} );   
    mesh.nElemU = size(elNodeCntU,1);               % number of elemens in u direction
    mesh.nElemV = size(elNodeCntV,1);               % number of elemens in v direction
    mesh.nElemW = size(elNodeCntW,1);               % number of elemens in w direction
    mesh.nElems = mesh.nElemU * mesh.nElemV * mesh.nElemW;      % total number of elements
    mesh.nElemCpts = geo.order(1) * geo.order(2) * geo.order(3);   % total number of nodes in one element
    mesh.elNodeCnt = zeros( mesh.nElems, mesh.nElemCpts );  % element node connectivity
    mesh.elDoma = zeros(mesh.nElems, 6);      % elements' parametric domain [u1 u2 v1 v2 w1 w2; u1 u2 v1 v2 w1 w2;...]    
    count = 0;
    for k = 1:mesh.nElemW
        for j = 1:mesh.nElemV
            for i = 1:mesh.nElemU
                count = count + 1;
                for ll = 1:geo.order(3)
                    for hh = 1:geo.order(2)
                        for gg = 1:geo.order(1)
                            qq = (ll-1)*geo.order(2)*geo.order(1) + (hh-1)*geo.order(1) + gg;
                            % build element-node connectivity
                            mesh.elNodeCnt(count,qq) = (elNodeCntW(k,ll)-1)*geo.number(2)*geo.number(1) + ...
                                            (elNodeCntV(j,hh)-1)*geo.number(1) + elNodeCntU(i,gg);  
                        end
                    end
                end
                % record element parametric domain
                mesh.elDoma(count,:) = [elDomaU(i,:) elDomaV(j,:) elDomaW(k,:)];  
            end
        end
    end
    
end


end

