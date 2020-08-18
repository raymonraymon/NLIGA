function plot_nurbs_basis( nrb )
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%  plot NURBS basis functions
%  nrb - nurbs structure
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

if length(nrb.number) == 1      % curve
    xx = linspace(nrb.knots(nrb.order),nrb.knots(nrb.number+1),1000);
    [B, N]  = nrbbasisfun (xx, nrb);
    y = cell(1,nrb.number);
    x = cell(1,nrb.number);
    for i = 1:nrb.number
        for j =1:size(N,1)
            for k = 1:nrb.order
                if N(j,k) == i
                    y{i} = [y{i},B(j,k)];
                    x{i} = [x{i},xx(j)];
                end

            end
        end
    end
    figure;
    for i = 1:length(y)
        hold on;
        plot(x{i},y{i},'-');
    end
elseif length(nrb.number) == 2  % surface
    crv = cell(1,2);
    coefs = nrb.coefs(:,:,1);
    knots = nrb.knots{1};
    crv{1} = nrbmak(coefs,knots);
    coefs = nrb.coefs(:,1,:);
    coefs = reshape(coefs,[4,nrb.number(2)] );
    knots = nrb.knots{2};
    crv{2} = nrbmak(coefs,knots);
    for ll = 1:2
        xx = linspace(crv{ll}.knots(crv{ll}.order),crv{ll}.knots(crv{ll}.number+1),1000);
        [B, N]  = nrbbasisfun (xx, crv{ll});
        y = cell(1,crv{ll}.number);
        x = cell(1,crv{ll}.number);
        for i = 1:crv{ll}.number
            for j =1:size(N,1)
                for k = 1:crv{ll}.order
                    if N(j,k) == i
                        y{i} = [y{i},B(j,k)];
                        x{i} = [x{i},xx(j)];
                    end

                end
            end
        end
        subplot(2,1,ll);
        for i = 1:length(y)
            hold on;
            plot(x{i},y{i},'-');
        end
        hold off;
    end
elseif length(nrb.number) == 3  % volume
    crv = cell(1,3);
    coefs = nrb.coefs(:,:,1,1);
    knots = nrb.knots{1};
    crv{1} = nrbmak(coefs,knots);
    coefs = nrb.coefs(:,1,:,1);
    coefs = reshape(coefs,[4,nrb.number(2)] );
    knots = nrb.knots{2};
    crv{2} = nrbmak(coefs,knots);
    coefs = nrb.coefs(:,1,1,:);
    coefs = reshape(coefs,[4,nrb.number(3)] );
    knots = nrb.knots{3};
    crv{3} = nrbmak(coefs,knots);

    for ll = 1:3
        xx = linspace(crv{ll}.knots(crv{ll}.order),crv{ll}.knots(crv{ll}.number+1),1000);
        [B, N]  = nrbbasisfun (xx, crv{ll});
        y = cell(1,crv{ll}.number);
        x = cell(1,crv{ll}.number);
        for i = 1:crv{ll}.number
            for j =1:size(N,1)
                for k = 1:crv{ll}.order
                    if N(j,k) == i
                        y{i} = [y{i},B(j,k)];
                        x{i} = [x{i},xx(j)];
                    end

                end
            end
        end
        subplot(3,1,ll);
        for i = 1:length(y)
            hold on;
            plot(x{i},y{i},'-');
        end
        hold off;
    end
end

