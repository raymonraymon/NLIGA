function [ stress, alpha, ep, Dtan ] = material_plasticity( D, dim, mat, deps, count )
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Main frame for nonlinear isogeometric analysis
% Inputs:
%   D 		  = elastic stiffness matrix
%   dim       = 2d or 3d
%   mat       = [index, E, nu, sigmay0, H, beta];
%   deps      = strain increment 
%   stressN   = [s11, s22, s33, t12, t23, t13];
%   alphaN    = back stress [a11, a22, a33, a12, a23, a13]; 
%   epN		  = effective plastic strain
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

global XQ  SIGMA; 

two3    = 2/3;
stwo3   = sqrt(two3);

E       = mat(2);
nu      = mat(3);
sigmay0 = mat(4);
H       = mat(5);
beta    = mat(6);
mu      = E/2/(1+nu);

ftol    = sigmay0*1E-6;%tolerance for yield

if dim == 2    % plane element
    
    ddeps = [deps(1,1) deps(2,2) deps(1,2)+deps(2,1)]';
    stressN = SIGMA(1:4,count);
    alphaN = XQ(1:4,count);
    epN = XQ(5,count);

    Iden = [1 1  0 1]'; 
    D_deps = zeros(4,1);
    D_deps(1:3) = D*ddeps;
    trstress = stressN + D_deps;
    trstress(4) = nu * (trstress(1) + trstress(2));%plane srain
    tr = trstress(1)+trstress(2)+trstress(4);      %trial stress
    shstress = trstress-alphaN-tr*Iden/3;          %shifted stress
    norm = sqrt(shstress(1)^2+shstress(2)^2+shstress(4)^2+2*shstress(3)^2);
    fyld = norm - stwo3*(sigmay0+(1-beta)*H*epN);       %yield function
    if fyld < ftol
        stress = trstress(1:4,:);               %update stress
        alpha = alphaN;                         %update back stress
        ep = epN;                               %update effective plastic strain
        Dtan = D;                               %elastic
    else
        Dtan1 = D(1,1)*eye(4);
        Dtan1(1:3,1:3) = D;
        trr_1 = stressN(1)+stressN(2)+stressN(4);
        stressr_1 = stressN-alphaN-trr_1*Iden/3;
        normr_1 = sqrt(stressr_1(1)^2+stressr_1(2)^2+stressr_1(4)^2+2*stressr_1(3)^2);
        fyldr_1 = normr_1 - stwo3*(sigmay0+(1-beta)*H*epN); 
        if fyldr_1 < -ftol
            R = (norm/stwo3-sigmay0)/(norm/stwo3-normr_1/stwo3);
        else
            R = 1;
        end
        SGTOT = stressN + (1-R)*D_deps;
        STRES = R * D_deps; 
        trSGTOT = SGTOT(1)+SGTOT(2)+SGTOT(4);
        DEVIA = SGTOT-alphaN-trSGTOT*Iden/3;
        STEFF = sqrt(0.5*DEVIA(1)^2+0.5*DEVIA(2)^2+0.5*DEVIA(4)^2+DEVIA(3)^2);
        YIELD = sqrt(3)*STEFF;
        AVECT = sqrt(3)* [DEVIA(1)/(2*STEFF) DEVIA(2)/(2*STEFF) DEVIA(3)/STEFF DEVIA(4)/(2*STEFF)]';
        Dtan1 = D(1,1)*eye(4);
        Dtan1(1:3,1:3) = D;
        Dtan1(1,4) = D(1,2);
        Dtan1(2,4) = D(1,2);
        Dtan1(4,1) = D(1,2);
        Dtan1(4,2) = D(1,2);
        d_D  = Dtan1*AVECT;
        AGASH = AVECT' * STRES;
        ABETA = 1/(H+d_D'*AVECT);
        DLAMD = AGASH*ABETA;
        BGASH = AVECT' * SGTOT;
        stress_r = SGTOT + STRES - DLAMD * d_D;
        ep1 = epN + DLAMD*BGASH/YIELD;
        tr_r = stress_r(1)+stress_r(2)+stress_r(4);
        DEVIA = stress_r-alphaN- tr_r*Iden/3;
        STEFF = sqrt(0.5*DEVIA(1)^2+0.5*DEVIA(2)^2+0.5*DEVIA(4)^2+DEVIA(3)^2);
        YIELD = sqrt(3)*STEFF;
        CURYS = sigmay0+(1-beta)*H*ep1;
        BRING = CURYS/YIELD;
        stress = BRING*stress_r;
        alpha = alphaN;
        ep = BRING*ep1;
        
        DT = Dtan1-(d_D*d_D'/(H+d_D'*AVECT));
        Dtan = DT(1:3,1:3);
    end
elseif dim == 3    % solid element
    
    ddeps = [deps(1,1) deps(2,2) deps(3,3)...
        deps(1,2)+deps(2,1) deps(2,3)+deps(3,2) deps(1,3)+deps(3,1)]';
    stressN = SIGMA(1:6,count);
    alphaN = XQ(1:6,count);
    epN = XQ(7,count);

    Iden = [1 1 1 0 0 0]';                      %unit tensor
                          
    trstress = stressN + D*ddeps;                %trial stress
    tr = trstress(1)+trstress(2)+trstress(3);   %trace of stress
    shstress = trstress-alphaN-tr*Iden/3;       %shifted stress
    norm = sqrt(shstress(1)^2+shstress(2)^2+shstress(3)^2+...
        2*(shstress(4)^2+shstress(5)^2+shstress(6)^2));
    fyld = norm - stwo3*(sigmay0+(1-beta)*H*epN);    %yield function
    
    if fyld < ftol
        stress = trstress;                      %update stress
        alpha = alphaN;                         %update back stress
        ep = epN;                               %update effective plastic strain
        Dtan = D;                                           %elastic
    else
        gamma = fyld / (2*mu+two3*H);           %consistency parameter
        N = shstress / norm;                    %unit deviatoric vector
        stress = trstress - 2*mu*gamma*N;       %update stress
        alpha = alphaN + two3*beta*H*gamma*N;   %update back stress
        ep = epN + stwo3*gamma;                 %update effective plastic strain
        
        
        var1 = 4*mu^2/(2*mu+two3*H);
        var2 = 4*mu^2*gamma/norm;                           %coefficients
        Dtan = D - (var1-var2)*N*N' + var2*Iden*Iden'/3;    %tangent stiffness
        Dtan(1,1) = Dtan(1,1) - var2;                       %contr. from 4th-order I
        Dtan(2,2) = Dtan(2,2) - var2;
        Dtan(3,3) = Dtan(3,3) - var2;
        Dtan(4,4) = Dtan(4,4) - .5*var2;
        Dtan(5,5) = Dtan(5,5) - .5*var2;
        Dtan(6,6) = Dtan(6,6) - .5*var2;
%         a=[3/2*shstress(1)/norm 3/2*shstress(2)/norm 3/2*shstress(3)/norm 3*shstress(4)/norm 3*shstress(4)/norm 3*shstress(4)/norm]';
%         Dtan=D-(D*a*a'*D/(H+a'*D*a));
        
    end

end

