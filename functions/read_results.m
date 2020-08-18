function vmesh = read_results( filename, etype, ncpts, nelems, ngp )
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% read the btained results in elastoplasticity
% Input: filename for saving,
%        etype - element type
%        ncpts - number of control points
%        nelems - number of elements
%        ngp - number of integration points
% Output: vmesh for visualization
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

fid = fopen(filename, 'r');
if fid == -1
    return;
end
if etype == 10
    dof = 2;  
    ns = 3;  % number of stress components
elseif etype == 20
    dof = 3;
    ns = 6;
end

vmesh.displacement = {};
vmesh.stress = {};

while ~feof(fid)
    tline = fgetl(fid);
    ln = sscanf(tline,'%s',1); %
    if strncmp(ln,'TIME',4)
        displacement = zeros(ncpts, dof);
        stress = zeros(ngp*nelems, ns );
        while 1
            tline = fgetl(fid);
            ln = sscanf(tline,'%s',2); %
            if strncmp(ln,'Node',4)  % read node displacement 
                for i = 1:ncpts
                    tline = fgetl(fid);
                    tmp = sscanf(tline,'%f')';
                    displacement(i,:) = tmp(2:end);
                end
            elseif strncmp(ln,'ElementStress',13)  % read element stress 
                count = 0;
                while 1
                    tline = fgetl(fid);
                    ln = sscanf(tline,'%s',2); %
                    if strncmp(ln,'Elem',4)
                        count = count+1;
                        for i = 1:ngp
                            tline = fgetl(fid);
%                             stress(i,(count-1)*ns+1:count*ns) = sscanf(tline,'%f')';
                            stress((count-1)*ngp+i,1:ns) = sscanf(tline,'%f')';
                        end
                    end
                    if count == nelems
                        vmesh.stress = [vmesh.stress, stress];
                        vmesh.displacement = [vmesh.displacement, displacement];
                        break;
                    end
                end  
                if count == nelems
                    break;
                end
            end        
        end

    end
end


end

