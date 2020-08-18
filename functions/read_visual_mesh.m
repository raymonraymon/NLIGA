function vmesh = read_visual_mesh( filename )
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% Read visual mesh file
% filename - inputed file name
% vmesh - output visualized mesh structure
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

fid = fopen(filename, 'r');
if fid == -1
    return;
end
vmesh.vertices = {};
vmesh.displacement = {};
vmesh.stress = {};
vmesh.strain = {};
vmesh.face = {};
vmesh.linmesh = {};
count = 0;
while 1
    tline = fgetl(fid);
    if tline == -1
        vmesh.vertices(count) = {vertices};
        vmesh.displacement(count) = {displacement};
        vmesh.stress(count) = {stress};
        vmesh.strain(count) = {strain};
        vmesh.face(count) = {face};
        vmesh.linmesh(count) = {linmesh};
        break;
    end
    ln = sscanf(tline,'%s',1); %
    if strncmp(ln,'STEP',4)
        if count >= 1
            vmesh.vertices(count) = {vertices};
            vmesh.displacement(count) = {displacement};
            vmesh.stress(count) = {stress};
            vmesh.strain(count) = {strain};
            vmesh.face(count) = {face};
            vmesh.linmesh(count) = {linmesh};
        end
        vertices = [];
        displacement = [];
        stress = [];
        strain = [];
        face = [];
        linmesh = [];
        count = count +1;        
    elseif strncmp(ln,'v',1)
        vertices = [vertices; sscanf(tline(2:end),'%f')'];
    elseif strncmp(ln,'d',1)
        displacement = [displacement; sscanf(tline(2:end),'%f')'];
    elseif strncmp(ln,'s',1)
        stress = [stress; sscanf(tline(2:end),'%f')'];
    elseif strncmp(ln,'t',1)
        strain = [strain; sscanf(tline(2:end),'%f')'];
    elseif strncmp(ln,'f',1)
        face = [face; sscanf(tline(2:end),'%d')'];
    elseif strncmp(ln,'l',1)
        linmesh = [linmesh; sscanf(tline(2:end),'%d')'];
    end
end


end

