function filename = get_output_file_name( filename )
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% get the address for file saving
% the output file is saved int the folder \output
%  ---------------------------------------
%  Please feel free to contact us with any questions! 
%  - Xiaoxiao Du, Beihang University
%  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

fname = fileparts(mfilename('fullpath')); % get current path
index_dir = strfind(fname,'\');
str_temp = fname(1:index_dir(end));
filename = [str_temp,'output\', filename, '.msh'];

end

