function image = loadAmOrNrrd(inputFile)
%Identifies & loads input file depending on it's file type: '.am' or '.nrrd'
%   Detailed explanation goes here
periods = strfind(inputFile,'.'); % this will look at the final period to determine filetype
filetype = char(inputFile(periods(length(periods)):length(inputFile))); 
switch filetype
    case '.am' % Load Amira file using am2mat
        amiraFile = am2mat(inputFile);
        data = amiraFile.data;
        voxel_size = amiraFile.voxel_size;
        start = amiraFile.start;
    case '.nrrd' % Load nrrd file and reformat voxel_size and start
        [data, PETmeta] = nrrdread(inputFile);
        voxel_size = PETmeta.spacedirections(2:end-1);
        voxel_size = strrep(voxel_size,',',' ');
        voxel_size = strrep(voxel_size,'(',' ');
        voxel_size = strrep(voxel_size,')',' ');
        voxel_size = textscan(voxel_size,'%f');
        voxel_size = (voxel_size{1})';
        voxel_size(voxel_size ==0) = [];

        start = PETmeta.spaceorigin(2:end-1);
        start = strrep(start,',',' ');
        start = textscan(start,'%f');
        start = (start{1})';
end 
image = struct('start', start, 'voxel_size', voxel_size, 'data', data);

return 

