function [v] = ZZ_read_reshape_tiff(fileName, volume_info)
% read a single tiff file, reshape into a N-d matrix based on volume_info
    d = loadtiff(fileName);
    disp('original shape of tiff is:');
    disp(size(d));
    v = reshape(d, volume_info.xRes, volume_info.yRes, volume_info.nChannel,...
            volume_info.zStacks, []);
    disp('reshaped as (xyczt):');
    disp(size(v));
    
    if ~isa(v, volume_info.dtype)
        disp('Warning: data type is not what specified! May cause clipping!');
        disp(strcat('Converted to', volume_info.dtype));
        v = cast(v, volume_info.dtype);
    end
    
end

