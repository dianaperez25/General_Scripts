function [path2] = convert_path(path1)
if ispc
    path2 = replace(path1, '\', '/');
elseif ismac
    path2 = replace(path1, '/', '\');
end
end