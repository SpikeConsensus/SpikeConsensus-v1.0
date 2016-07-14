function obj=resetAppendfile(obj)
% resets the file path list to its first path, which corresponds to the
% loaded file
    obj.Params.detect.filepathgroup=obj.Params.detect.filepathgroup(1);
end