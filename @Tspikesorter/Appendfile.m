function obj=Appendfile(obj)
%simply append the path of another file into Params.detect.filepathgroup.
%All file path will be processed as a single file. The spikedata file will
%be saved at the first path of the list. The sorted spike data will be
%saved at all path of the list.
    fpath=uigetdir(obj.filepath);
    obj.Params.detect.filepathgroup=[obj.Params.detect.filepathgroup {fpath}];
end