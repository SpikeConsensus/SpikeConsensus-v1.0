function obj=saveSpkdataTemp(obj,fpath,num)
% saves the detected spikes of the segment of data that has been processed 
% in a temporary folder (usually saved under the same folder as the 
% recording file).
    if ~isdir([fpath filesep 'temp'])
        mkdir([fpath filesep 'temp']);
    end
    fname=[fpath filesep 'temp' filesep 'Spkdata' num2str(num) '.mat'];
    spkevent=obj.Spkevent;
    spkeventUM=obj.SpkeventUM;
    spkwave=obj.Spkwave;
    save(fname,'spkevent','spkeventUM','spkwave','-v7.3');
    obj.Spkevent=[];
    obj.SpkeventUM=[];
    obj.Spkwave=[];
end