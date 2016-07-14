function obj=catSpkdataTemp(obj,fpath,nbfiles)
% concatenates the spike data that were saved by saveSpkdataTemp during the
% detection procedure.
    if ~isdir([fpath filesep 'temp'])
        error('temp folder for spike data doesn''t exist');
    end
    obj.Spkevent=[];
    obj.SpkeventUM=[];
    obj.Spkwave=[];
    for num=1:nbfiles
        fname=[fpath filesep 'temp' filesep 'Spkdata' num2str(num) '.mat'];
        load(fname,'spkevent');
        obj.Spkevent=cat(2,obj.Spkevent,spkevent);
        load(fname,'spkeventUM');
        obj.SpkeventUM=cat(2,obj.SpkeventUM,spkeventUM);
    end
    
    oldNspike = 0;    
    for num=1:nbfiles
        fname=[fpath filesep 'temp' filesep 'Spkdata' num2str(num) '.mat'];
        load(fname,'spkwave');
        newNspike = size(spkwave,2);
        if num == 1
            obj.Spkwave = zeros(obj.Params.fileinfo.Channelcount, size(obj.Spkevent,2), size(spkwave,3), 'int16');
        end
        obj.Spkwave(:, oldNspike+1:oldNspike+newNspike,:) = spkwave;
        oldNspike = oldNspike + newNspike;
%         obj.Spkwave=cat(2,obj.Spkwave,spkwave);        
    end
    rmdir([fpath filesep 'temp'],'s');
end