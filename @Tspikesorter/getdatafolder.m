function datafolder=getdatafolder(obj,fpath)
%for Neuralynx recording, data are in a seperate folder named with date and
%time. This function gets the path to this folder.
    if isequal(obj.Params.fileinfo.acqsystem,'Neuralynx')        
        Nlxfolder=[];
        list=dir(fpath);
        for l=1:size(list,1)
            if length(list(l).name)>4
                if list(l,1).isdir==1 && ~isnan(str2double(list(l).name(1:4)))
                    Nlxfolder=list(l).name;
                end
            end
        end
        datafolder = [fpath filesep Nlxfolder];
    else
        datafolder = obj.datafolder;
    end
end