function Epstartime = getEpstartime(obj,length2load)
% fills the EpStartime property with the starting timestamps of the chunks
% of raw data which will be successively processed for spike detection. 
% These segments are defined by dividing the whole recording into 
% length2load second long non-overlapping segment. 
% Eventually, EpStartime must be in the same units as required by the file 
% reader (ie in micro seconds for Neuralynx files and in samples for MCD 
% files)
    switch obj.Params.fileinfo.acqsystem
        case 'Neuralynx'
            fname=[obj.datafolder filesep 'Events.nev'];
            if isunix
                [Timestamps, TTLs] = Nlx2MatEV_v3(fname, [1 0 1 0 0], 0, 1,[]);
            else
                [Timestamps, TTLs] = Nlx2MatEV(fname, [1 0 1 0 0], 0, 1,[]);
            end
            ttlzero=find(TTLs==0);
            Epstartime = Timestamps(ttlzero(1:(end-1))+1);
            tt=1;
            while tt+2<=numel(Epstartime)
                if Epstartime(tt+2)-Epstartime(tt)<=length2load*10^6
                   Epstartime(tt+1)=[];
                else
                    tt=tt+1;
                end
            end 
        case 'Mcdfile'
            fname=[obj.filepath filesep obj.filename '.mcd']; 
            Mcdlibpath=which('nsMCDLibrary.so');
            [~] = ns_SetLibrary(Mcdlibpath);
            [~, hfile] = ns_OpenFile(fname);
            [~,info] = ns_GetFileInfo(hfile); 
            nbSamples = info.TimeSpan/info.TimeStampResolution;
            nbchunk = info.TimeSpan/length2load;
            Epstartime = 1:floor(nbSamples/nbchunk):nbSamples;
        case 'matfile_raw'
            load([obj.filepath filesep obj.filename],'fileinfo');
            Epstartime = fileinfo.HeaderParams.Epstartime;
        case 'BlackRock'
            global pepNEV;
            [tmp,SamplingRateInKHZ,nchan] = nsopen([obj.filepath filesep obj.filename '.ns5']);
            nbSamples = size(pepNEV.ns.Data.data,2);
            nbchunk = (size(pepNEV.ns.Data.data,2)/(SamplingRateInKHZ*1000))/length2load;
            Epstartime = 1:floor(nbSamples/nbchunk):nbSamples;
        case 'datfile'
            datfile = fopen([obj.filepath filesep obj.filename '.dat']);
            fseek(datfile, 0, 'eof');
            nbSamples = ftell(datfile)/(obj.Params.fileinfo.Channelcount*2);
            nbchunk = ftell(datfile)/(obj.Params.fileinfo.Channelcount*2)/(obj.Params.fileinfo.samplingrate)/length2load;
            Epstartime = 1:floor(nbSamples/nbchunk):nbSamples;
    end  
end