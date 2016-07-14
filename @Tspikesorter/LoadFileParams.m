function obj=LoadFileParams(obj,filepath)
% Loads the paramter file and read the path of the file. The paramter file is a
% mat file with the same name as the recording. It contains one structure
% called fileinfo with fields indicating the acquisition system, the 
% sampling rate, the number of channels, number of electrodes, the position
% of the channels. filepath is the path to the folder containing the
% recording (ie not the path to the recording). This folder must be called
% with the same name as the recording.
    if nargin<2
        if ischar(obj.filepath) && ~ischar(obj.filepath)
            filepath=uigetdir(obj.filepath);
        else
            filepath=uigetdir;
        end
    end
    
    obj.Params.detect.filepathgroup={filepath};
    
    nameposit=strfind(filepath,filesep);
    if ~isempty(nameposit)
        k=length(nameposit);
        filname=filepath(nameposit(k)+1:size(filepath,2));
    end

    obj.filepath=filepath;
    obj.filename=filname; 
    obj.Params.fileinfo=[];
    
    fname=[obj.filepath filesep obj.filename '.mat'];
    if exist(fname)
        load([obj.filepath filesep obj.filename],'fileinfo');
        obj.Params.fileinfo.acqsystem=fileinfo.acqsystem;%string, acquisition system type (use to read data later)
        obj.Params.fileinfo.samplingrate=fileinfo.samplingrate;%sampling rate
        obj.Params.fileinfo.Channelcount=fileinfo.Channelcount;%total number of channels
        obj.Params.fileinfo.Channelnum=fileinfo.Channelnum;%physical #ID of the channels analog #ID
        
        obj.Params.fileinfo.ChannelPosX=fileinfo.ChannelPosX;%X positions of the channels on the array, in inter-site distance Unit
        obj.Params.fileinfo.ChannelPosY=fileinfo.ChannelPosY;%X position of the channels  on the array in  inter-site distance
        obj.Params.fileinfo.electrode=fileinfo.electrode;%number of channels per probes (each probes will be sorted seperately)
        
        obj.Params.fileinfo.ADgain=fileinfo.ADgain;%gain of the AD converter
        obj.Params.fileinfo.H0gain=fileinfo.H0gain;%amplification gain (usually equals 1)        
        obj.Params.fileinfo.InputRange=fileinfo.InputRange;%input range of the amplifier. Later used to detect when signal saturates
        
        if isfield(fileinfo,'FGroundtruth')
            obj.Params.detect.FGroundtruth = fileinfo.FGroundtruth;
        else
            obj.Params.detect.FGroundtruth = false;
        end
    else
        fileinfo.acqsystem = 'Neuralynx';
        fileinfo.samplingrate = 32000;
        fileinfo.ADgain = 0.0001;
        fileinfo.H0gain = 1;
        fileinfo.InputRange = 2000;
        fileinfo.Channelcount = 64;
        fileinfo.electrode = [32 32];
        fileinfo.Channelnum = 1:64;  
        fileinfo.ChannelPosX = [zeros(1,32) 1000*ones(1,32)];
        fileinfo.ChannelPosY = [0:31 0:31]; 
        fileinfo.FGroundtruth = false;
        fprintf('First save in a file (same folder) a structure with the following fields : \n')
        assignin('base','fileinfo',fileinfo);
        disp(fileinfo);
    end
    
    if isequal(obj.Params.fileinfo.acqsystem,'Neuralynx')
        obj.datafolder=obj.getdatafolder(obj.filepath);
    else
        obj.datafolder=obj.filepath;
    end 
    
    %if the sampling rate is below 20kHz, we'll do a spline interpolation
    %of the recording
    if obj.Params.fileinfo.samplingrate<20000
        obj.Params.detect.Fspline=true;
        obj.Params.fileinfo.samplingrate=4*obj.Params.fileinfo.samplingrate;
    else
        obj.Params.detect.Fspline=false;
    end    
end