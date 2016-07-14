classdef Tspikesorter < handle
    
    properties        
        % 3-dimensional array (ch x Nspk x time) containing the detected 
        % spike waveforms
        Spkwave
        
        % 4 x NSpk array containing the spatial extent of the detected 
        % spikes (Xmin,Xmax,Ymin,Ymax). It includes all channels for
        % which a threshold crossing was detected at a distance less than 3
        % channels away.
        Spkmask
        
        % 12 x NSpk array containing information about the detected spikes:
        % Spkevent(1,:) = file number (currently used only for Spiracle files)
        % Spkevent(2,:) = spike times in microseconds
        % Spkevent(3,:) = channel of maximal amplitude
        % Spkevent(4,:) = electrode on which the spike was detected
        % Spkevent(5,:) = maximal amplitude of the spike waveform
        % Spkevent(6,:) = cluster label
        % Spkevent(7,:) = not used anymore
        % Spkevent(8,:) = chi square value
        % Spkevent(9,:) = square sum value
        % Spkevent(10,:) = information on how the spike was sorted (0:
        % clustered during consensus clustering phase; 1: fitted as a
        % single spike; 2: processed as a putative overlapping spike; 3:
        % detected as an overlapping spike; 4: added overlapping spike
        % (only for SpkeventUM); 5,6: spikes with refractory period
        % violation.
        Spkevent
        
        % same as Spkevent but also contain the detected overlapping
        % spikes, reodordered according to their spike times. Therefore, 
        % SpkeventUM doesn't necessarily have the same length as Spkevent 
        % and Spkwave.
        SpkeventUM
        
        % cell array of cluster templates, defined as the average spike 
        % waveform. SpkwaveclustAve is only used for display purpose.
        SpkwaveclustAve
        
        % 2-dimensional array of correspondance between cluster and single
        % unit labels. Each row of SpkclustIDMulti correspond to 1 single 
        % unit and can contain multiple cluster label.
        SpkclustIDMulti 
        
        %Pmiss threshold applied to each cluster so that all clusters below
        %this threshold are automatically merged into one single unit.
        SpkclustPmisMulti
        
        % channel of maximal amplitude for each identified cluster
        MunitChID 
        
        % channel of maximal amplitude for each identified unit
        MclusterChID
        % same dimensions as Spkwave but Noisewave contains snippets of
        % noise waveforms, ie waveforms which didn't cross a threshold
        % corresponding to 4 x (mean absolute deviation) on any of the
        % channels
        Noisewave
        
        % cell array (1 x nb electrode) containing the spatial covaraince
        % matrix of the noise for each electrode. (computed from Noisewave)
        WhiteningMatAllSpatial
        % 3-dimensional array (NPC x ch x time). It contains the principal
        % components used for all iterations of the K-means clustering
        % algorithm.
        SpkFeatures
        
        % 2-dimensional matrix (Niteration x NSpk) containing the K-means 
        % cluster labels of each spike for every iteration of the K-means 
        % algorithm. Computed by the SpkClusterMultiCh method.
        IDXglobalAll
        
        % 2-dimensional matrix (Niteration x NSpk) containing the chi 
        % square value of how well the K-means cluster explained the spike
        % waveform at every iteration of the K-means algorithm.
        ErrorglobalAll
        
        % structure containing quality measures computed from the consensus
        % clustering. Fields are:
        % ConfusionCMat (resp. ConfusionUMat) = Pairwise probabilities of 
        % spike misclassification between pairs of clusters (resp. units),
        % computed from IDXglobalAll as the average number of times spikes
        % from cluster (or unit) #x have been co-localized in the same
        % K-means cluster as spikes from cluster (or unit) #y.
        % Note that this probabbility is computed by considering the total 
        % number of spikes in clusters #x and #y.
        % FalsePosCMat (resp. FalsePosUMat): same as ConfusionCMat but for
        % proababilities of false positive misclassification. This
        % probability os computed by considering only the number of spikes
        % contained in cluster #x (rows of FalsePosCMat).
        % FalseNegCMat (resp. FalseNegUMat) = same as FalsePosCMat but for
        % false negative misclassifications.
        % FalsePosC (resp. FalsePosU) = total rate of false positive
        % misclassifications for each cluster (resp. unit)
        % FalseNegC (resp. FalseNegU) = total rate of false negative
        % misclassifications for each cluster (resp. unit)
        Quality 
        
        % Structure containing parameters related to the file details, the 
        % detection procedure and the clustering. Fields are:
        % - fileinfo = structure with file information.
        % - detect = structure with parameters related to the detection
        % procedure.
        % - MCluster = structure with parameters related to the clustering
        % and classification phase.
        % See help in DefaultParams for mpore details on related subfields.
        Params
        
        % name of the loaded file
        filename
        % path of the loaded file
        filepath
        
        %
        Vmintra
        APevent
    end
    
    properties (Hidden = true)
        % path of the folder containing the raw data. Default to the path
        % of the specified file except for Neuralynx data for which the 
        % data are in a subfolder called with the date and time of the 
        % recording  
        datafolder        
        % 2-dimensional array (ch x time) containing the raw data. RawData
        % is band-pass filtered and then used for spike detection.
        RawData
        % vector containing the starting timestamps of the segments of raw data
        % to load succesively during the spike detection. Defined by
        % calling the getEpstartime method. Epstartime is in microseconds
        % for Neuralynx recording and in sample unit for MCD files.
        Epstartime
        % start timestamps (in microseconds) of the chunk of data from
        % which spikes are detected.
        timestamp1
        % end timestamps (in microseconds) of the chunk of data from
        % which spikes are detected.
        timestamp2
        % butterworth high-pass filter
        butterHB
        % butterworth high-pass filter
        butterHA
        % butterworth low-pass filter
        butterLB
        % butterworth low-pass filter
        butterLA  
        % voltage value of the threshold used for each channel
        detectthresh        
        
        % GUI figures, panels and menus
        SpkSortoolbox
        loadmenu
        paramsmenu
        paramspanel
        detectionmenu
        clustermenu 
        classificationmenu
        pagedetect 
        pagecluster 
        pageclassify
        dialogdetect
        dialogcluster
        dialogclassify
        pagebuttondetect 
        pagebuttoncluster 
        pagebuttonclassify
        dialogunit
        dialogClustquality
        dialogprofile
        dialogUnitquality                        
    end
    
    properties (Dependent = true)
        % down sampling factor defined from the Nyquist frequency of the
        % low-pass filter used to filter the raw data. 
        dwsamplefactor
        % 2 x Nelectrode array. Contains the # of the first and
        % last channels on each electrode.
        elconfig
        % total number of identified units
        nbUnitMultiCh
        % total number of identified clusters
        nbClusterMultiCh
    end
    methods
        function obj=Tspikesorter(Fgui)
            assignin('base','SpikeSorter',obj);
            if nargin<1
                Fgui = true;
            end
            %default values            
            obj.filename='';
            obj.filepath='';
            obj.datafolder=[];
            
            obj.InitSpkdata();                        
            obj.DefaultParams();
            
            
            obj.SpkSortoolbox='not yet created';
            obj.paramsmenu=-1;
            obj.loadmenu=-1;
            obj.paramspanel=-1;
            obj.detectionmenu=-1;
            obj.clustermenu=-1;
            obj.classificationmenu=-1;
            obj.addlistener('ObjectBeingDestroyed',@deleteobj);
            function deleteobj(src,event)
                if ishandle(obj.SpkSortoolbox)
                    delete(obj.SpkSortoolbox);
                end                          
            end
            if Fgui
                obj.SpikeSortToolbox(); 
            end
        end                
        
        %% GET METHODS
        function dwsamplefactor=get.dwsamplefactor(obj)
            dwsamplefactor=max(1,2*floor(round(obj.Params.fileinfo.samplingrate/(2*obj.Params.detect.FcHigh))/2));           
        end
        
        function elconfig=get.elconfig(obj)
            elconfig=[1;obj.Params.fileinfo.electrode(1)];
            if size(obj.Params.fileinfo.electrode,2)>1
                for el=2:size(obj.Params.fileinfo.electrode,2)
                    elconfig=[elconfig [1+sum(obj.Params.fileinfo.electrode(1:el-1));sum(obj.Params.fileinfo.electrode(1:el-1))+obj.Params.fileinfo.electrode(el)]];
                end  
            end            
        end
                
        function nbunit=get.nbUnitMultiCh(obj)
            nbunit=zeros(1,obj.Params.fileinfo.Channelcount);
            if ~isempty(obj.SpkclustIDMulti)
                [allunit,~]=find(obj.SpkclustIDMulti>0);
                nbunit=numel(unique(allunit));                 
            end            
        end 
        
        function unitID=getUnitIDMultiCh(obj,clusterID)
            unitID=zeros(1,numel(clusterID));
            if ~isempty(obj.SpkclustIDMulti)
                [nonzero,~]=find(obj.SpkclustIDMulti>0);
                nonzero=unique(nonzero);
                [unitID,~]=find(ismember(obj.SpkclustIDMulti(nonzero,:),clusterID));                
            end            
        end
        function clusterID=getClusterIDMultiCh(obj,unitID)  
            [nonzero,~]=find(obj.SpkclustIDMulti>0);
            nonzero=unique(nonzero);
            clusterID=obj.SpkclustIDMulti(nonzero(unitID),obj.SpkclustIDMulti(nonzero(unitID),:)>0);                        
        end
            
        function nbcluster=get.nbClusterMultiCh(obj)
            nbcluster = numel(unique(obj.Spkevent(6,obj.Spkevent(6,:)>0)));
        end
        
        function spkclustIDmulti=get.SpkclustIDMulti(obj)
            if ~isempty(obj.SpkclustIDMulti)
                spkclustIDmulti=obj.SpkclustIDMulti;
            else
                spkclustIDmulti=single(zeros(min(20,obj.Params.MCluster.nbmaxunits),min(20,obj.Params.MCluster.nbmaxunits)));
            end
        end
        
        %% DEFAULT AND REINITIALIZATION 
        obj = InitSpkdata(obj);
        obj = DefaultParams(obj);
                
        %% LOAD AND SAVE SPIKES
        obj=SaveSpikeData(obj,fname);
        obj=LoadSpikeData(obj,fnickname);        
        obj=SaveSortedSpikes(obj,fnickname,fcomment);            
        obj=LoadSortedSpikes(obj,fnickname,Fask);
        
        function obj=SaveSortedData(obj)
            option=[];
            option.WindowStyle='modal';
            option.Resize='on';
            option.Interpreter='none';
            filesubname=inputdlg({'file subname','comment'},'Sorted Data : file description',[1 40 ; 2 40],{'',''},option);
            if ~isempty((filesubname))
                fname=filesubname{1};
                fcomment=filesubname{2};
                if ~isempty(fname)
                    obj.SaveSortedSpikes(fname,fcomment);
                else
                    obj.SaveSortedSpikes('',fcomment);
                end
            end            
        end
        
        %% LOAD FILE PARAMS AND RAW DATA
        obj = LoadFileParams(obj,fpath);
        datafolder = getdatafolder(obj,fpath);
        Epstartime = getEpstartime(obj,length2load);
        obj = LoadNeuralynxCSC(obj,num);
        obj = LoadMcdFile(obj,num);
        obj = obj.LoadMatFile(obj,num);
        obj = LoadNSFile(obj,num);
        obj = Appendfile(obj);
        obj = resetAppendfile(obj);
        
        function obj=LoadRecording(obj,num)
            switch obj.Params.fileinfo.acqsystem
                case 'Neuralynx'
                    obj.LoadNeuralynxCSC(num);
                case 'Mcdfile'
                    obj.LoadMcdFile(num);
                case 'matfile_raw'
                    obj.LoadMatFile(num);
                case 'BlackRock'
                    obj.LoadNSFile(num);
                case 'datfile'
                    obj.LoadDATFile(num);
                otherwise
                    %add your Load method function here
                    %and modify getEpstartime method accordingly
            end            
        end
        
        
        %% FILTERING AND SPIKES DETECTION                                  
        obj = createRawDataFilter(obj);
        obj = saveSpkdataTemp(obj,fpath,num);
        obj = catSpkdataTemp(obj,fpath,nbfiles);
                
        obj = FilterRawData(obj);
        obj = executeSpkdetection(obj,num);
        
        obj = updateNoisewave(obj,rawdata,winsize,dwfactor);            
        obj = BuildNoiseCovMatrix(obj);
        spkwave = NoiseWhitening(obj,SpkWave);                
        
        function obj=executedetection(obj,Fparams)
            if  nargin<2
                Fparams=false;
            end
            if Fparams
                obj.Spkdetectdialog;                
            else
                obj.Spkdetection;
                obj.SaveSpikeData;
            end
        end        
        
        function obj=Spkdetection(obj)
            obj.InitSpkdata();
            obj.createRawDataFilter();
            datapath = obj.datafolder;
            numtemp = 0;
            fpathgroup = obj.Params.detect.filepathgroup;
            for f = 1:numel(obj.Params.detect.filepathgroup)
                obj.LoadFileParams(obj.Params.detect.filepathgroup{f});
                obj.Params.detect.filepathgroup = fpathgroup;
                obj.Epstartime = obj.getEpstartime(obj.Params.detect.length2load);
                for num = 1:size(obj.Epstartime,2)
                    numtemp = numtemp+1;
                    fprintf(['processing Data segment# ' num2str(num) ' out of ' num2str(size(obj.Epstartime,2)) ' Data segments\n']);
                    fprintf('Loading...\n');
                    obj.LoadRecording(num);
                    fprintf('Filtering...\n');
                    obj.FilterRawData();
                    fprintf('Spike Detection...\n');
                    obj.executeSpkdetection(num);
                    fprintf('Almost done...\n');
                    obj.saveSpkdataTemp(datapath,numtemp);
                end   
            end
            obj.catSpkdataTemp(datapath,numtemp);
            obj.LoadFileParams(obj.Params.detect.filepathgroup{1});
            obj.Params.detect.filepathgroup = fpathgroup;
            if ~ischar(obj.Noisewave)
                obj.BuildNoiseCovMatrix();
            end
            fprintf(['Spike detection done\n']);
        end
        
        %% SPIKES CLUSTERING AND CLASSIFICATION
        [spkwave,cleanspktt] = ClustPreProcessing(obj,chfocus,timeresol,wnbptsinit,wnbpts);
        matwspk = globalPCA(obj,matwspk,cleanspktt,nbPCmin,nbPCmax,chfocus,prevPC);
        obj = SpkClusterMultiCh(obj,FKitespace);
        obj = SpkMetaClusterMultiCh(obj,Frefine);
        obj = SpkMetaClassifyMultiCh(obj);
        obj = SpkMetaUnitQuality(obj,Funit);
        obj = ReorderClusters(obj);
        Prefractory = RefractViolations(obj,Trefractory,unitsel);
        
        %% GUI
        obj = SpikeSortToolbox(obj);
        obj = SpkClusterdialog(obj);
        obj = Spkdetectdialog(obj);
        obj = UpdateParamsPanel(obj);
        obj = InstallDialogdetect(obj);
        obj = Updatepagedetect(obj,ChSel);
        obj = InstallDialogcluster(obj);
        obj = editUnitID(obj,clustSel);
        obj = UnitRating(obj);
        obj = ShowClusterquality(obj,clustSel);
        obj = ShowUnitquality(obj,clustSel);
        obj = clustprofile(obj,Usel);
        obj = Updatepagecluster(obj,ChSel,UnitSel,Funit,Foversampled,Frealign,Fwhitened,Fsorted,Fdispstd,Fdispchauto);
        obj = InstallDialogclassify(obj);
        obj = Updatepageclassify(obj,UnitSel,Fwhitened);  
        
        %% GROUND TRUTH TEST
        [falsespk,goodclust] = GroundTruthTest(obj);               
    end    
end