function obj = DefaultParams(obj)
% reinitializes to default the parameter values contained in the
% Params structure property

    %% file info parameters
    %type of acquisition system (currently supported: 'Neuralynx','Mcd')
    obj.Params.fileinfo.acqsystem='Neuralynx';
    %sampling rate
    obj.Params.fileinfo.samplingrate=32000;
    %gain of the analog to digital converter
    obj.Params.fileinfo.ADgain=0.0001;
    % gain of the amplifier(usually = 1)
    obj.Params.fileinfo.H0gain=1;
    %input range of the amplifier. Later used to detect when signal saturates
    obj.Params.fileinfo.InputRange=2000;
    %total number of channels
    obj.Params.fileinfo.Channelcount=64;
    %physical #ID of the channels analog #ID
    obj.Params.fileinfo.Channelnum=1:64;
    %number of channels per probes (each probes will be sorted seperately)
    obj.Params.fileinfo.electrode=[32 32];
    %X positions of the channels on the array, in inter-site distance Unit
    obj.Params.fileinfo.ChannelPosX=[zeros(1,32) 1000*ones(1,32)];
    %Y positions of the channels on the array, in inter-site distance Unit
    obj.Params.fileinfo.ChannelPosY=[0:31 0:31];
    
    
    %% detection parameters
    %approximate duration of the data segments to load successively during
    %the spike detection procedure (in seconds)
    obj.Params.detect.length2load=64;
    %High cut-off of the band pass filter
    obj.Params.detect.FcHigh=4000;
    %Low cut-off of the band pass filter
    obj.Params.detect.FcLow=200;
    %Order of the Butterworth filter
    obj.Params.detect.butterOrder=6;
    %If true, we'll perform a spatial median filter
    obj.Params.detect.Fsubmedian=true;
    %If true, the raw data are interpolated 4 times with a spline 
    %interpolation 
    obj.Params.detect.Fspline=false;
    %Z-threshold for spike detection. Crossing threshold is defined as
    %zthresh x median of the absolute deviations
    obj.Params.detect.zthresh=5;
    %time window of the spike waveforms before spike time. The spike time
    %is defined as the time of the trough following threshold crossing
    obj.Params.detect.wtime1=-3;
    %time window of the spike waveforms after spike time. The spike time
    %is defined as the time of the trough following threshold crossing
    obj.Params.detect.wtime2=3;
    %AD gain used to store the spike waveforms
    obj.Params.detect.SpkADgain=32767/1000;
    %censored period in fraction of the spike waveform time window. Spike
    %times closer than the censored period are lumped together if their
    %spatial extent is less than 3 electrode sites apart.
    obj.Params.detect.Tcensoredfactor=0.3;
    %if true, we'll detect the AP in the intracellular recording (currently 
    %valid for matfile_raw and mcd recording only)
    obj.Params.detect.FGroundtruth=false;
    
    
    %% Clustering parameters
    %Maximal number of spikes which is used to perform the feature
    %extraction and the initial K-means clustering at each iteration.
    obj.Params.MCluster.TrainingMaxnbSpk=100000;
    %maximal number of PC per channel used by the K-means clustering
    obj.Params.MCluster.nbPCmax=3;
    %maximal number of PC per channel used by the K-means clustering
    obj.Params.MCluster.nbPCmin=0;
    %During the classification phase, the cluster templates are updated
    %every UpdateTime seconds
    obj.Params.MCluster.UpdateTime=300;
    %maximal range of the scaling factor of the cluster template amplitude
    obj.Params.MCluster.ProjRange=0.2;
    %maximal number of K-means iteration
    obj.Params.MCluster.NbIteration = 100;
    %maximal number of clusters identified by each iterated K-means
    obj.Params.MCluster.nbmaxunits=floor(50000^0.5);
    %Values of the minimal size threshold tested during the consensus
    %clustering phase.
    obj.Params.MCluster.minUsize=3:15;
    %Best values of the minimal size threshold, identified as the size 
    %threshold which led to the largest number of clusters.
    obj.Params.MCluster.bestminUsize=0;
    %Probability of misclassification threshold. Clusters identified by
    %consensus clustering will have no more than Pmisthreshold fraction of
    %misclassified spikes.
    obj.Params.MCluster.Pmisthreshold=0.1;
    %if true, the consensus based clusters are identified using a greedy
    %method instead of cutting the dendrogram at Pmisthreshold in one step.
    obj.Params.MCluster.Fgreedy=false;
    %The Chi square threshold is defined from the cumulative distribution
    %of chi square values. Spikes with a chi square value higher than this
    %threshold won't be used in the consensus clustering and will be fitted
    %as overlaps.
    obj.Params.MCluster.cdfChi2thresh=0.95;
    %Spikes with a chi square value higher than this threshold will be
    %conidered as noise.
    obj.Params.MCluster.cdfChi2threshnoise=0.95;
    %seed of the random number generator initialiazed before strating the
    %K-means clustering iterations
    obj.Params.MCluster.Seed = 1;
    %number of Kmeans clusters used at each Kmeans iteration
    obj.Params.MCluster.NbKcluster = 80;
    %number of Replicate for the K-means algorithm; see Matlab help
    obj.Params.MCluster.KmeansReplicates = 1;
    %online update phase option for Matlab K-means algorithm
    obj.Params.MCluster.FKmeansOnline = 'off';
    %method to define the core clusters: 'Consistency' defines core
    %clusters bas on a hierarchical clustering using the inconsistency
    %coefficient; 'Best' uses the K-clusters found for the best iteration
    %(in the least squares sense). For large datasets, prefer 'Best' to
    %'Consistency as the latter could lead to memory issues.
    obj.Params.MCluster.CoreClustersMeth = 'Consistency';
end