function obj = InitSpkdata(obj)
% reinitialize properties of the spikesorter object related to the sorted 
% spikes  
    obj.Spkevent=[];
    obj.Spkwave=[];
    obj.Spkmask=[];
    obj.SpkeventUM=[];
    obj.IDXglobalAll=[];
    obj.ErrorglobalAll=[];
    obj.SpkclustIDMulti=[];
    obj.Quality=[];
    
    obj.MunitChID=[];
    obj.MclusterChID=[];
    obj.Noisewave=[];
    obj.WhiteningMatAllSpatial=[];
    obj.SpkFeatures=[];
    
    obj.APevent=[];
end