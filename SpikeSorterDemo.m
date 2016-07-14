function SpikeSorterDemo
%this is the sequence of methods to call to detect and sort spikes using
%the Tspikesorter object. Alternatively, one can use the GUI (by calling
%the creator with Fgui = true or the SpikeSortToolbox method once the
%object is created).
    FdispDialogs = false;
    %alternatively, one can use the GUI to perform the same commands and
    %visualize the isolated clusters. Use Fgui = true or call method 
    %SpikeSorter.SpikeSortToolbox after object creation.
    Fgui = false;
    %Creating TSpikeSorter object
    SpikeSorter = Tspikesorter(Fgui);
    assignin('base','SpikeSorter',SpikeSorter);
        
    
    %Loading file details
    fpath = '/storage/laur/Data_2/FournierJulien/SpikeSorterDemo/ctrl0031';
    SpikeSorter.LoadFileParams(fpath);
    
    %Spike detection
    SpikeSorter.executedetection(FdispDialogs);
    
    %Spike clustering (phase 1: K-means iterations)
    SpikeSorter.SpkClusterMultiCh();
    %Spike clustering (phase 2: consensus clustering)
    SpikeSorter.SpkMetaClusterMultiCh();
    %Spike clustering (phase 3: classification of outliers and overlapping 
    %spikes by adaptive template matching)
    SpikeSorter.SpkMetaClassifyMultiCh();
    
    %Cluster quality measures based on K-means iterations
    SpikeSorter.SpkMetaUnitQuality(false);
    %Reordering the clusters based on their proximity oin the dendrogram
    SpikeSorter.ReorderClusters();
    
    %validation
    % modify matrix SpikeSorter.SpkclustIDMulti such as clusters which should be
    % merged into the same putative single unit are on the same line of the
    % matrix
    
    %Unit quality measures based on K-means iterations
    SpikeSorter.SpkMetaUnitQuality(true);
    
    %Save the sorted spikes
    SpikeSorter.SaveSortedSpikes();
end