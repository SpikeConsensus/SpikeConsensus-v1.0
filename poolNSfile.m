function poolNSfile(fnamelist)
global pepNEV;
RawData = [];
for f = 1:numel(fnamelist)    
    fname = fnamelist{f};
    [~,~,nchan] = nsopen(fname);
    RawData = [RawData nchan];
end