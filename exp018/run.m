restoredefaultpath
addpath(genpath('../../neutral_common'))
addpath(genpath('.'))

path_to_gammanc='/home/nfs/z3439823/mymatlab/neutral_common/gamma_n/gamma_JackettMcDougall96';
[s,ct,p,lon,lat]=get_input_jackett96(path_to_gammanc);

[s,ct,p,lon,lat]=crop_to_gamma_n(s,ct,p,lon,lat);
[s,ct]=only_keep_largest_region(s,ct); % remove everything except largest region
    
    
%field = ew99_modified(s,ct,p,lon,lat,false);
%field = v2_ew99_modified(s,ct,p,lon,lat,false);

%field=ew99_modified(s,ct,p,lon,lat,true); % last arg is true-> modified
field=v2_ew99_modified(s,ct,p,lon,lat,true); % last arg is true-> modified