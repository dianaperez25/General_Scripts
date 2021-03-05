function [LH_files, RH_files, middle_files] = makingdatafiles(dataLoc, fileName, LH, middle, RH, type)
%% Script to make list of paths to files
%type = type of map, 1 = variant map, 2 = overlap map

LH_files = [];
middle_files = [];
RH_files = [];

if type == 1
    
    for x = 1:length(LH)
        file = [dataLoc num2str(LH(x)) fileName];
        LH_files{x,1} = file;
    end
    for x = 1:length(middle)
        file = [dataLoc num2str(middle(x)) fileName];
        middle_files{x,1} = file;
    end
    for x = 1:length(RH)
        file = [dataLoc num2str(RH(x)) fileName];
        RH_files{x,1} = file;
    end
    
elseif type == 2
    
   for x = 1:length(LH)
        file = [dataLoc '/Left/' num2str(LH(x)) fileName];
        LH_files{x,1} = file;
    end
    for x = 1:length(middle)
        file = [dataLoc '/Middle/' num2str(middle(x)) fileName];
        middle_files{x,1} = file;
    end
    for x = 1:length(RH)
        file = [dataLoc '/Right/' num2str(RH(x)) fileName];
        RH_files{x,1} = file;
    end
    
end
