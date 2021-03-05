clear all

%allSubs = xlsread('RESTRICTED_HCP');
allSubs = table2cell(readtable('RESTRICTED_HCP.xls'));
goodSubs = table2cell(readtable('HCP_GoodSubs.xlsx'));
goodSubs_col = cell2mat(goodSubs(1:385,1));
goodSubs_all = {};
count = 1;
for x = 1:size(allSubs,1)
    sub = cell2mat(allSubs(x,1));
    if ismember(sub, goodSubs_col)
        goodSubs_all(count,:) = allSubs(x,:);
        count = count + 1;
    end 
end 