function [peptides_nash, binds_nash] = import_Nash2012(SH2_Domains)
datapath = 'Input\BindingData\';

[~,~, raw] = xlsread([datapath 'Nash2012.xls']);

binds_nash = zeros(size(SH2_Domains(:,1),1),size(raw(:,1),1));
peptides_nash = raw(:,2);
notmatch = cell(500,1);
l = 1;

map = {'BRDG1', 'STAP1';
    'GADS', 'GRAP2';
    'BRK', 'PTK6';
    'MIST', 'CLNK';
    'PI3K1_C', 'PIK3R1-C';
    'PI3K1_N', 'PIK3R1-N';
    'PLCG2_C', 'PLCG2-C';
    'PLCG2_N', 'PLCG2-N';
    'PTPN11_C', 'PTPN11-C';
    'PTPN11_N', 'PTPN11-N';
    'RASA1_C', 'RASA1-C';
    'SHIP1', 'INPP5D';
    'SLNK', 'SH2D6';
    'SYK_C', 'SYK-C';
    'Sh2b', 'SH2B1';
    'YES', 'YES1';
    'ZAP70_N', 'ZAP70-N'};
    
    
    
    
for i = 1:size(raw,1)
    peptides_nash{i} = [peptides_nash{i}(1:4) 'y' peptides_nash{i}(7:end)];
    if(~isnan(raw{i,7}))
        if(~isempty(strfind(raw{i,7}, ',')))
            bDomains =  regexp(raw{i,7}, ',', 'split');
        else
            bDomains = raw(i,7);
        end
    else
        bDomains = cell(0);
    end
    
    if(~isnan(raw{i,8}))
        if(~isempty(strfind(raw{i,8}, ',')))
            wDomains = regexp(raw{i,8}, ',', 'split');
        else
            wDomains = raw(i,8);
        end
    else
        wDomains = cell(0);
    end
    
    if(~isempty(bDomains))
        for j = 1:length(bDomains)
            [tf, loc] = ismember(bDomains{j}, map(:,1));
            if(tf)
                bDomains{j} = map{loc,2};
            end
            idx = strcmp(bDomains{j}, SH2_Domains(:,2));
            if(~isempty(find(idx, 1)))
                binds_nash(idx,i) = 2;
            else
                notmatch{l} = bDomains{j};
                l = l+1;
            end
        end
    end
    if(~isempty(wDomains))
        for j = 1:length(wDomains)
            [tf, loc] = ismember(wDomains{j}, map(:,1));
            if(tf)
                wDomains{j} = map{loc,2};
            end
            idx = strcmp(wDomains{j}, SH2_Domains(:,2));
            if(~isempty(find(idx, 1)))
                binds_nash(idx,i) = 1;
            else
                notmatch{l} = wDomains{j};
                l = l+1;
            end
        end
    end
end

for i = 1:size(SH2_Domains,1)
    if(sum(binds_nash(i,:)) == 0)
        binds_nash(i,:) = ones(1,size(binds_nash,2)) * -1;
    end
end

%Filters Array non-specific interactions
isNS = strcmp(raw(:,6), 'Y');
peptides_nash(isNS) = []; binds_nash(:,isNS) = [];

for i = 1:size(peptides_nash,1)
    peptides_nash{i} = ['---',peptides_nash{i},'-'];
end

clear notmatch wDomains bDomains i j l raw tf loc idx isNS map

