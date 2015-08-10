function [peptides_jones, binds_jones, quant_jones] = import_Jones2012(SH2_Domains)
datapath = 'Input\BindingData\';
[~,~,raw] = xlsread([datapath 'Jones2012.xlsx']);

raw(1,:) = []; 
raw(1396:end,:) = [];
raw(:,4) = strrep(raw(:,4), 'pY', 'y');
raw(:,4) = strrep(raw(:,4), 'd', '');

map = {'SH2B', 'SH2B1'};

%%Reorganizes into matrix with three categories -1 (not analyzed), 0
%%(not bound) and 2 (bound)
peptides_jones = unique(raw(:,4)); peptides_jones = strrep(peptides_jones, ' ', ''); %removes spaces
binds_jones = ones(size(SH2_Domains,1), size(peptides_jones,1));
binds_jones = binds_jones .* -1; %Initializes as not analyzed
quant_jones = ones(size(SH2_Domains,1), size(peptides_jones,1)); 
quant_jones = quant_jones .* -1; %Initializes as not analyzed
%max_jones = quant_jones;

for i = 1:size(raw,1)
    x = strcmp(SH2_Domains(:,2), raw{i,6}); y = strcmp(peptides_jones, raw{i,4});
    if(~isempty(find(x, 1)))
        binds_jones(x,y) = 2;
        quant_jones(x,y) = raw{i,7}*1000; %Stores the kD
        %max_jones(x,y) = raw{i,9}; %Stores maximal polarization
    end
end

Affinity_threshold = 2500;

isAnalyzed = sum(binds_jones,2) ~= -size(binds_jones,2); %finds domains which have binding events
quant_jones(repmat(isAnalyzed, 1, size(binds_jones,2)) & binds_jones == -1) = 0;
%max_jones(repmat(isAnalyzed, 1, size(binds_jones,2)) & binds_jones == -1) = 0;

binds_jones(repmat(isAnalyzed, 1, size(binds_jones,2)) & binds_jones == -1) = 0;
binds_jones(quant_jones > Affinity_threshold) = 1;

for i = 1:size(peptides_jones,1)
    peptides_jones{i} = normalize_peptide(peptides_jones{i});
end

clear Affinity_threshold i isAnalyzed map raw x y
end
