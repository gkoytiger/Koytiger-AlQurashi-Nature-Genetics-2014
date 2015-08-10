function [peptides_macbeath, binds_macbeath, quant_macbeath] = import_MacBeath2013(SH2_Domains)

%Greg Koytiger
%Sep 2012
%Imports SH2 Domain array info
datapath = 'Input\BindingData\';


Affinity_threshold = 1000; %Sets the affinity threshold for calling bound

load([datapath 'MacBeath2013.mat']);
[peptides_macbeath, m] = unique(Peptide_list(:,4));

quant_macbeath = ones(size(SH2_Domains,1), size(peptides_macbeath,1)); 
quant_macbeath = quant_macbeath .* -1; %Initializes as not analyzed

for i = 1:size(Domain_list,1)
    x = strcmp(SH2_Domains(:,2), Domain_list(i,3));
    if(sum(x) == 1)
        quant_macbeath(x,:) = matrix_Kd(i,m);
    end
    
end


binds_macbeath = double(quant_macbeath > 0 & quant_macbeath < Affinity_threshold)*2;
binds_macbeath(quant_macbeath > Affinity_threshold) = 1;
binds_macbeath(quant_macbeath == -1) = -1;

for i = 1:size(peptides_macbeath,1)
    peptides_macbeath{i} = peptides_macbeath{i}(3:end);
end

clear Domain_list Peptide_list i m matrix_Kd x Affinity_threshold