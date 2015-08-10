function [peptides_cesarini, binds_cesarini, quant_cesarini] = import_Cesarini2013_HTP(SH2_Domains)
%Greg Koytiger
%June 2013
%Imports Cesarini Peptide Array data and builds a binding matrix

datapath = 'Input\BindingData\Cesarini2013_HTP\';

L = dir([datapath '*.seam']);
filenames = {L(:).name};

%Initializes peptide list
A = importdata([datapath filenames{1}],'\t', 1);
A.textdata(1,:) = [];
goodpep = ~strcmp(A.textdata(:,2), 'BAD_PTYR');
peptides_cesarini = A.textdata(goodpep,3);

%Initializes binding data
binds_cesarini = -1 .* ones(length(SH2_Domains), length(peptides_cesarini));
quant_cesarini =  binds_cesarini;

for i = 1:length(filenames)
    A = importdata([datapath filenames{i}],'\t', 1);
    
    binds_cesarini(strcmp(SH2_Domains(:,2), filenames{i}(1:end-5)),:) = A.data(goodpep,5)';
    quant_cesarini(strcmp(SH2_Domains(:,2), filenames{i}(1:end-5)),:) = A.data(goodpep,3)';
    
end

%Normalizes peptide
for i = 1:size(peptides_cesarini,1)
    isnormal = peptides_cesarini{i}(7) == 'Y';
   
    if(isnormal)
        peptides_cesarini{i}(7) = 'y';
        peptides_cesarini{i} = normalize_peptide(peptides_cesarini{i});
    else
        %If not regular peptide, then find closest Y to where we expect it
        %to be and "phosphorylate" it and then normalize it
        idx = find(peptides_cesarini{i} == 'Y');
        ismin = abs(8-idx) == min(abs(8-idx));
        peptides_cesarini{i}(idx(ismin)) = 'y';
        peptides_cesarini{i} = normalize_peptide(peptides_cesarini{i});
    end
     
     
end

clear isAnalyzed notChip data i A K datapath filenames goodpep idx ismin isnormal L
end