function [peptides_ltp, binds_ltp] = import_Cesarini2013_LTP(SH2_Domains)

%Greg Koytiger
%Aug 2013
%Imports curated low throughput data from Cesarini 2013 paper

datapath = 'Input\BindingData\';
[~,~,raw] = xlsread([datapath 'Cesarini2013_LTP.xls']);

raw(1,:) = [];
overlap = ismember(raw(:,2), SH2_Domains(:,3));
map = unique(raw(~overlap,2)); map = [map, cell(size(map))];
peptides_ltp = unique(raw(:,9));

for i = 1:length(map)
    try
        [hdr,~] = fastaread(['http://www.uniprot.org/uniprot/' map{i,1} '.fasta']);
        idx = strfind(hdr,'GN='); idx = idx + 3; sidx = strfind(hdr,' ');
        idxend= sidx(find(sidx > idx, 1, 'first')); %finds the first space following gene
        map{i,2} = upper(hdr(idx:idxend-1));
    catch err
        %Changed some of the non-sprot entries to gene name, enters them
        %here
        map{i,2} = map{i,1};
    end
end

overlap2 = ismember(map(:,2), SH2_Domains(:,1));
map(~overlap2,:) = [];

binds_ltp = -1 * ones(length(SH2_Domains), length(peptides_ltp)); %Initializes list as untested

for i = 1:length(raw)
    if(overlap(i))
        idx = strcmp(SH2_Domains(:,3), raw{i,2}); idx = find(idx,1,'last'); %find protein, assigns multi-domain data to tandem
        binds_ltp(idx, strcmp(peptides_ltp, raw{i,9})) = 2;
    elseif (ismember(raw{i,2}, map(:,1)))
        idx = strcmp(SH2_Domains(:,1), raw{i,2}); idx = find(idx,1,'last');
        binds_ltp(idx, strcmp(peptides_ltp, raw{i,9})) = 2;
    end
end

%Normalizes peptide
for i = 1:size(peptides_ltp,1)
    isnormal = peptides_ltp{i}(7) == 'Y';
   
    if(isnormal)
        peptides_ltp{i}(7) = 'y';
        peptides_ltp{i} = normalize_peptide(peptides_ltp{i});
    else
        %If not regular peptide, then find closest Y to where we expect it
        %to be and "phosphorylate" it and then normalize it
        idx = find(peptides_ltp{i} == 'Y');
        ismin = abs(8-idx) == min(abs(8-idx));
        peptides_ltp{i}(idx(ismin)) = 'y';
        peptides_ltp{i} = normalize_peptide(peptides_ltp{i});
    end
     
     
end

clear err hdr i idx idxend isnormal map overlap overlap2 raw sidx
end