data_path = 'Input/';
fileToRead = 'Phosphorylation_site_dataset';
header_lines = 4;

try
    newData1 = importdata([data_path fileToRead], '\t', header_lines);
catch err
    gunzip('http://www.phosphosite.org/downloads/Phosphorylation_site_dataset.gz', data_path)
    newData1 = importdata([data_path fileToRead], '\t', header_lines);
end

% Create new variables in the base workspace from those fields.
data = newData1.data;
textdata = newData1.textdata;

headers = textdata(header_lines,:);
textdata(1:header_lines,:) = [];
data(isnan(data)) = 0;

clear fileToRead header_lines newData1 vars i

%%Finds the indexes for the information (Future proofs to changes in
%%indexing)

MOD_RSD = find(strcmp('MOD_RSD',headers));
ACC_ID = find(strcmp('ACC_ID',headers));
GENE_SYMB = find(strcmp('GENE',headers));
ORG = find(strcmp('ORGANISM',headers));
IN_DOMAIN = find(strcmp('DOMAIN',headers));
MODSITE_SEQ = find(strcmp('SITE_+/-7_AA',headers));
PROTEIN =  find(strcmp('PROTEIN',headers));

%Filters out all non-pY records 
eliminate = cellfun(@isempty, strfind(textdata(:,MOD_RSD), 'Y'));  
textdata(eliminate, :) = [];
data(eliminate,:) = [];


%Finds other species data that corresponds to residues in human and sums
%the supporting information
eliminate = false(size(textdata,1),1);
isHuman = strcmp(textdata(:,ORG), 'human');
textdata(:,MODSITE_SEQ) = upper(textdata(:,MODSITE_SEQ));
textdata(:,GENE_SYMB) = upper(textdata(:,GENE_SYMB));


for i = 1:size(textdata,1)

    if(~isHuman(i))
        %Checks if a simple data merge can be performed (Human site exists
        %with same pY sequence
        
        homolog = strcmp(textdata{i,PROTEIN}, textdata(:,PROTEIN)) & isHuman;

        tf = homolog & strcmp(textdata{i,MODSITE_SEQ}, textdata(:,MODSITE_SEQ));
             
        if ~sum(tf)==0
            try
                data(tf,:) = data(tf,:) + data(i,:);
            catch err
                data(tf,:) = data(tf,:) + repmat(data(i,:), sum(tf), 1); %Assigns support to all when ambiguity
            end
            eliminate(i) = true;
        elseif ~sum(homolog)==0
            humanid = unique(textdata(homolog,ACC_ID));
            try %Tries to get sequence data from uniprot; eliminates if deprecated or not uniprot id
                [~,seqhum] = fastaread(strcat('http://www.uniprot.org/uniprot/', humanid{:}, '.fasta'));
                [score,~,start] = swalign(textdata{i,MODSITE_SEQ}, seqhum);
                orthoseq = seqhum(start(2):start(2)+14);
                
                tf = strcmp(textdata{i,GENE_SYMB}, textdata(:,GENE_SYMB)) & isHuman &  strcmp(orthoseq, textdata(:,MODSITE_SEQ));
                
                %If conversion results in a redundant entry, sum data
                if ~sum(tf) == 0
                    data(tf,:) = data(tf,:) + data(i,:);
                    eliminate(i) = true;
                    
                elseif ~sum(homolog) == 0 && score > 0 && orthoseq(8) == 'Y'
                    %If conversion results in a new site, convert to ortholog
                    
                    textdata{i,ACC_ID} = humanid{:};
                    tmp = unique(textdata(homolog,GENE_SYMB));
                    textdata{i,GENE_SYMB} = tmp{:};
                    textdata{i,MOD_RSD} = strcat('Y', num2str(start(2)+7));
                    textdata{i,ORG} = 'human';
                    textdata{i,MODSITE_SEQ} = orthoseq;
                        
                else
                    eliminate(i) = true; %If protein never phosphorylated in humans, eliminate
                end
            catch err
                eliminate(i) = true;
            end


        end
        
        
    end
end
eliminate = eliminate | ~(cellfun('length', textdata(:,2)) == 6); %length to eliminate not-Uniprot id
textdata(eliminate,:) = []; data(eliminate,:) = [];

clear other_textdata other_data i newdata isHuman eliminate tf
    
% Converts all modified residues from sequence to their WT except for the target site
textdata(:,MODSITE_SEQ) = strrep(textdata(:,MODSITE_SEQ), '_', '-');
for i = 1:length(textdata)
    textdata{i,MODSITE_SEQ}(8) = 'y';
    textdata{i,MOD_RSD} = sscanf(textdata{i,MOD_RSD}(2:end), '%i');
end

clear i  valid_phosphosite

% Filters phosphosites that do not meet the criteria of at least two literature references or having a CST antibody against
% epitope
data(data(:,4) > 0,4) = 1;
references = sum(data(:,1:2), 2) + data(:,4);

%Only keeps data that will be useful for downstream processing
headers = [headers(GENE_SYMB) headers(ACC_ID) headers(MOD_RSD) headers(MODSITE_SEQ) headers(IN_DOMAIN)];
Phosphosites = [textdata(:, GENE_SYMB) textdata(:, ACC_ID),  textdata(:, MOD_RSD),  textdata(:, MODSITE_SEQ), textdata(:, IN_DOMAIN)];
eliminate = cellfun('isempty', Phosphosites(:,1));
references(eliminate) = [];
Phosphosites(eliminate,:) = [];

%Removes non-unique entries. Seems to occur when same gene is mapped by
%multiple Uniprot entries

[~, unique_loc] = unique(strcat(Phosphosites(:,1), Phosphosites(:,3)));

Phosphosites = Phosphosites(unique_loc,:);
references = references(unique_loc);

save('Output/PhosphositeData.mat', 'Phosphosites', 'headers', 'references');
clear textdata data loc pep_end pep_start score seq seqhum score
