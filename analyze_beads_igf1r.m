%Greg Koytiger Sept 2013
%Process Bead Images of PI3K and IGF1R Binding

im_path = 'N:\SYSTEMS BIOLOGY SHARED\Greg\2014_05_08\2014_05_08_anti_rfp__2014-05-08T16_24_45-Measurement1\Images\';
L = dir([im_path '*.tiff']);
filenames = {L(:).name};

rows = zeros(numel(filenames),1);
columns = rows; fields = rows; channels  = rows;

for i = 1:numel(filenames)
    rows(i) = sscanf(filenames{i}(2:3),'%u');
    columns(i) = sscanf(filenames{i}(5:6),'%u');
    fields(i) = sscanf(filenames{i}(8:9),'%u');
    channels(i) = sscanf(filenames{i}(16),'%u');
end

data = zeros(length(unique(rows)), length(unique(columns)), length(unique(fields)), length(unique(channels)));
r0 = min(rows); c0 = min(columns);

for i = 1:numel(filenames)
        I = imread([im_path filenames{i}]);
        data(rows(i)+1-r0, columns(i)+1-c0, fields(i), channels(i)) = median(I(:));
end

unique_chan = unique(channels);
for i = 1:numel(unique_chan)
    assignin('base', ['channel' num2str(unique_chan(i))], mean(data(:,:,:,i),3));
%    assignin('base', ['channel' num2str(unique_chan(i)) '_e'], std(data(:,:,:,i),0,3));
end

clear im_path L filenames rows columns channels fields r0 c0 I i unique_chan

%norm_sig=channel1(4,:)./(channel1(1,:).*channel2(4,:));
%norm_sig=channel2(4,:)./(channel2(1,:).*channel1(4,:));
%norm_sig=channel2(2,:)./(channel2(1,:).*channel1(2,:));