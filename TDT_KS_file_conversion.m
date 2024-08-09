%% TDTtoKSfileconversion

%% combining files the easy way

BLOCKPATH = 'F:\LAKHANI\240705\WT-240705-112721';
BLOCKPATH2 = 'F:\LAKHANI\240705\WT-240705-115002';
BLOCKPATH3 = 'F:\LAKHANI\240705\WT-240705-120445';
data = TDTbin2mat(BLOCKPATH, 'TYPE', {'streams'});
data2 = TDTbin2mat(BLOCKPATH2, 'TYPE', {'streams'});
data3 = TDTbin2mat(BLOCKPATH3, 'TYPE', {'streams'});
wav = data.streams.Raws.data;
wav2 = data2.streams.Raws.data;
wav3 = data3.streams.Raws.data;
finalwav = horzcat(wav, wav2, wav3);
fid = fopen('D:\Data\240705\file_240705.dat', 'w');
fwrite(fid, finalwav, 'int16');
fclose(fid);

%% the loop that works for 64chan probe if I want to combine all files in original form in a folder

expday = '230220';
folder = strcat('F:\LAKHANI\', expday);
list = dir(folder);  
list(1:2) = [];

wav = struct;
%need to make sure I'm in that folder on the top left in matlab 
for i = 1:length(list) 
    BLOCKPATH = list(i).name;
    data = TDTbin2mat(BLOCKPATH, 'TYPE', {'streams'});
    wav(i).data = data.streams.Raws.data;
end
finalwav = [];
for i = 1:length(wav)
    finalwav = horzcat(finalwav, wav(i).data);
end

fid = fopen(fullfile('C:\Users\lakha\Documents\Data\', expday, 'file_230220.dat'), 'w');
fwrite(fid, finalwav, 'int16');
fclose(fid);
