% Process Telesto raw data into a few different files:
% For each XX.raw:
%   Create a new folder XX
%   XX_interf.mha: header file with interferogram
%   XX_interf.raw: interferometric data 
%   XX_volume.mha: header file of volume
%   XX_volume.raw: image volume data

datapaths = {'../../data/FingerTip/','../../data/Mirror/','../../data/NDFilter/'};
datapaths = {'C:/Temp/'};

for m = 1:length(datapaths)
	clear volume interferogram interferogram0

    datapath = datapaths{m};

    files = dir([datapath '*.raw']);

    % Options for telesto processing
    options.dataDimension = 3;
    options.saveInterferogram = 0;
    options.saveImage = 0;
    options.showBScan = 0;

    for n =2 % 1:numel(files)
       fprintf('%s\n',files(n).name); 
       tic;
       % Create directory for raw file
       filePath = [datapath files(n).name(1:end-4) '/'];
       if ~exist(filePath,'dir')
            mkdir(filePath);
       end

       % 

       % Save interferogram as int16 
        [interferogram,info] = getInterferograms([datapath files(n).name],options);
        interferogram = int16(interferogram);
        fileName_interf = [files(n).name(1:end-4) '_interf_int16'];
        info.ElementType = 'int16';
        info.ElementDataFile = [fileName_interf '.raw'];
        info.ImagingSystem = 'Telesto';
        mha_write_header(info,[filePath fileName_interf '.mha']);
        mha_write_volume(interferogram,[filePath fileName_interf '.raw'],info.ElementType);
        fprintf('\tWrote int16 interferogram: %.2f\n',toc); tic;

       % Save interferogram at uint8 - by maximizing range...
        interferogram0 = interferogram; % Save original inteferogram (int16)
        interferogram = double(interferogram);
        interferogram = interferogram-min(interferogram(:));
        interferogram = uint8(255/max(interferogram(:))*interferogram);

        fileName_interf = [files(n).name(1:end-4) '_interf'];
        info.ElementType = 'uint8';
        info.ElementDataFile = [fileName_interf '.raw'];
        info.ImagingSystem = 'Telesto';
        mha_write_header(info,[filePath fileName_interf '.mha']);
        mha_write_volume(interferogram,[filePath fileName_interf '.raw'],info.ElementType);   

        fprintf('\tWrote uint8 interferogram: %.2f\n',toc); tic;

       % Save volume at uint8...
       %interferogram0 = interferogram0(1:10,:,:);
       volume = getBScans(double(interferogram0),options);
  
       % Convert to uint8
       volume(~isfinite(volume)) = 0;
       volume = volume-min(volume(:));
       volume = uint8(255/max(volume(:))*volume);

       fileName_vol = [files(n).name(1:end-4) '_volume'];
      % info.DimSize(3) = info.DimSize(3)/2;
       info.DimSize = size(volume);
       info.ElementType = 'uint8';
       info.ElementDataFile = [fileName_vol '.raw'];
       info.ImagingSystem = 'Telesto';
       fprintf('\tWrote uint8 volume: %.2f\n',toc); tic;

       mha_write_header(info,[filePath fileName_vol '.mha']);
       mha_write_volume(volume,[filePath fileName_vol '.raw'],info.ElementType);   

        % Move inteferogram to folder
        system(['cp ' datapath files(n).name ' ' filePath files(n).name(1:end-4) '_origInterf.raw']);
  % testTelestoFiles;
    end
    
end