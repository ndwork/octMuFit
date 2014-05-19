% LoadAndDisplayImage.m 
% processes a .raw file from the Telesto 
% - (no mha file associated with it) 
% and creates a Bscan"
% Strictly meant as a test script to show example processing

% KLL 6/2013

% Set data location, options
dataPath = '../../data/';
fileName = 'FingerTip_6_origInterf.raw'; 

options.dataDimension = 2;
options.showBscan = 0;
options.saveImage = 0;
options.nBScans = 400;
options.dataDimensions = 2;
options.saveInterferogram = 0;

saveImages = 1;
% Add file paths we'll use
    addpath ../mha_manipuluation/

% Load interferograms
    if ~exist('interf','var');
        [interf,info] = getInterferograms([dataPath fileName],options);
    end

% Load Bscans
if ~exist('bscan','var')
    bscan = getBScans(squeeze(interf),options);
end

subplot(2,1,1);
im{1} = squeeze(bscan(200,:,:))';
imagesc(im{1});
title([fileName ' x=200']);

subplot(2,1,2);
im{2} = squeeze(bscan(:,100,:))';
imagesc(im{2})
title([fileName ' y=100']);

save('ExampleImages.mat','im');
