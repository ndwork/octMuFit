%datafile = '~/Desktop/Bladder_1_8_0.raw';
% datafile = '/Users/Genna/Documents/AttenuationProject/data/2D_10avg_FP14x14_6.raw';
datafile = '/Users/Genna/Documents/AttenuationProject/data/ACslices14APRtest_6.raw';
options.dataDimension = 2;
options.saveInterferogram = 0;

[interf,info] = getInterferograms(datafile,options);

bscans = getBScans(interf);

save( 'octDataLayers.mat', 'bscans' );


%%

% for n = 1:512
%     maxVal = 50;
%     minVal = 0;
% 
%     bscan = squeeze(bscans(n,:,:))';
%     bscanAdj = (bscan-20)/(maxVal-minVal);
%     bscanAdj(bscanAdj>1) = 1;
%     bscanAdj(bscanAdj<0) = 0;
%     imagesc(bscan)
%     title(num2str(n));
%     pause(0.5);
% end