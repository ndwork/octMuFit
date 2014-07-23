
close all; clear;


z = 0:0.01:2;
z0 = 0.6;
zR = 0.2;  


h = ( ( (z-z0)./zR ).^2 + 1 ).^(-1/2);


plot( fliplr(h), z, 'k', 'LineWidth', 2 );
%xlabel('Confocal Intensity', 'FontSize', 16, 'FontWeight', 'bold');
%ylabel('Depth', 'FontSize', 16, 'FontWeight', 'bold');

ah = gca;   % get the current axes handle
set(ah,'xtick',[]);
set(ah,'xticklabel',[]);
set(ah,'ytick',[]);
set(ah,'yticklabel',[]);
