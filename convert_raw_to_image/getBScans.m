function BScan = getBScans(interf)
    % function BScan = getBScans(interf,options)
    % Process interferogram datasets to get BScans
    % Input: fileAndPath = full file and path
    %       options.dataDimension = {2,3}
    %       options.showBScan = {1,0} - will plot B scans
    %       options.saveImage = {1,0} - save image
    %       nBScans = number of B-scans per C-scan
    %       BScan output
        load('NDFT_coeff.mat','NW','coeff');

% Process data to get interferogram
    if ndims(interf)==2
        BScan = get2DScans(interf);
        
    else
        BScan = get3DScans(interf,NW);
%         if options.saveInterferogram
%             for n = 1:length(nBScans)
%                 dlmwrite([fileAndPath(1:end-10) 'BScan_' num2str(n) '.txt'],BScan);
%             end
%         end
    
    end


end

function bscan = get2DScans(interf)
   % interf = dlmread(fileAndPath);
    FFT_length = 1024; % Is actually in original header file, but this should be correct
    load FFTM	%This file actually contains both the forward and inverse matrices, but usually we only need the forward (so the filesize can be reduced if necessary)
    NDFT = TempFFTM(1:FFT_length/2,:)*interf; %compute DFT for first 512 freq indices

    bscan=20*log10(abs(NDFT));

end

function bscan = get3DScans(interf,NW)
   % interf = dlmread(fileAndPath);
    FFT_length = 1024;

    %compute DFT for first 512 freq indices
    sz = size(interf);
    bscan = zeros(sz(1),sz(2),sz(3)/2);
    if ndims(interf)==2
        bscan = zeros(sz(2),1,sz(1)/2);
    end
    
    tic
    for n = 1:sz(1)
        tic;
            NDFT = NW(1:FFT_length/2,:)*squeeze(interf(n,:,:))';
            bscan(n,:,:)=20*log10(abs(NDFT'));

%         for m = 1:sz(2)
%             NDFT = NW(1:FFT_length/2,:)*double(squeeze(
%             NDFT = NW(1:FFT_length/2,:)*double(squeeze(interf(n,m,:)));
%             bscan(n,m,:)=20*log10(abs(NDFT'));
%         end
        
        if mod(n,10)==0
            fprintf('%d %.2f\n',n,toc);
        end
    end
end



