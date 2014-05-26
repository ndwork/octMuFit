function [interf,info] = getInterferograms(fileAndPath,options)
    % Takes raw data file from Telesto and extracts interferogram doing
    % apodization correction
    % Input: fileAndPath = full raw file and path 
    %        options.dataDimension = {2,3}
    %        options.saveInterferogram = {1,0}
    % Output: interf = apodization corrected interferograms in a matrix
    % NOTE: The raw data is in int16 fmt; apodization correction is in
    % short int
    
% Process data to get interferogram and write to file
    if options.dataDimension == 2
        [interf,info] = get2DInterferograms(fileAndPath);
        if options.saveInterferogram
            dlmwrite([fileAndPath(1:end-4) '_interf.txt'],interf);
        end
    else
        [interf,info] = get3DInterferograms(fileAndPath);
                
        if options.saveInterferogram
            dlmwrite([fileAndPath(1:end-4) '_interf.txt'],interf(:,:,1));
            for n = 2:size(interf,3)
               dlmwrite([fileAndPath(1:end-4) '_interf.txt'],interf(:,:,n),'-append');
            end
        end
    end


end

function [interfL,info] = get2DInterferograms(fileAndPath)
    % Get interferograms from a 2D dataset
    fclose('all');
    error =0;
    fp=fopen(fileAndPath, 'r', 'l'); %'b' for big endian
    fseek(fp,7,'bof');
    dim=fread(fp,1,'*char');
    fseek(fp,12,'cof');
    lat_x=fread(fp,1,'uint32');
    lat_z=fread(fp,1,'uint32');
    lat_y=fread(fp,1,'uint32');
    num_frame=fread(fp,1,'uint32');
    FFT_length=fread(fp,1,'uint32');
    fseek(fp,24,'cof');
    lambda_c=fread(fp,1,'float32');
    bandwidth=fread(fp,1,'float32');
    scan_width=fread(fp,1,'float32');
    scan_length=fread(fp,1,'float32');
    fseek(fp,4,'cof');
    apod_cycles=fread(fp,1,'short');
    ascan_ave=fread(fp,1,'short'); 
    bytes_per_data=fread(fp,1,'uint8'); 
    fseek(fp,463,'cof'); 
    frame0=fread(fp,apod_cycles*FFT_length*ascan_ave,'int16'); %ok
    frame=fread(fp,lat_x*FFT_length*ascan_ave,'int16'); %ok
    apod = reshape(frame0,FFT_length,apod_cycles);
    apod = flipud(apod);
    apod = mean((apod),2);
    interf=reshape(frame,FFT_length,lat_x);
    offset_error_vect=fread(fp,FFT_length,'single=>double');
    chirp_vect=fread(fp,FFT_length,'single=>double');
    fseek(fp,104,'bof');
    spacing_z=fread(fp,1,'float32');
    fclose(fp);
    
	info.DimSize = [lat_x lat_y 2*lat_z];
    info.ElementSpacing = [scan_length scan_width 2.54*512/lat_z];
    info.NDims = 2;

    
    % Use apodization data
    interf = flipud(interf);
    interf = interf - repmat(apod,1,lat_x);

    %The Bscan sometimes looks more accurate if you multiply the interferogram with a Gaussian. Currently disabled.
    if (0)
        gaussy = 0.1*sqrt(2*log(2));	%FWHM width of Gaussian
        interfL = interf.*repmat(400*gausswin(1024,gaussy)/sum(gausswin(1024,gaussy)),1,size(interf,2));
    else
        interfL = interf;
    end

    %Enable this if you want to modify the interferogram and write it to a file, so that it can be opened with the Telesto.
    if (0)
        interfL = interfL + repmat(apod,1,lat_x);
        interfL = flipud(interfL);
        frame=reshape(interfL,FFT_length*lat_x,1);

        fclose('all');
        error =0;
        file_name='Output Name.raw';	%This file must already exist.  Copy the input file and rename it if necessary. Interferogram in file will be overwritten.
        fp=fopen(file_name, 'r+', 'l'); %'b' for big endian
        fseek(fp,7,'bof');
        dim=fread(fp,1,'*char');
        fseek(fp,12,'cof');
        lat_x=fread(fp,1,'uint32');
        lat_z=fread(fp,1,'uint32');
        lat_y=fread(fp,1,'uint32');
        num_frame=fread(fp,1,'uint32');
        FFT_length=fread(fp,1,'uint32');
        fseek(fp,24,'cof');
        lambda_c=fread(fp,1,'float32');
        bandwidth=fread(fp,1,'float32');
        scan_width=fread(fp,1,'float32');
        scan_length=fread(fp,1,'float32');
        fseek(fp,4,'cof');
        apod_cycles=fread(fp,1,'short');
        ascan_ave=fread(fp,1,'short'); 
        bytes_per_data=fread(fp,1,'uint8'); 
        fseek(fp,463,'cof'); 
        fseek(fp,25600*2,'cof');
        fwrite(fp,frame,'int16'); %ok
        fclose(fp);
    end
end

function [interfL,info] = get3DInterferograms(fileAndPath)
    % Get interferograms from a 3D dataset
    
    fclose('all');
    error =0;

    %fileID = fopen(filename, permission, machineformat)
    fp=fopen(fileAndPath, 'r', 'l'); %'b' for big endian

    fseek(fp,7,'bof');
    dim=fread(fp,1,'*char');
    fseek(fp,12,'cof');

    lat_x=fread(fp,1,'uint32');
    lat_z=fread(fp,1,'uint32');
    lat_y=fread(fp,1,'uint32');

    num_frame=fread(fp,1,'uint32');
    FFT_length=fread(fp,1,'uint32');

    fseek(fp,24,'cof');
    lambda_c=fread(fp,1,'float32');
    bandwidth=fread(fp,1,'float32');
    scan_width=fread(fp,1,'float32');
    scan_length=fread(fp,1,'float32');
    fseek(fp,4,'cof');
    if(ftell(fp) ~= 84)
        error = 1;
    end
    apod_cycles=fread(fp,1,'short');
    ascan_ave=fread(fp,1,'short');
    bytes_per_data=fread(fp,1,'uint8');

    % Populate header structure
    info.DimSize = [lat_x lat_y 2*lat_z];
    info.ElementSpacing = [scan_length scan_width 2.57*512/lat_z];
    info.NDims = 3;
    fseek(fp,512,'bof'); %this lands you on frame header

    ph_mat = zeros(lat_y,lat_x);
    prin_ph_mat = zeros(lat_y,lat_x);
    interfL = zeros(lat_y,lat_x,2*lat_z);
    for n1=1:lat_y    
        junk = fread(fp,40,'uint8');
        
        % Get apodization for this frame, reshape, and average 25 apods
        apod=fread(fp,apod_cycles*FFT_length*ascan_ave,'int16'); %ok
        try
            apod = reshape(apod,FFT_length,apod_cycles);
        catch
            disp('Problem reshaping apodization');
        end
        apod = flipud(apod);
        apod = mean((apod),2);
        apod = round(apod); % keep in integer format
        frame=fread(fp,lat_x*FFT_length*ascan_ave,'int16'); %ok

        try
            interf=reshape(frame,FFT_length,lat_x);
        catch
            disp('Problem reshaping interfereogram');
        end

        interf = flipud(interf);   

        if max(max(interf)) == 4095 %if saturation
            is_sat = 1;
            disp('Interferogram is saturated');
        end
        
        % Do apodization correction
        interf = interf - repmat(apod,1,lat_x);

        interfL(n1,:,:) = interf';
        
    end

end
