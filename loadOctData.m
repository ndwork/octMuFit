function [I, z, dx, z0, zR, alpha, beta, L0, trueMu] = ...
  loadOctData(dataCase, plotIt)

  if nargin<2 plotIt=false; end;

  trueMu = [];
  cols2Del = [];
  rows2Del = [];
  alpha = [];
  beta = [];
  L0 = [];
  switch dataCase
      case 0
          N = 100;
          dx = 13d-3;
          [I,z,z0,zR,alpha,beta,L0,trueMu] = makePhantom2D(N);
          return
      case 1 %Low
          z0 = 1; %(mm)
          zR = 0.105905; %(mm)
          dx = 13d-3;
          imgDepth = 2.57;
          numPix = 512;
          datafileParts = {'..','20140428','low'};
          cols2Del = 1:52;
      case 2 %High
          z0 = 1; %(mm)
          zR = 0.105905; %(mm)
          dx = 13d-3;
          imgDepth = 2.57;
          numPix = 512;
          datafileParts = {'..','20140428','high'};
          cols2Del = 1:52;
      case 3 %Layered regular
          z0 = 1; %(mm)
          zR = 0.105905; %(mm)
          dx = 13d-3;
          imgDepth = 2.57;
          numPix = 512;
          datafileParts = {'..','20140428','layered', 'regular'};
          cols2Del = 1:52;
          rows2Del = 1:25;
      case 4 %Layered bladder
          z0 = 1; %(mm)
          zR = 0.105905; %(mm)
          dx = 13d-3;
          imgDepth = 2.57;
          numPix = 512;
          datafileParts = {'..','20140428','layered', 'bladder'};
          cols2Del = 1:52;
      otherwise
          error('Invalid data case');        
  end


  if ispc
      slash = '\';
  else
      slash = '/';
  end
  dataDir = cell2mat(strcat( datafileParts, slash ));

  z = linspace(0, imgDepth, numPix)';
  options.dataDimension = 2;
  options.saveInterferogram = 0;

  filenames=dir(dataDir);

  Isum = zeros(numPix);
  numSamps = 0;

  for i = 1:length(filenames)
      datafile = strcat(dataDir, filenames(i).name);
      [pathstr, name, ext] = fileparts(datafile);
      if(~strcmp(ext, '.raw'))  %If the file is not *.raw then skip
          continue;
      end

      [interf,info] = getInterferograms(datafile,options);
      bscans = getBScans(interf);
      I_i = 10.^(bscans./20); % Convert from db
      %I_i = bscans;
      Isum = Isum + I_i;
      numSamps = numSamps+1;
      if(plotIt && numSamps == 1)
        IsumDec = 20*log10(Isum);
        figure, plot(IsumDec(:,308),'r'), drawnow, hold on
        %figure, plot(Isum(:,256),'r'), drawnow, hold on
      end
  end

  I = Isum/numSamps;
  if(plotIt)
    Idec = 20*log10(I);
    plot(Idec(:,308),'b');
    %plot(I(:,256), 'b')
    xlabel('pixel')
    legend('10 Averaged', '100 Averaged')
    ylabel('Intensity in dB')
  end
  
  I(rows2Del, :) = [];
  z(rows2Del, :) = [];
  I(:, cols2Del) = [];
  if numel(trueMu)>0
    trueMu(rows2Del, :) = [];
    trueMu(:,cols2Del) = [];
  end

end
