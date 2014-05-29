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
          averageFiles = 1;
      case 2 %High
          z0 = 1; %(mm)
          zR = 0.105905; %(mm)
          dx = 13d-3;
          imgDepth = 2.57;
          numPix = 512;
          datafileParts = {'..','20140428','high'};
          cols2Del = 1:52;
          averageFiles = 1;
      case 3 %Layered regular
          z0 = 1; %(mm)
          zR = 0.105905; %(mm)
          dx = 13d-3;
          imgDepth = 2.57;
          numPix = 512;
          datafileParts = {'..','20140428','layered', 'regular'};
          cols2Del = 1:52;
          rows2Del = 1:25;
          averageFiles = 1;
      case 4 %Layered bladder
          z0 = 1; %(mm)
          zR = 0.105905; %(mm)
          dx = 13d-3;
          imgDepth = 2.57;
          numPix = 512;
          datafileParts = {'..','20140428','layered', 'bladder'};
          cols2Del = 1:52;
          averageFiles = 1;
      case 5 %Bladder
          z0 = 1;
          zR = 0.1059;
          dx = 0.00502;
          imgDepth = 2.358;
          numPix = 512;
          datafileParts = {'..', '20140519', 'Bladder', '1avg'};
          cols2Del = [];
          averageFiles = 0;
      case 6 % Retina
          z0 = 0.9;
          zR = 0.1059;
          dx = 0.00502;
          imgDepth = 2.358;
          numPix = 512;
          datafileParts = {'..', '20140519', 'Retina', '1avg'};
          cols2Del = [];
          averageFiles = 0;    
      case 7 % Layered 20140519
          z0 = 1;
          zR = 0.1059;
          dx = 0.00502;
          imgDepth = 2.358;
          numPix = 512;
          datafileParts = {'..', '20140519', 'TiO2_phantoms', 'layered', '1avg'};
          cols2Del = [];
          rows2Del = 1:50;
          averageFiles = 0;
      case 8 % Sclera
          z0 = 1;
          zR = 0.1059;
          dx = 0.00502;
          imgDepth = 2.358;
          numPix = 512;
          datafileParts = {'..', '20140519', 'Sclera', '1avg'};
          cols2Del = [];
          averageFiles = 0;
      case 9 % Skin
          z0 = 1;
          zR = 0.1059;
          dx = 0.00502;
          imgDepth = 2.358;
          numPix = 512;
          datafileParts = {'..', '20140519', 'Skin', '1avg'};
          cols2Del = [];
          averageFiles = 0;
      case 10 % Colon
          z0 = 0.8;
          zR = 0.1059;
          dx = 0.00502;
          imgDepth = 2.358;
          numPix = 512;
          datafileParts = {'..', '20140519', 'Colon', '1avg'};
          cols2Del = [];
          averageFiles = 0;
      case 11 % Intralipid 1.25
          z0 = 0.75;
          zR = 0.1059;
          dx = 0.00502;
          imgDepth = 2.57;
          numPix = 512;
          datafileParts = {'..', '20140501', 'IntralipidPhantoms', '1-25', '1avg'};
          cols2Del = [];
          averageFiles = 0;
      case 12 % Intralipid 2.5
          z0 = 0.65;
          zR = 0.1059;
          dx = 0.00502;
          imgDepth = 2.57;
          numPix = 512;
          datafileParts = {'..', '20140501', 'IntralipidPhantoms', '2-5', '1avg'};
          cols2Del = [];
          averageFiles = 0;
      case 13 % Intralipid 5
          z0 = 0.65;
          zR = 0.1059;
          dx = 0.00502;
          imgDepth = 2.57;
          numPix = 512;
          datafileParts = {'..', '20140501', 'IntralipidPhantoms', '5', '1avg'};
          cols2Del = [];
          averageFiles = 0;
      case 14 % Intralipid 10
          z0 = 0.65;
          zR = 0.1059;
          dx = 0.00502;
          imgDepth = 2.57;
          numPix = 512;
          datafileParts = {'..', '20140501', 'IntralipidPhantoms', '10', '1avg'};
          cols2Del = [];
          averageFiles = 0;
      case 15 % Intralipid 15
          z0 = 0.7;
          zR = 0.1059;
          dx = 0.00502;
          imgDepth = 2.57;
          numPix = 512;
          datafileParts = {'..', '20140501', 'IntralipidPhantoms', '15', '1avg'};
          cols2Del = [];
          averageFiles = 0;
      case 16 % Intralipid 20
          z0 = 0.5;
          zR = 0.1059;
          dx = 0.00502;
          imgDepth = 2.57;
          numPix = 512;
          datafileParts = {'..', '20140501', 'IntralipidPhantoms', '20', '1avg'};
          cols2Del = [];
          averageFiles = 0;
      case 17 % Intralipid 1.25
          z0 = 1.5;
          zR = 0.1059;
          dx = 0.00502;
          imgDepth = 2.57;
          numPix = 512;
          datafileParts = {'..', '20140519', 'Intralipid', '1p25', '1avg'};
          cols2Del = [];
          averageFiles = 0;
      case 18 % Intralipid 15
          z0 = 1.1;
          zR = 0.1059;
          dx = 0.00502;
          imgDepth = 2.57;
          numPix = 512;
          datafileParts = {'..', '20140519', 'Intralipid', '15', '1avg'};
          cols2Del = [];
          averageFiles = 0;
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

  if(averageFiles)
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
  else
    for i = 1:length(filenames)  % Find the first .raw file
      datafile = strcat(dataDir, filenames(i).name);
      [pathstr, name, ext] = fileparts(datafile);
      if(~strcmp(ext, '.raw'))  %If the file is not *.raw then skip
        continue;
      else
        break;
      end
    end
    [interf,info] = getInterferograms(datafile,options);
    bscans = getBScans(interf);
    I = 10.^(bscans./20); % Convert from db    
  end



  I(rows2Del, :) = [];
  z(rows2Del, :) = [];
  if numel(rows2Del) > 0
    zRemove = z(max(rows2Del));
    z0 = z0 - zRemove;
    z = z - zRemove;
  end
  I(:, cols2Del) = [];
  if numel(trueMu)>0
    trueMu(rows2Del, :) = [];
    trueMu(:,cols2Del) = [];
  end


  if dataCase ~= 0
    zR = 2*zR;    % Multiply by 2 to convert from zR to apparent zR
                  % See Faber paper for more details
  end


end
