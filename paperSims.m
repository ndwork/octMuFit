
function paperSims

  clear; close all; rng(1);
  addpath(genpath(pwd))

  outDir = 'paperSims_out';
  mkdir('paperSims_out');
  
  fid = fopen([outDir,'/sims.csv'],'w');
  fprintf(fid, ['name, algorithm, z0, z0_true, mus, thicks, ', ...
    'meanEtbDepth, meanVDepth \n']);


  simMus{1} = [1,4];
  simMus{2} = [1,2];
  simMus{3} = [1,1.5];
  
  simThicks{1} = [1];
  simThicks{2} = [0.3];
  simThicks{3} = [1.7];
  
  z0s = 0:0.25:2;
  
  noiseProportion = 1d-5;
  
  N = 100;

  simIndx = 0;



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%  Changing Mus and Thicknesses  %%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  for z0=z0s
    for i=1:numel(simMus)

      muStrings = strtrim(cellstr(num2str(simMus{i}'))');
      muString = strjoin( muStrings, '_' );

      for j=1:numel(simThicks)

        disp([ 'Working on z0=', num2str(z0), ', i=', num2str(i), ...
          ', j=', num2str(j) ]);

        [I, z, z0, zR, muAlpha, muBeta, muL0, trueMu] = makePhantom2D( ...
          N,-1,z0,noiseProportion,simMus{i},simThicks{j});
        dz = z(2) - z(1);
        dx = dz;

        thickStrings = strtrim(cellstr(num2str(simThicks{j}'))');
        thickStrings{end+1} = num2str( max(z) - sum(simThicks{j}) );
        thickString = strjoin( thickStrings, '_' );

        simIndx = evalCodes( simIndx, I, z, z0, zR, fid, outDir, ...
          z0, muString, thickString, noiseProportion, trueMu );

      end
    end
  end



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%  Changing the Focal Length  %%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  disp(['Working on changing focal lengths']);
  for z0=z0s
    disp(['Working on z0=', num2str(z0)]);
    [I, z, ~, zR, muAlpha, muBeta, muL0, trueMu] = makePhantom2D( ...
      N,2,z0,noiseProportion);  % bladder phantom
    dz = z(2) - z(1);
    dx = dz;

    simIndx = evalCodes( simIndx, I, z, z0, zR, fid, outDir, ...
      z0, 'bladder', 'bladder', noiseProportion, trueMu );

  end



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%  Focal Length Error  %%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  disp(['Working on focal length error']);
  
  z0_true = 0.5;

  simMus{1} = [1,4];
  simMus{2} = [1,2];
  simThicks{1} = [1];

  for i=1:numel(simMus)
    disp(['Working on i=', num2str(i)]);

    muStrings = strtrim(cellstr(num2str(simMus{i}'))');
    muString = strjoin( muStrings, '_' );

    [I, z, z0, zR, muAlpha, muBeta, muL0, trueMu] = makePhantom2D( ...
          N,-1,z0_true,noiseProportion,[1 2],[1]);
    dz = z(2) - z(1);
    dx = dz;

    thickStrings = strtrim(cellstr(num2str(simThicks{1}'))');
    thickStrings{end+1} = num2str( max(z) - sum(simThicks{1}) );
    thickString = strjoin( thickStrings, '_' );

    for z0=z0_true-0.4:0.1:z0_true+0.4;
      simIndx = evalCodes( simIndx, I, z, z0, zR, fid, outDir, ...
        z0_true, muString, thickString, noiseProportion, trueMu );
    end
  end



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%  Changing the Amount of Noise  %%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  disp(['Working on changing amount of noise']);
  
  simMus{1} = [1,4];
  simMus{2} = [1,2];
  thickString = '1_1';
  dNoise = 1d-5 / 3;
  z0 = 0.5;
  for i=1:numel(simMus)
    
    muStrings = strtrim(cellstr(num2str(simMus{i}'))');
    muString = strjoin( muStrings, '_' );
    
    for noiseProportion=0:dNoise:2d-5
      disp(['Working on i=', num2str(i), ', noise=', num2str(noiseProportion) ]);

      [I, z, z0, zR, muAlpha, muBeta, muL0, trueMu] = makePhantom2D( ...
        N, -1, z0, noiseProportion, [1 2], [1] );
      dz = z(2) - z(1);
      dx = dz;

      simIndx = evalCodes( simIndx, I, z, z0, zR, fid, outDir, ...
        z0, muString, thickString, noiseProportion, trueMu );
    end
  end


  fclose(fid);
end



function simIndx_out = evalCodes( simIndx, I, z, z0, zR, fid, outDir, ...
  z0_true, muString, thickString, noiseProportion, trueMu )

  mask = ones( size(I) );

  %TV Parameters
  maxIter = 1000;
  tau = .02;  % Actually, the CP codes choose step sizes tau = sigma = 1/nrmK.
  overRelax = 1.9;
  epsilon = 1d-1;
  weight = 1 ./ ( I + epsilon );
  gamma = .5*weight;
  showTrigger = pi;  %by making irrational, never shows
  paramsCP = struct('gamma',gamma,'tau',tau,'maxIter',maxIter,...
            'showTrigger',showTrigger,'theta',1,'overRelax',overRelax);

  muFit_mVer = muFit2D_mVer( I, z, z0, zR );
  [meanEtbDepth, meanVDepth] = findErrorMetrics(muFit_mVer,trueMu,z);
  filename = ['sim_', num2str(simIndx,'%5.5i'),'.mat'];
  save( [outDir,'/',filename], 'muFit_mVer' );
  fprintf( fid, [filename, ', MVM, %8.3f, %8.3f, ', muString, ...
        ', ', thickString, ', %5.2e\n, %8.3f, %8.3f'], ...
        z0, z0_true, noiseProportion, meanEtbDepth, meanVDepth );
  simIndx = simIndx + 1;

  muFit_gBlur = muFit2D_mVer_gBlur( I, z, z0, zR );
  [meanEtbDepth, meanVDepth] = findErrorMetrics(muFit_gBlur,trueMu,z);
  filename = ['sim_', num2str(simIndx,'%5.5i'),'.mat'];
  save( [outDir,'/',filename], 'muFit_gBlur' );
  fprintf( fid, [filename, ', gBlur, %8.3f, %8.3f, ', muString, ...
        ', ', thickString, ', %5.2e\n'], z0, z0_true, noiseProportion );
  simIndx = simIndx + 1;

  muFit_TV = muFit2D_TV( I, z, z0, zR );
  [meanEtbDepth, meanVDepth] = findErrorMetrics(muFit_TV,trueMu,z);
  filename = ['sim_', num2str(simIndx,'%5.5i'),'.mat'];
  save( [outDir,'/',filename], 'muFit_TV' );
  fprintf( fid, [filename, ', TV, %8.3f, %8.3f, ', muString, ...
        ', ', thickString, ', %5.2e\n'], z0, z0_true, noiseProportion );
  simIndx = simIndx + 1;

  muFit_whTV = muFit2D_whTV( I, z, z0, zR, mask );
  [meanEtbDepth, meanVDepth] = findErrorMetrics(muFit_whTV,trueMu,z);
  filename = ['sim_', num2str(simIndx,'%5.5i'),'.mat'];
  save( [outDir,'/',filename], 'muFit_whTV' );
  fprintf( fid, [filename, ', whTV, %8.3f, %8.3f, ', muString, ...
        ', ', thickString, ', %5.2e\n'], z0, z0_true, noiseProportion );
  simIndx = simIndx + 1;


  %noiseLevel = 0;
  %[M N] = size( I );
  %%I = max( I - noiseLevel, 0 );
  %%I = I .* mask;
  %%I = I ./ 1000;  % don't divide by 1000 for simulation

  %col = round(N/2);
  %lineI = I(:,col);
  %mask = trueMu > 0;
  %lineMask = mask(:,col);
  
  
  %muFit_RICR = muFitCVX( lineI, lineMask, z, z0, zR, eta );
  %muFit_wTV = weightedTvDenoise_CP( muFit_mVer, paramsCP );
  
  
  simIndx_out = simIndx;
  
end


function [meanEtbDepth, meanVDepth] = findErrorMetrics(muFit,trueMu,z)
  offsetThreshPercent = 0.05;
  [M N] = size(muFit);
  
  etbDepths = zeros(1,N);
  vDepths = zeros(1,N);
  
  for i=1:N
    thisETB = findErrorEnergyTooBigDepth( offsetThreshPercent, ...
      trueMu, muFit, z );
    thisV = findViolationDepth( offsetThreshPercent, ...
      trueMu, muFit, z );
    
    etbDepths(i) = thisETB;
    vDepths(i) = thisV;
  end
  
  meanEtbDepth = mean(etbDepths);
  meanVDepth = mean(vDepths);
end




