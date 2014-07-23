
function paperSims

  clear; close all; rng(1);
  addpath(genpath(pwd))

  outDir = 'paperSims_out';
  mkdir('paperSims_out');
  
  fid = fopen([outDir,'/sims.csv'],'w');
  fprintf(fid, 'name, algorithm, z0, z0_true, mus, thicks \n');


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

        thickStrings = strtrim(cellstr(num2str(simThicks{i}'))');
        thickStrings{end+1} = num2str( max(z) - sum(simThicks{i}) );
        thickString = strjoin( thickStrings, '_' );

        [I, z, z0, zR, muAlpha, muBeta, muL0, trueMu] = makePhantom2D( ...
          N,-1,z0,noiseProportion,simMus{i},simThicks{j});
        dz = z(2) - z(1);
        dx = dz;

        simIndx = evalCodes( simIndx, I, z, z0, zR, fid, outDir, ...
          z0_true, muString, thickString );

      end
    end
  end



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%  Changing the Focal Length  %%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  for z0=z0s
    [I, z, z0, zR, muAlpha, muBeta, muL0, trueMu] = makePhantom2D( ...
      N,2,z0,noiseProportion);  % bladder phantom
    dz = z(2) - z(1);
    dx = dz;

    simIndx = evalCodes( simIndx, I, z, z0, zR, fid, outDir, ...
      z0_true, 'bladder', 'bladder' );

  end



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%  Focal Length Error  %%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%

  z0_true = 0.6;

  simMus{1} = [1,4];
  simMus{2} = [1,2];
  simThicks{1} = [1];

  for i=1:numel(simMus)

    muStrings = strtrim(cellstr(num2str(simMus{i}'))');
    muString = strjoin( muStrings, '_' );

    thickStrings = strtrim(cellstr(num2str(simThicks{1}'))');
    thickStrings{end+1} = num2str( max(z) - sum(simThicks{1}) );
    thickString = strjoin( thickStrings, '_' );

    [I, z, z0, zR, muAlpha, muBeta, muL0, trueMu] = makePhantom2D( ...
          N,-1,z0_true,noiseProportion,[1 2],[1]);
    dz = z(2) - z(1);
    dx = dz;
    
    for z0=z0_true-0.4:0.1:z0_true+0.4;
      simIndx = evalCodes( simIndx, I, z, z0, zR, fid, outDir, ...
        z0_true, muString, thickString );
    end
  end



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%  Changing the Amount of Noise  %%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  simMus{1} = [1,4];
  simMus{2} = [1,2];
  thickString = '1_1';
  dNoise = 1d-5 / 3;
  for i=1:numel(simMus)
    
    muStrings = strtrim(cellstr(num2str(simMus{i}'))');
    muString = strjoin( muStrings, '_' );
    
    for noiseProportion=0:dNoise:2d-5
      [I, z, z0, zR, muAlpha, muBeta, muL0, trueMu] = makePhantom2D( ...
        N, -1, z0_true, noiseProportion, [1 2], [1] );
      dz = z(2) - z(1);
      dx = dz;

      simIndx = evalCodes( simIndx, I, z, z0, zR, fid, outDir, ...
        z0_true, muString, thickString );
    end
  end
  

  fclose(fid);
end



function simIndx_out = evalCodes( simIndx, I, z, z0, zR, fid, outDir, ...
  z0_true, muString, thickString )


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
  filename = ['sim_', num2str(simIndx,'%5.5i'),'.mat'];
  save( [outDir,'/',filename], 'muFit_mVer' );
  fprintf( fid, [filename, ', gBlur, %8.3f, %8.3f, ', muString, ...
        ', ', thickString, '\n'], z0, z0_true );
  simIndx = simIndx + 1;

  muFit_gBlur = muFit2D_mVer_gBlur( I, z, z0, zR );
  filename = ['sim_', num2str(simIndx,'%5.5i'),'.mat'];
  save( [outDir,'/',filename], 'muFit_gBlur' );
  fprintf( fid, [filename, ', gBlur, %8.3f, %8.3f, ', muString, ...
        ', ', thickString, '\n'], z0, z0_true );
  simIndx = simIndx + 1;

  muFit_TV = muFit2D_TV( I, z, z0, zR );
  filename = ['sim_', num2str(simIndx,'%5.5i'),'.mat'];
  save( [outDir,'/',filename], 'muFit_TV' );
  fprintf( fid, [filename, ', gBlur, %8.3f, %8.3f, ', muString, ...
        ', ', thickString, '\n'], z0, z0_true );
  simIndx = simIndx + 1;

  muFit_whTV = muFit2D_whTV( I, z, z0, zR, mask );
  filename = ['sim_', num2str(simIndx,'%5.5i'),'.mat'];
  save( [outDir,'/',filename], 'muFit_whTV' );
  fprintf( fid, [filename, ', gBlur, %8.3f, %8.3f, ', muString, ...
        ', ', thickString, '\n'], z0, z0_true );
  simIndx = simIndx + 1;

  
  
  %mask = ones( size(I) );
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

