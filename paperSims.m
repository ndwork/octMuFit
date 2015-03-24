
function paperSims
  clear; close all; rng(1);
  addpath(genpath(pwd))

  outDir = 'paperSims_out';
  mkdir(outDir);

  fid = fopen([outDir,'/sims.csv'],'w');
  fprintf(fid, ['name, algorithm, z0, z0_true, zR, mus, thicks, ', ...
    'noiseProportion, meanEtbDepth, meanVDepth, time taken (s) \n']);

  z0s = -0.5:0.25:2;
  noiseProportion = 1d-5;
  N = 100;
  simIndx = 0;

  [lambda,deltaLambda,dLambda] = getTelestoFalloffParams();


  %%  Changing Focal Plane
  disp('Changing Focal Plane');
  fprintf(fid,'Changing Focal Plane  \n');

  simMus{1} = [1,2];
  simThicks{1} = [1];

  for z0=z0s
    for i=1:numel(simMus)

      muStrings = strtrim(cellstr(num2str(simMus{i}'))');
      muString = strjoin( muStrings, '_' );

      for j=1:numel(simThicks)

        disp([ 'Working on z0=', num2str(z0), ', i=', num2str(i), ...
          ', j=', num2str(j) ]);

        [I, z, z0, zR, muAlpha, muBeta, muL0, trueMu] = makePhantom2D( ...
          N,-1,lambda, deltaLambda, dLambda,z0,noiseProportion,...
          simMus{i},simThicks{j});
        dz = z(2) - z(1);
        dx = dz;

        thickStrings = strtrim(cellstr(num2str(simThicks{j}'))');
        thickStrings{end+1} = num2str( max(z) - sum(simThicks{j}) );
        thickString = strjoin( thickStrings, '_' );

        simIndx = evalCodes( simIndx, I, z, z0, zR, lambda, ...
          deltaLambda, dLambda, fid, outDir, z0, muString, thickString, ...
          noiseProportion, trueMu );

      end
    end
  end
  clear simMus;
  clear simThicks;


  %%  Changing Mus and Thicknesses
  disp('Changing Mus and Thicknesses');
  fprintf(fid,'Changing Mus and Thicknesses  \n');
  simMus{1} = [1,4];
  simMus{2} = [1,2];
  simMus{3} = [1,1.5];

  simThicks{1} = [1];
  simThicks{2} = [0.3];
  simThicks{3} = [1.7];

  z0 = 0.5;
  
  for i=1:numel(simMus)

    muStrings = strtrim(cellstr(num2str(simMus{i}'))');
    muString = strjoin( muStrings, '_' );

    for j=1:numel(simThicks)

      disp([ 'Working on z0=', num2str(z0), ', i=', num2str(i), ...
        ', j=', num2str(j) ]);

      [I, z, z0, zR, muAlpha, muBeta, muL0, trueMu] = makePhantom2D( ...
        N,-1,lambda, deltaLambda, dLambda,z0,noiseProportion,...
        simMus{i},simThicks{j});
      dz = z(2) - z(1);
      dx = dz;

      thickStrings = strtrim(cellstr(num2str(simThicks{j}'))');
      thickStrings{end+1} = num2str( max(z) - sum(simThicks{j}) );
      thickString = strjoin( thickStrings, '_' );

      simIndx = evalCodes( simIndx, I, z, z0, zR, lambda, ...
        deltaLambda, dLambda, fid, outDir, z0, muString, thickString, ...
        noiseProportion, trueMu );

    end
  end
  clear simMus;
  clear simThicks;


  %%  Focal Length Error
  disp('Focal length error');
  fprintf( fid, 'Focal length error  \n');

  z0_true = 0.5;

  simMus{1} = [1,2];
  simThicks{1} = [1];

  zRs = [0.027, 0.1059, 0.2384] * 2 * 1.37;   % n = 1.37

  for i=1:numel(zRs)
    disp(['Working on i=', num2str(i), ' of ', num2str(numel(zRs)) ]);
    for j=1:numel(simMus)

      muStrings = strtrim(cellstr(num2str(simMus{i}'))');
      muString = strjoin( muStrings, '_' );

      [I, z, z0, zR, muAlpha, muBeta, muL0, trueMu] = makePhantom2D( ...
            N,-1,lambda, deltaLambda, dLambda,z0_true,noiseProportion,...
            simMus{i},[1],zRs(j));
      dz = z(2) - z(1);
      dx = dz;

      thickStrings = strtrim(cellstr(num2str(simThicks{1}'))');
      thickStrings{end+1} = num2str( max(z) - sum(simThicks{1}) );
      thickString = strjoin( thickStrings, '_' );

      for z0=z0_true-0.1:0.01:z0_true+0.1;
        simIndx = evalCodes( simIndx, I, z, z0, zR, lambda, deltaLambda, ...
          dLambda, fid, outDir, z0_true, muString, thickString, ...
          noiseProportion, trueMu );
      end
    end
  end
  clear simMus;
  clear simThicks;


  %%  Changing the Amount of Noise
  disp('Changing amount of noise');
  fprintf(fid,'Changing amount of noise  \n');

  simMus{1} = [1,2];
  thickString = '1_1';
  z0 = 0.5;
  noises = [ 0.00006 0.00003 0.00002 0.00001 0.000006 0.000003 ...
    0.000002 0.000001 0.0000006 0.0000003 0.0000002 0.0000001 ...
    0.00000006 0.00000003 ];

  for i=1:numel(simMus)

    muStrings = strtrim(cellstr(num2str(simMus{i}'))');
    muString = strjoin( muStrings, '_' );

    for noiseProportion=noises
      disp(['Working on i=', num2str(i), ', noise=', ...
        num2str(noiseProportion) ]);

      [I, z, z0, zR, muAlpha, muBeta, muL0, trueMu] = makePhantom2D( ...
        N, -1, lambda, deltaLambda, dLambda, z0, noiseProportion, ...
        simMus{i}, [1] );
      dz = z(2) - z(1);
      dx = dz;

      simIndx = evalCodes( simIndx, I, z, z0, zR, lambda, deltaLambda, ...
        dLambda, fid, outDir, z0, muString, thickString, ...
        noiseProportion, trueMu );
    end
  end
  clear simMus;
  clear simThicks;


  fclose(fid);
end


%% Support functions
function simIndx_out = evalCodes( simIndx, I, z, z0, zR, lambda, ...
  deltaLambda, dLambda, fid, outDir, z0_true, muString, thickString, ...
  noiseProportion, trueMu )

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

  tic;
  muFit_Ver = muFit2D_ver( I, z );
  timeTaken = toc;
  [meanEtbDepth, meanVDepth] = findErrorMetrics(muFit_Ver,trueMu,z);
  filename = ['sim_', num2str(simIndx,'%5.5i'),'.mat'];
  save( [outDir,'/',filename], 'muFit_Ver' );
  fprintf( fid, [filename, ', Vermeer, %8.3f, %8.3f, %8.3f, ', muString, ...
        ', ', thickString, ', %5.2e, %8.3f, %8.3f, %16.11f\n'], ...
        z0, z0_true, zR, noiseProportion, meanEtbDepth, meanVDepth, ...
        timeTaken );
  simIndx = simIndx + 1;

  tic;
  muFit_mVer = muFit2D_mVer( I, z, z0, zR, lambda, deltaLambda, dLambda );
  timeTaken = toc;
  [meanEtbDepth, meanVDepth] = findErrorMetrics(muFit_mVer,trueMu,z);
  filename = ['sim_', num2str(simIndx,'%5.5i'),'.mat'];
  save( [outDir,'/',filename], 'muFit_mVer' );
  fprintf( fid, [filename, ', MVM, %8.3f, %8.3f, %8.3f, ', muString, ...
        ', ', thickString, ', %5.2e, %8.3f, %8.3f, %16.11f\n'], ...
        z0, z0_true, zR, noiseProportion, meanEtbDepth, meanVDepth, ...
        timeTaken );
  simIndx = simIndx + 1;

  tic;
  muFit_gBlur = muFit2D_mVer_gBlur( I, z, z0, zR, lambda, ...
    deltaLambda, dLambda );
  timeTaken = toc;
  [meanEtbDepth, meanVDepth] = findErrorMetrics(muFit_gBlur,trueMu,z);
  filename = ['sim_', num2str(simIndx,'%5.5i'),'.mat'];
  save( [outDir,'/',filename], 'muFit_gBlur' );
  fprintf( fid, [filename, ', gBlur, %8.3f, %8.3f, %8.3f, ', muString, ...
        ', ', thickString, ', %5.2e, %8.3f, %8.3f, %16.11f\n'], ...
        z0, z0_true, zR, noiseProportion, meanEtbDepth, meanVDepth, ...
        timeTaken );
  simIndx = simIndx + 1;

  tic;
  muFit_TV = muFit2D_TV( I, z, z0, zR, lambda, deltaLambda, dLambda );
  timeTaken = toc;
  [meanEtbDepth, meanVDepth] = findErrorMetrics(muFit_TV,trueMu,z);
  filename = ['sim_', num2str(simIndx,'%5.5i'),'.mat'];
  save( [outDir,'/',filename], 'muFit_TV' );
  fprintf( fid, [filename, ', TV, %8.3f, %8.3f, %8.3f, ', muString, ...
        ', ', thickString, ', %5.2e, %8.3f, %8.3f, %16.11f\n'], ...
        z0, z0_true, zR, noiseProportion, meanEtbDepth, meanVDepth, ...
        timeTaken );
  simIndx = simIndx + 1;

  tic;
  muFit_whTV = muFit2D_whTV( I, z, z0, zR, mask, lambda, deltaLambda, ...
    dLambda );
  timeTaken = toc;
  [meanEtbDepth, meanVDepth] = findErrorMetrics(muFit_whTV,trueMu,z);
  filename = ['sim_', num2str(simIndx,'%5.5i'),'.mat'];
  save( [outDir,'/',filename], 'muFit_whTV' );
  fprintf( fid, [filename, ', whTV, %8.3f, %8.3f, %8.3f, ', muString, ...
        ', ', thickString, ', %5.2e, %8.3f, %8.3f, %16.11f\n'], ...
        z0, z0_true, zR, noiseProportion, meanEtbDepth, meanVDepth, ...
        timeTaken );
  simIndx = simIndx + 1;
  
  simIndx_out = simIndx;
end


