
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

  %[lambda,deltaLambda,dLambda] = getTelestoFalloffParams();


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

        [I, z, z0, zR, ~, ~, ~, trueMu] = makePhantom2D( ...
          N,-1,z0,noiseProportion,simMus{i},simThicks{j});
        dz = z(2) - z(1);
        dx = dz;

        thickStrings = strtrim(cellstr(num2str(simThicks{j}'))');
        thickStrings{end+1} = num2str( max(z) - sum(simThicks{j}) );
        thickString = strjoin( thickStrings, '_' );

        simIndx = evalCodes( simIndx, I, z, z0, zR, fid, outDir, z0, ...
          muString, thickString, noiseProportion, trueMu );

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

      [I, z, z0, zR, ~, ~, ~, trueMu] = makePhantom2D( ...
        N, -1, z0, noiseProportion, simMus{i}, simThicks{j} );
      dz = z(2) - z(1);
      dx = dz;

      thickStrings = strtrim(cellstr(num2str(simThicks{j}'))');
      thickStrings{end+1} = num2str( max(z) - sum(simThicks{j}) );
      thickString = strjoin( thickStrings, '_' );

      simIndx = evalCodes( simIndx, I, z, z0, zR, fid, outDir, z0, ...
        muString, thickString, noiseProportion, trueMu );

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

  zRs = [0.05, 0.1059, 0.2384] * 2 * 1.37;   % n = 1.37

  for i=1:numel(zRs)
    disp(['Working on i=', num2str(i), ' of ', num2str(numel(zRs)) ]);
    for j=1:numel(simMus)

      muStrings = strtrim(cellstr(num2str(simMus{j}'))');
      muString = strjoin( muStrings, '_' );

      [I, z, z0, zR, muAlpha, muBeta, muL0, trueMu] = makePhantom2D( ...
            N, -1, z0_true, noiseProportion, simMus{j}, [1], zRs(i) );
      dz = z(2) - z(1);
      dx = dz;

      thickStrings = strtrim(cellstr(num2str(simThicks{1}'))');
      thickStrings{end+1} = num2str( max(z) - sum(simThicks{1}) );
      thickString = strjoin( thickStrings, '_' );

      for z0=z0_true-0.1:0.01:z0_true+0.1;
        simIndx = evalCodes( simIndx, I, z, z0, zR, fid, outDir, ...
          z0_true, muString, thickString, noiseProportion, trueMu );
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

      [I, z, z0, zR, ~, ~, ~, trueMu] = makePhantom2D( ...
        N, -1, z0, noiseProportion, simMus{i}, [1] );
      dz = z(2) - z(1);
      dx = dz;

      simIndx = evalCodes( simIndx, I, z, z0, zR, fid, outDir, z0, ...
        muString, thickString, noiseProportion, trueMu );
    end
  end
  clear simMus;
  clear simThicks;


  fclose(fid);
end

