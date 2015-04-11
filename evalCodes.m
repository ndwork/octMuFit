
function simIndx_out = evalCodes( simIndx, I, z, z0, zR, fid, outDir, ...
  z0_true, muString, thickString, noiseProportion, trueMu )

  mask = ones( size(I) );

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
  muFit_mVer = muFit2D_mVer( I, z, z0, zR );
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
  muFit_gBlur = muFit2D_mVer_gBlur( I, z, z0, zR );
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
  muFit_TV = muFit2D_TV( I, z, z0, zR );
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
  muFit_whTV = muFit2D_whTV( I, z, z0, zR, mask );
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

