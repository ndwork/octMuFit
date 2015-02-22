
function mu = muFitFaber( I, faberPts, z, z0, zR )

  mu = zeros( numel(I), 1 );
  k = zeros( numel(I), 1 );

  dy = abs( I(2:end) - I(1:end-1) );
  dy = [dy; dy(end)];
  minSig = mean( dy( dy~=0 ) );

  h = makeConfocalFunction( z, z0, zR );

  fitOptions = optimoptions( 'lsqnonlin', ...
    'Algorithm', 'levenberg-marquardt', ...
    'display', 'off', ...
    'TolFun', 1d-14, ...
    'TolX', 1d-14, ...
    'MaxFunEvals', 10000, ...
    'MaxIter', 10000 ...
  );

  function errors = faberCost( vars )
    lmK = vars(1);
    lmMu = vars(2);
    dz = thisZ-thisZ(1);
    modelValues = lmK * exp( -2*lmMu .* dz ) .* thisH;
    %errors = ( modelValues - data ) ./ max( thisSig, minSig );
    errors = modelValues - data;
    %disp( mean( errors.^2 ) );
  end


  lastFaberPt = 1;
  ks = zeros(numel(faberPts)+1,1);
  for i=1:numel(faberPts)+1
    
    if i == numel(faberPts)+1
      data = I(lastFaberPt:end);
      thisZ = z(lastFaberPt:end);
      thisSig = dy(lastFaberPt:end);
      thisH = h(lastFaberPt:end);
    else
      faberPt = faberPts(i);
      data = I(lastFaberPt:faberPt);
      thisZ = z(lastFaberPt:faberPt);
      thisSig = dy(lastFaberPt:faberPt);
      thisH = h(lastFaberPt:faberPt);
    end

    lmOut = lsqnonlin( @faberCost, [1 0], [], [], fitOptions );
    lmK = lmOut(1);
    lmMu = lmOut(2);

    %close all;
    %dz = thisZ - thisZ(1);
    %resultValues = lmK * exp( -lmMu.*dz ) .* thisg;
    %plot( resultValues, 'r' );
    %hold on; plot( data, 'b' );

    if lmK < 1d-3
      lmMu = 0;
    end

    if i == numel(faberPts)+1
      mu(lastFaberPt:end) = lmMu;
      k(lastFaberPt:end) = lmK;
    else
      mu(lastFaberPt:faberPt) = lmMu;
      k(lastFaberPt:faberPt) = lmK;
      lastFaberPt = faberPt;
    end
  end

  
  %dz = z - z(1);
  %modelValues = k .* exp( -2 .* mu .* dz ) .* g;
  %plot( modelValues, 'r-' );
  %hold on; plot( I, 'k:' );
  
  
  %figure;
  %plot( muFit, 'r' );
  %hold on; plot(mu, 'k' );
  %axis([1 numel(I) 0 1]);

end

