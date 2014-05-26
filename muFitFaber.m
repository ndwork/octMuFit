
function mu = muFitFaber( I, faberPts, z, z0, zR, muFit )

  mu = zeros( numel(I), 1 );

  dy = abs( I(2:end) - I(1:end-1) );
  dy = [dy; dy(end)];
  minSig = mean( dy( dy~=0 ) );

  g = ( ( (z-z0)/zR ).^2 + 1 ).^(-1/2);

  fitOptions = optimoptions( 'lsqnonlin', ...
    'Algorithm', 'levenberg-marquardt', ...
    'MaxFunEvals',200);

  function errors = faberCost( vars )
    lmK = vars(1);
    lmMu = vars(2);
    dz = thisZ-thisZ(1);
    modelValues = lmK * exp( -2 * lmMu.*dz ) .* thisg;
    errors = ( modelValues - data ) ./ max( thisSig, minSig );
    disp( mean( errors.^2 ) );
  end


  lastFaberPt = 1;
  ks = zeros(numel(faberPts)+1,1);
  for i=1:numel(faberPts)+1
    
    if i == numel(faberPts)+1
      data = I(lastFaberPt:end);
      thisZ = z(lastFaberPt:end);
      thisSig = dy(lastFaberPt:end);
      thisg = g(lastFaberPt:end);
    else
      faberPt = faberPts(i);
      data = I(lastFaberPt:faberPt);
      thisZ = z(lastFaberPt:faberPt);
      thisSig = dy(lastFaberPt:faberPt);
      thisg = g(lastFaberPt:faberPt);
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
    else
      mu(lastFaberPt:faberPt) = lmMu;
    end

    lastFaberPt = faberPt;
  end

  %figure;
  %plot( muFit, 'r' );
  %hold on; plot(mu, 'k' );
  %axis([1 numel(I) 0 1]);

end

