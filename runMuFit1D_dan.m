
function runMuFit1D
  clear; close all; rng(1);
  addpath(genpath(pwd))

  muAlpha = 0;    % Note, if we know true value than problem is better

  dataCase = 0;
  [I, z, dx, z0, zR, muAlpha, muBeta, muL0, ALA, trueMu ] = loadOctData( dataCase, false );

  if ~isvector(I)
    if dataCase==0
      trueMu = trueMu(:,1);
      I = I(:,1);
    else
      I = I(:,250:274);
      I = mean( I, 2 );
    end
  end  

  if dataCase == 0
    mask = trueMu > 0;
mask(end-120:end) = 0;
%I = I .* mask;
    dMu = trueMu(2:end) - trueMu(1:end-1);
    faberPts = find( abs(dMu) > 0 );
    eta = 1d-3;
  else
    if dataCase == 7
      faberPts = [ 150, 189, 229, 270, 308 ];
    elseif dataCase == 17
      faberPts = [ 210, 350 ];
    end

    [mask, noiseLevel] = makeTheMask( I, ALA );
    I = I - noiseLevel;
    I = max( I, 0 );
    I(mask==0) = 0;
    I = I ./ 1000;
    %eta = 1;
    %eta = 1d-1;  % Best answer so far
  end



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  %%%%%%%%%% My code begins here. %%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % check that I seems correct.
  alpha = 0.2;
  beta = 1.0;
  L0 = 2000;
  dz = z(2) - z(1);
  muIntegral = cumsum(trueMu)*dz;
  h = (((z - z0)/zR).^2 + 1).^(-.5);
  Icheck = dz*alpha*beta*L0*(trueMu.*h).*exp(-2*muIntegral);  
  
%   I = dz*I; % QUESTION 1: Doesn't this line need to be here? (See eq. 7 in paper.)  
    
  fnd = find(mask ~= 0);
  Md = fnd(end);  

  A = zeros(Md - 1, Md);
  b = zeros(Md - 1,1);

  for m = 1:(Md - 1)
      A(m,m) = I(m)*(1 + log(h(m+1)/h(m)));
      A(m,Md) = -I(Md);
      for n = (m+1):(Md - 1)
          A(m,n) = I(n)*log(h(n+1)/h(n));
      end
      b(m) = sum(I(m:Md - 1));
  end
  A = (.5/dz)*A;

  trueGamma = 1./trueMu;
  res = A*trueGamma(1:Md) - b; % QUESTION 2: Why is this not small ?  
                               %It gets much smaller with the line I = dz*I, but
                               %not small enough.
  plot( res );  title('A*trueGamma - b');


  D = zeros(Md,Md);
  for m = 1:(Md-1)
      D(m,m) = -1;
      D(m,m+1) = 1;
  end

  cvx_begin

    cvx_precision best
    variable x(Md)
    minimize .5*sum((A*x - b).^2) + (eta/dz)*sum(abs(D*x))    
    subject to 
        x >= 0
    
  cvx_end % NOTE that x agrees closely with gammaStar below.

  muCheck = 1./x;
  

%   res = res./b; % Should I look at this relative residual instead?  Still doesn't seem that small.                           

  % IDEA: When deep values of I are nearly zero, modifying deep values of gamma
  % doesn't help much with the data term, despite the fact that deep values of gamma
  % appear in almost all subterms in the data term.  However, if deep values of I are
  % not nearly zero, then modifying deep values of gamma may help a lot
  % with the data term.  Perhaps that sheds some light on why deep values
  % of gamma are being corrupted when we do mask(end-120:end) = 0;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%% My code ends here. %%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  [muStar,gammaStar] = muFitCVX_dan( I, mask, z, z0, zR, eta );  
  res = A*gammaStar - b;
  figure, plot( A*gammaStar - b, 'b' );
  title('A*gammaStar - b');

%load muStar.mat
  tic
  muFit = muStar;
  %[muFit, fos, relFos] = muFit1D_CP( I, mask, z, z0, zR, eta, muStar );
  %[muFit, fos, relFos] = muFit1D_CP2( I, mask, z, z0, zR, eta, muStar );
  %[muFit, fos, relFos] = muFit1D_LADMM(I, mask, z, z0, zR, eta, muStar );
  %[muFit, fos, relFos] = muFit1D_ADMM(I, mask, z, z0, zR, eta, muStar );
  timeTaken = toc;
  disp(['Time taken: ', num2str(timeTaken)]);
  if exist('fos')
    figure, semilogy(fos);  title('fos');
  end
  if exist('relFos')
    figure, semilogy(relFos, 'LineWidth', 2);
    xlabel('Iteration', 'FontSize', 14)
    ylabel('Relative Error', 'FontWeight', 'bold', 'FontSize', 14)
  end


  figure;
  plot( z, muFit, 'b', 'LineWidth', 3 );
  maxMu2Display = 6.0;
  a = [0 max(z) 0 maxMu2Display];
  axis(a)
  hold on;
  if numel( trueMu ) > 0
    plot( z, trueMu, 'r--', 'LineWidth', 3 );
    muFaber = muFitFaber( I, faberPts, z, z0, zR );
    plot( z, muFaber, 'k:', 'LineWidth', 3 );
    h_legend = legend( 'RICR', 'True \mu', 'Faber', 'Location', 'NorthWest' );
  elseif dataCase==7 || dataCase==17
    muFaber = muFitFaber( I, faberPts, z, z0, zR );
    plot( z, muFaber, 'g--', 'LineWidth', 3 );
    h_legend = legend( 'RICR', 'Faber Fit', 'FontSize', 12, 'Location', 'NorthWest' );
  end
  axis(a)
  set(h_legend,'FontSize',14);
  xlabel('Depth (mm)','FontSize',16, 'FontWeight', 'bold');
  ylabel('Attenuation (mm^{-1})', 'FontSize', 14, 'FontWeight', 'bold');

  gamFit = 1 ./ muFit;
  figure;
  plot( z, gamFit, 'r', 'LineWidth', 2 );
  hold on;
  if numel( trueMu ) > 0
    trueGam = 1 ./ trueMu;
    plot( z, trueGam, 'b', 'LineWidth', 2 );
    legend( 'CVX Fit', 'True \gamma' );
  end
  ylim([0 2.5]);
  xlabel('Depth (mm)');
  ylabel('Gamma (mm)');

  figure;
  plot( z, I );
  fitI = dz*mu2I( muFit, z, z0, zR, muAlpha, muBeta, muL0 );
  scale = median(I(mask~=0)) / median( fitI(mask~=0) );
  if scale==0
    scale = mean(I(mask~=0)) / mean( fitI(mask~=0) );
  end
  fitI = fitI * scale;
  hold on;
  plot( z, fitI, 'r--', 'LineWidth', 3 );
  if numel( trueMu )
    trueI = dz*mu2I( trueMu, z, z0, zR, muAlpha, muBeta, muL0 );
    plot( z, trueI, 'k', 'LineWidth', 2 );
    h_legend = legend('data','RICR fit','actual');
    set(h_legend,'FontSize',14);
  end
  xlabel('Depth (mm)', 'FontSize', 14, 'FontWeight', 'bold');
  ylabel('Intensity', 'FontSize', 14, 'FontWeight', 'bold');
  
  zless1mm = z < z(1)+0.75; %Find depths for which we are less than 0.75mm away from top surface 
  diffI = norm(I(zless1mm) - fitI(zless1mm));
  disp(['Norm difference between I and fit of I: ' num2str(diffI)])


  dz = z(2)-z(1);
  gam = 1 ./ muFit;
  gam(mask==0) = 0;
  m = numel(I);
  ut = triu( ones(m-1,m-1) );
  tmp = (z-z0)/zR;
  g = 1 ./ sqrt( tmp.^2 + 1 );
  dz = z(2:end) - z(1:end-1);
  dz = [ dz; dz(end) ];
  fitf = 1/2 * I(1:m-1) .* gam(1:m-1) ./ dz(1:m-1) ...
         - 1/2 * I(m) * gam(m) / dz(m) ...
         + 1/2 * ut * ( ...
           I(1:m-1) ./ dz(1:m-1) .* gam(1:m-1) .* log( g(2:m)./g(1:m-1) ) ...
         );
  f = ut * I(1:m-1);
  figure;
  plot( z(1:m-1), f );
  hold on;
  plot( z(1:m-1), fitf, 'r', 'LineWidth', 2 );
  if numel( trueMu )
    trueI = dz.*mu2I( trueMu, z, z0, zR, muAlpha, muBeta, muL0 );
    truef = flipud(cumsum(flipud(trueI)));
    plot( z, truef, 'k', 'LineWidth', 2 );
    legend('data (b)','fit (A\gamma)','actual (A*true\gamma)');
  else
    legend('data (b)','fit (A\gamma)');
  end
  xlabel('Depth (mm)');
  ylabel('f');


end

