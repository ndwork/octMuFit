
function testModules
  clear; close all; rng(1);
  addpath(genpath(pwd))
 
  alpha = 0;    % Note, if we know true value than problem is better

  dataCase = 0;
  [I, z, dx, z0, zR, ~, ~, ~, trueMu ] = loadOctData( dataCase, false );
  dz = z(2) - z(1);
  
  % Extract smaller square of data
  if dataCase > 0
    firstRow = 138; %for dataCase 2
    firstCol = 168;
    numPix = 64;
    I = I(firstRow:firstRow+numPix-1 , firstCol:firstCol+numPix-1);
    z = z(firstRow:firstRow+numPix-1);
  end

  
  %mask = findNonZeroMus(I);
  mask = ones( size(I) );


  % Test applyA and applyAdjointA
  testApplyA_adjA(I, dz, z, z0, zR );


  % Test applyD and applyAdjointD
  testApplyD_adjD(I, dz);  
  
  % Test applyK and applyAdjointK
  testApplyK_adjK(I, z, dz, dx, z0, zR);


  % Test applyK and applyAdjointK in 1D
  testK_adjointK_1D(I, dz, z, z0, zR);  

  % Test Prox of G
  temp = rand(size(I)) - 0.5;
  prox_G = proxG(temp, mask);
  minProxG = min(prox_G(:));
  disp(['Minimium after proxG:' num2str(minProxG)]);
  
  % Test of proxFconj1D
  testProxFconj_1D(I);
  
  % Test proxF1D
  testProxF_1D(I, mask);
  
  % Compare makeA_1D to applyA
  testMakeA_1D(I, dz, z, z0, zR);
  
  % Compare makeK_1D to applyK_1D
  testMakeK_1D(I, dz, z, z0, zR);
  
  % Compare applyAdjointK1D to makeK' *x
  testApplyAdjointK1D(I, z, dz, z0, zR);

  % Compare Agamma to applyA
  testAGam(I, z, dz, z0, zR);

  % Compare objective function computation
  testObjFunc(I, mask, z, dz, z0, zR);

end


function testAGam( I2D, z, dz, z0, zR )
  I = I2D(:,1);
  tmp = (z-z0)/zR;
  g = 1 ./ sqrt( tmp.^2 + 1 );
  gam = rand(size(I));
  m = numel(I);
  A = makeA_1D(I, dz, g);
  
  Agam = A*gam;
  
  applyAgam = applyA(gam, I, dz, z, z0, zR);
  
  diff = norm(Agam - applyAgam);
  disp(['Difference between Agam and applyAgam:' num2str(diff)]);
end


function testObjFunc(I, mask, z, dz, z0, zR)
  I1D = I(:,1);
  mask1D = mask(:,1);
  b = makeb(I1D);
  tmp = (z-z0)/zR;
  g = 1 ./ sqrt( tmp.^2 + 1 );
  gamma = rand(size(I1D));
  D = makeD_1D(length(gamma), dz);
  A = makeA_1D(I1D, dz, g);
  
  eta = 100;

  objF = objFunction(gamma, I1D, mask1D, dz, g, D, eta);
  
  objF_matrx = (1/2)*norm( mask1D .* (A*gamma - b) ).^2 + ...
               eta*norm(mask1D.*(D*gamma), 1);
  
  relDiff = norm(objF - objF_matrx) / objF_matrx;
  disp(['Difference between objFunction and matrix objective Func:' num2str(relDiff)]);  
  
end


function testApplyAdjointK1D(I, z, dz, z0, zR)
  I1D = I(:,1);
  x1 = rand(size(I1D));
  x2 = rand(size(I1D));
  x = [x1; x2];
  
  tmp = (z-z0)/zR;
  g = 1 ./ sqrt( tmp.^2 + 1 );
  K = makeK_1D( I1D, dz, g );
  adjKx = K'*x;
  
  appAdjKx = applyAdjointK1D( x1, x2, I1D, z, dz, z0, zR );
  
  diff = norm(adjKx - appAdjKx);
  disp(['Difference between K''*x and applyAdjointK1D to x: ' num2str(diff)]);
  
end

function testMakeK_1D(I, dz, z, z0, zR)
  I1D = I(:,1);
  x = rand(size(I1D));

  tmp = (z-z0)/zR;
  g = 1 ./ sqrt( tmp.^2 + 1 );
  K = makeK_1D(I1D, dz, g);
  Kx = K*x;
  
  [appKx1, appKx2] = applyK1D(x, I1D, z, dz, z0, zR);
  appKx = [appKx1; appKx2];
  diff = norm(appKx - Kx);
  disp(['Difference between K*x and applyK_1D to x: ' num2str(diff)]);  
end

function testMakeA_1D(I, dz, z, z0, zR)
  I1D = I(:,1);
  x = rand(size(I1D));
  tmp = (z-z0)/zR;
  g = 1 ./ sqrt( tmp.^2 + 1 );
  A = makeA_1D(I1D, dz, g);
  Ax = A*x;
  
  appAx = applyA(x, I1D, dz, z, z0, zR);
  
  diff = norm(appAx - Ax);
  disp(['Difference between A*x and applyA_1D to x: ' num2str(diff)]);  
end

function testApplyA_adjA(I, dz, z, z0, zR )
  % check to see if <A*x,y> = <x,adjA*y>
  x = rand(size(I));
  y = rand(size(I));
  Ax = applyA( x, I, dz, z, z0, zR );
  Axy = Ax .* y;
  innerProd1 = sum( Axy(:) );
  adjAy = applyAdjointA( y, I, z, dz, z0, zR );
  xAdjAy = x .* adjAy;
  innerProd2 = sum( xAdjAy(:) );
  errorA = abs( innerProd1 - innerProd2 );
  disp(['Adjoint A error: ', num2str(errorA) ]);
end

function testApplyD_adjD(I, dz)
  dim = 1;
  x = rand( size(I) );
  y = rand( size(I) );
  Dx = applyD( x, dz, dim);
  adjDy = applyAdjointD(y, dz, dim);
  Dxy = Dx .* y;
  adjDyx = adjDy .* x;
  innerProd1 = sum( Dxy(:) );
  innerProd2 = sum( adjDyx(:) );
  errorA = abs( innerProd1 - innerProd2 );
  disp(['Adjoint D error: ', num2str(errorA) ]);
  
end

function testApplyK_adjK(I, z, dz, dx, z0, zR)
  y1 = rand( size(I) );
  y2 = rand( size(I) );
  y3 = rand( size(I) );
  x = rand ( size(I) );
  
  [Kx1, Kx2, Kx3] = applyK( x, I, z, dz, dx, z0, zR );
  tmp1 = Kx1.*y1 + Kx2.*y2 + Kx3.*y3;
  innerProd1 = sum( tmp1(:) );
  adjKy = applyAdjointK( y1, y2, y3, I, z, dz, dx, z0, zR );
  tmp2 = x .* adjKy;
  innerProd2 = sum( tmp2(:) );
  errorK = abs( innerProd1 - innerProd2 );
  disp(['Adjoint K error: ', num2str(errorK) ]);
end


function testK_adjointK_1D(I, dz, z, z0, zR)
  I1D = I(:,1);
  y1 = rand(size(I1D));
  mask1D = ones(size(I1D));
  x1D = rand(size(I1D));
  Ax = applyA(x1D,I1D, dz, z, z0, zR );
  innerProd1 = sum( Ax .* y1 );
  adjAy = applyAdjointA(y1,I1D,z,dz,z0,zR);
  innerProd2 = sum( x1D .* adjAy );
  errorK = abs( innerProd1 - innerProd2 );
  disp(['Adjoint A1D error: ', num2str(errorK) ]);
  
  % Test applyD and applyAdjointD in 1D
  Dx = applyD(x1D,dz,1 );
  innerProd1 = sum( Dx .* y1 );
  adjDy = applyAdjointD(y1,dz,1);
  innerProd2 = sum( x1D .* adjDy );
  errorK = abs( innerProd1 - innerProd2 );
  disp(['Adjoint D1D error: ', num2str(errorK) ]);
  
  % Test applyK1D and applyAdjointK1D
  y2 = rand(size(I1D));
  [Kx1, Kx2] = applyK1D( x1D, I1D, z, dz, z0, zR );
  tmp1 =  Kx1.*y1 + Kx2.*y2;
  innerProd1 = sum( tmp1(:) );
  adjKy = applyAdjointK1D( y1, y2, I1D, z, dz, z0, zR );
  tmp2 = x1D.*adjKy;
  innerProd2 = sum( tmp2(:) );
  errorK = abs( innerProd1 - innerProd2 );
  disp(['Adjoint K1D error: ', num2str(errorK) ]);  
end

function testProxFconj_1D(I)
  I1D = I(:,1)./max(I(:,1));
  b = makeb(I1D);
  
  m = length(b);
  x1 = rand(m,1);
  x2 = rand(m,1);
    
  sigma = 1;
  eta = 10;
  [prx1, prx2] = proxFconj1D(x1, x2, b, sigma, eta);
  
  % Determine proxF1conj using cvx
  cvx_begin quiet
    cvx_precision high
    variable u1(m,1)
    minimize ((1/2)*sum(u1.^2) + u1'*b + (1/2/sigma)*sum((u1 - x1).^2))
  cvx_end

  error1 = norm(u1 - prx1);
  disp(['Error with proxF1conj1D:' num2str(error1)]);
  
  % Determine proxF2conj using cvx
  cvx_begin quiet
    variable u2(m,1)
    minimize ((1/2/sigma)*sum((u2 - x2).^2))
    norm(u2, inf) <= eta
  cvx_end
  
  error2 = norm(u2 - prx2);
  disp(['Error with proxF2conj1D:' num2str(error2)]);
  
end



function testProxF_1D( I, mask )
  I1D = I(:,1)./max(I(:,1));
  mask1D = mask(:,1);
  b = makeb(I1D);
  
  m = length(b);
  x1 = rand(m,1);
  x2 = rand(m,1);
  
  lambda = 1;
  eta = 10;
  
  [prxF1, prxF2] = proxF1D(x1, x2, mask1D, b, lambda, eta);
  
  % Determine proxF1 using cvx
  cvx_begin quiet
    variable u1(m,1)
    minimize ((1/2)*sum((u1 - b).^2) + (1/2/lambda)*sum((u1 - x1).^2))
  cvx_end
  
  err1 = norm(prxF1 - u1);
  disp(['Error with proxF1_1D:' num2str(err1)]);
  
  % Determine proxF2 using cvx
  cvx_begin quiet
    variable u2(m,1)
    minimize (eta * norm(u2, 1) + (1/2/lambda)*sum((u2 - x2).^2))
  cvx_end
  
  err2 = norm(prxF2 - u2);
  disp(['Error with proxF2_1D:' num2str(err2)]);  
  
end
