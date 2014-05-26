

function out = applyA( gamma, I, dz, g )
  % I is the intensity data
  % dz is a scalar representing the size of each pixel in the z (depth)
  %   direction
  % gamma the inverse attenuation coefficient matrix, the same size as I

  s = size(I);
  ndimsI = sum( s > 1 );

  switch ndimsI
    case 1
      M = numel(I);
      ut = triu( ones(M-1,M-1) );
      out = applyA1D( gamma, I, dz, g, ut );
    case 2
      out = applyA2D( gamma, I, dz, g );
    case 3
      out = applyA3D( gamma, I, dz, g );
    otherwise
      error('improper use of applyA');
  end

end



function out = applyA3D( gamma, I, dz, g )
  [M, N, K] = size(I);
  out = zeros( M, N, K );

  Lg = repmat( log( g(2:M)./g(1:M-1) ), [1 N K] );
  utArg = I(1:M-1,:,:) .* gamma(1:M-1,:,:) .* Lg;

  summed = sum( utArg, 1 );
  summedM = repmat( summed, [M-1, 1, 1] );
  cumSummed = cumsum( utArg, 1 );
  utOut = summedM - cumSummed + utArg;

  out(1:M-1,:,:) = 1/(2*dz) * ( ...
    I(1:M-1,:,:) .* gamma(1:M-1,:,:) ...
    - repmat( I(M,:,:) .* gamma(M,:,:), [M-1,1] ) ...
    + utOut );
end




function out = applyA2D( gamma, I, dz, g )
  [M, N] = size(I);
  out = zeros( M, N );

  Lg = repmat( log( g(2:M)./g(1:M-1) ), [1 N] );
  utArg = I(1:M-1,:) .* gamma(1:M-1,:) .* Lg;

  summed = sum( utArg, 1 );
  summedM = repmat( summed, [M-1, 1] );
  cumSummed = cumsum( utArg, 1 );
  utOut = summedM - cumSummed + utArg;

  out(1:M-1,:) = 1/(2*dz) * ( ...
    I(1:M-1,:) .* gamma(1:M-1,:) ...
    - repmat( I(M,:) .* gamma(M,:), [M-1,1] ) ...
    + utOut );
end


% function out = applyA2D( gamma, I, dz, g )
%   [M,N] = size(I);
%   out = zeros(M,N);
%   ut = triu( ones(M-1,M-1) );
% 
%   tmp = repmat( log( g(2:M)./g(1:M-1) ), [1 N] );
%   
%   out(1:M-1,:) = 1/2 * ( ...
%       I(1:M-1,:) .* gamma(1:M-1,:) ./ dz ...
%       - repmat( I(M,:) .* gamma(M,:) ./ dz, [M-1,1] ) ...
%       + ut * ( I(1:M-1,:) ./ dz .* gamma(1:M-1,:) .* tmp ) ...
%     );
% end

% function out = applyA2D( gamma, I, dz, g )
%   [M, N] = size(I);
%   out = zeros(M, N);
%   
%   ut = triu( ones(M-1,M-1) );
%   
%   for j=1:N
%     out(:,j) = applyA1D( gamma(:,j), I(:,j), dz, g, ut );
%   end
% end



function out = applyA1D( gamma, I, dz, g, ut )

  M = numel(I);

  out = zeros( M, 1 );
  out(1:M-1) = 1/2 * I(1:M-1) .* gamma(1:M-1) ./ dz ...
    - 1/2 * I(M) * gamma(M) / dz ...
    + 1/2 * ut * ( ...
      I(1:M-1) ./ dz .* gamma(1:M-1) .* log( g(2:M)./g(1:M-1) ) ...
    );

end

