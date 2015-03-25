
%function [I, z, z0, zR, alpha, beta, L0, trueMu] = makePhantom2D( ...
%  N, phantomType, z0, noiseProportion, muValues, muThicks, zR )


function [I, z, z0, zR, alpha, beta, L0, trueMu] = makePhantom2D( ...
  N, phantomType, varargin )

  M = 201;
  I = zeros(M, N);
  trueMu = zeros(M, N);
  for i = 1:N
    [I_col, z, z0, zR, alpha, beta, L0, trueMu_col] = makePhantom( ...
      phantomType, varargin{:} );
    I(:,i) = I_col;
    trueMu(:,i) = trueMu_col;
  end

end
