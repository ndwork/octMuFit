
function K = makeK_1D( I, dz, g )

  M = numel(I);

  A = makeA_1D( I, dz, g );
%A = eye( numel(I) );
  D = makeD_1D( M, dz );
%D = zeros( M, M );

  K = [ A; D ];

end

