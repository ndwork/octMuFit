
function out = henyeyGreenstein
  clear; close all;

  theta = -pi:pi/100:pi;

  p1 = findP( theta, 0.3 );
  p2 = findP( theta, 0.7 );
  p3 = findP( theta, 0.9 );

  semilogy( theta, p1, 'LineWidth', 4 );
  hold all;
  semilogy( theta, p2, 'LineWidth', 4 );
  semilogy( theta, p3, 'LineWidth', 4 );
  set( gcf );

  legend( 'g=0.3', 'g=0.7', 'g=0.9' );

  ah = gca;   % get the current axes handle
  set(ah,'xtick',[]);
  set(ah,'xticklabel',[]);
  %set(ah,'ytick',[]);
  %set(ah,'yticklabel',[]);
  set(ah,'FontSize',16,'FontWeight','bold');
  
end

function p = findP( theta, g )
  p = 1/(4*pi) * (1-g^2) ./ ( 1 + g^2 - 2*g*cos(theta)).^(3/2);
  p = p / sum(p);
end