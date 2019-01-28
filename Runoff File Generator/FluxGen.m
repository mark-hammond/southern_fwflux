function [ flux ] = FluxGen( nx , ny , mask , x , y , V , facing , opening , X , Y )
%FluxGen Generate a flux distribution for one source
%   Takes the position of a source,t he angles through which it is
%   allowed to propagate, and the overall mask. Produces a normally
%   distributed flux out into the allowed area.

% call MaskGen to exclude an additional area defined by the opening and
% facing angles
mask = MaskGen( facing , opening , X , Y , x , y , mask );

% produce an array of x and y distances from the sources
dist = zeros(nx,ny);

Xdist = abs(X-x);
Ydist = abs(Y-y);

Xmax = max(max(X));
halfXmax = Xmax / 2;

Xdistwrap = abs(Xmax - Xdist);      % remember that it wraps at 0E!


for i=1:(nx*ny)     % calculate absolute distances
    
    dist(i) = sqrt((Xdist(i)*Xdist(i) + Ydist(i)*Ydist(i)));
    
    % take wrapping into account
    if(dist(i) > halfXmax)
        dist(i) = sqrt((Xdistwrap(i)*Xdistwrap(i) + Ydist(i)*Ydist(i)));
    end
    
end

dist = reshape(dist,nx,ny);

% find minimum distance, i.e. closest location on grid
[M,I] = min(dist(:));
[I_y, I_x] = ind2sub(size(dist),I);


flux = zeros(nx,ny);

% plot normally distributed flux from this location
for m=1:nx
    for n=1:ny
        flux(m,n) = exp(-0.00000000005*dist(m,n)*dist(m,n));
    end
end

% apply mask
flux = flux .* mask;

% scale to total volume flux of source
scale = V / sum(sum(flux));
flux = flux .* scale;



end

