function [ flux ] = RiverFluxMap( lats , lons , V , mask , DXC , DYC , XC , YC)
%FluxMap Plots normally decaying flux distribution from each shelf
%   Uses ShelfGen to pick a number of "sources" which share the total
%   volume flux from the shelf. Each source produces its own normal
%   distribution of flux; these are summed to give a distribution which
%   hopefully reflects the geometry of the shelf.


nx = size(DXC,1);
ny = size(DXC,2);


flux = zeros(nx,ny);

for i = 1:length(lats)
    
    % find nearest grid position to each lat and lon value
    [latval latindex] = min(abs(YC(:)-lats(i)));
    [lat_row, latfrac] = ind2sub(size(YC),latindex);
    closestlat = YC(latindex);

    [lonval lonindex] = min(abs(XC(:)-180-lons(i)));
    [lonfrac, lon_col] = ind2sub(size(XC),lonindex);
    closestlon = XC(lonindex)-180;

    
    % make a temporary grid in metres using DXC and DYC
    X = zeros(nx,ny);
    Y = zeros(nx,ny);
    
    for i=1:nx-1
        X(i+1,:) = X(i,1) + DXC(i,1);
    end

    for i=1:ny-1
        Y(:,i+1) = Y(:,i) + DYC(1,i);
    end
    
    % grid wraps around at 0E so take this into account
    if(lonfrac < 1080)
        y = Y(lonfrac+1080,latfrac);
        x = X(lonfrac+1080,latfrac);
    else
        y = Y(lonfrac-1080,latfrac);
        x = X(lonfrac-1080,latfrac);
    end
    
    % sum flux from all sources; FluxGen generates for one source
    flux = flux + FluxGen(nx,ny,mask,x,y,V,-pi, 6, X,Y);
end


end

