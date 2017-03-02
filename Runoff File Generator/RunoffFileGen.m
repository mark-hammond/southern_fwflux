function [ totalfluxmyr totalvol a b c ] = RunoffFileGen( mask,DXC,DYC,XC,YC,smallfrac,basalvol,calvingvol,cmap)
%RunoffFileGen Produces a binary Antarctic freshwater runoff flux file
%   1) XC/YC: Grid cell centers
%   2) DXC/DYC: X/Y cell sizes
%   3) mask: Mask defining land (0) and ocean (1) on above grid
%   4) smallfrac: Fraction of "small" (normally distributed from coast
%   rather than following tracking database distribution) icebergs
%   5) basalvol/calvingvol: Total yearly volume of freshwater from basal
%   melting and iceberg calving, in Gt/yr
%   6) cmap: colourmap

% You'll need the "Bedmap2 ToolBox for Matlab" installed for ice shelf
% locations and shapes.

% This function should be run in a directory containing a folder 
% "Iceberg Tracks" with location files from 
% http://www.scp.byu.edu/data/iceberg/database1.html; I used QSCAT.

% The iceberg track files need to have their extensions changed to .dat 
% to work nicely with MATLAB's import function.

% It also need the riverrunoff.txt file (data from
% http://www.cgd.ucar.edu/cas/catalog/surface/dai-runoff/)



% calculate basal melting volume flux
basalflux = RunoffGen(mask,DXC,DYC,XC,YC,'jet');

% calculate iceberg  volume flux following tracking database distribution
bigbergflux = BigBergGen(XC,YC,mask,1);

% calculate iceberg volume flux normally distributed from coast
smallbergflux = SmallBergGen(mask,1);

% weight iceberg fluxes to match Rignot (2013) data
[newbig newsmall] = BigSmallWeight(bigbergflux,smallbergflux,smallfrac,mask);

% sort and rescale (a,b,c are individual fluxes in order above)
[totalvol a b c] = FreshFluxGen(basalflux,newbig,newsmall,basalvol,calvingvol,smallfrac);


% divide volume by area to get actual flux, then convert to m/yr
totalflux = totalvol ./ (DXC.*DYC);
totalfluxmyr = totalflux .* 1E9;

% add river runoff
figure
totalfluxmyr = totalfluxmyr + RiverRunoffGen(mask,DXC,DYC,XC,YC,cmap);


figure
imagesc((totalfluxmyr));
colormap(cmap);

% write final flux product to file
runofffile = totalfluxmyr(:);
fileID = fopen('runoff.bin','w');
fwrite(fileID,runofffile,'float32','ieee-be');
fclose(fileID);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ flux ] = RunoffGen( mask,DXC,DYC,XC,YC,cmap )
%RunoffGen Build a distribution of volume flux from basal melt
%   Reads in shelf data and passes it to fluxmap which generates the
%   distributions for that shelf. Sums these and prints shelf names as it
%   goes.

[shelfdata names] = xlsread('shelfdata.csv');       % read in shelf data

vols = shelfdata(:,1);
facingangles = shelfdata(:,2);
openingangles = shelfdata(:,3);                         % put data in arrays

flux = zeros(size(XC,1),size(XC,2));

for i=1:length(vols)        % sum fluxes due to each individual shelf

flux = flux + fluxmap(names{i} , vols(i) , mask , facingangles(i) , openingangles(i) ,DXC,DYC,XC,YC);
names{i};
colormap(cmap);

imagesc(flux)

pause(0.1);
%k = waitforbuttonpress     % can turn this on to wait for click each time

end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ iceberggrid ] = BigBergGen( XC , YC , mask , totalbigbergflux )
%BigBergGen Plot a volume flux distribution due to "big" icebergs
%   Reads in tracks of icebergs from the NASA SCP Iceberg Tracking Database
%   and plots them on the grid, weighting each grid point by the number of
%   occurences. Then applies a Gaussian filter to approximate a
%   climatology.


files = dir('Iceberg Tracks/*.dat');           % search for all .dat files
bergs = zeros(0,1);

for file = files'
                                % use 'try' as formatting can vary
    try
    T=readtable(file.name,'Delimiter',' ','ReadVariableNames',false);
    pos =  table2array(T(:,[2 4]));
    bergs = vertcat(bergs,pos);         % read in and add to array
    end
    
end

bergs(logical(sum(bergs~=bergs,2)),:)=[];   % remove null values

result = roundtowardvec(bergs(:,1),YC(1,:));
result2 = roundtowardvec(bergs(:,2)+180,XC(:,1));   % round to grid values

A=[result result2];

[Auniq,~,IC] = unique(A,'rows');
cnt = accumarray(IC,1);                 % count occurrences at each grid point

latindex = zeros(length(cnt),1);
lonindex = zeros(length(cnt),1);

iceberggrid = zeros(size(XC,1),size(XC,2));

for i=1:length(cnt) % fill grid with count values

    latindex(i)=find(round(100*YC(1,:))==round(100*Auniq(i,1)));
    lonindex(i)= find(round(100*XC(:,1))==round(100*Auniq(i,2)));

    iceberggrid(lonindex(i),latindex(i)) = cnt(i);
    
end

iceberggrid =  circshift(iceberggrid,1080,1);   % shift lon values because I'm dumb

iceberggrid = mask .* imgaussfilt(iceberggrid,20); % apply mask

iceberggrid = totalbigbergflux .* iceberggrid ./ sum(sum(iceberggrid)); % rescale to correct total flux

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ flux ] = SmallBergGen( mask , totalsmallicebergflux )
%SmallBergGen Generate rough climatology for small untrackable icebergs
%   Generates normally distributed volume flux from coastline for icebergs
%   which are too small for the tracking database. Standard deviation is
%   estimated to be 5E5m so that max extent is about 60S as in
%   Holland(2014).



flux = edge(mask);      % edge-detect coastline from mask
flux(:,121:end) = 0;    % get rid of non-Antarctic coastline

for j=1:320
for i =1:length(mask)
    newflux(i,j) = flux(i,j) ./ (sum(flux(i,1:120)));   % divide by y-count to avoid concentration on peninsulas etc
end
end

flux = mask .* imgaussfilt(double(newflux),50); % sigma of 5E5 is roughly 50 grid cells


flux = totalsmallicebergflux .* flux ./ sum(sum(flux)); % rescale to correct flux

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ newbigfluxdata newsmallfluxdata ] = BigSmallWeight( bigfluxdata , smallfluxdata , smallfrac , mask)
%BigSmallWeight Weight big (tracked) and small (normal) berg distributions
%   This splits the coastline into four sections which appear to send bergs
%   into different areas from the Tracking Database tracks. It weights
%   these so that they are in line with the total volume fluxes from
%   Rignot. There's probably a much better way of doing this but it's only
%   making small adjustments and isn't totally necessary.

% rignot calving flux data in four sectors
rigflux = [280 325 169 197];

% tracking database location fractions in four sectors
bigfrac = [0.582 0.295 0.102 0.022];

bigflux = zeros(4,1);
for i=1:4
    bigflux(i) = sum(rigflux) .* (1-smallfrac) .* bigfrac(i);
end

smallflux = zeros(4,1);
for i=1:4
    smallflux(i) = rigflux(i) - bigflux(i);
end

newbigfluxdata = bigfluxdata;
newsmallfluxdata = smallfluxdata;

% avert your eyes
%%%%
newbigfluxdata(1741:2040,:) = bigflux(1) .* bigfluxdata(1741:2040,:) ./ sum(sum(bigfluxdata(1741:2040,:)));

A=newbigfluxdata(2041:end,:);
B=newbigfluxdata(1:990,:);
C=vertcat(A,B);
C = bigflux(2) .* C ./ sum(sum(C));

newbigfluxdata(2041:end,:) = C(1:120,:);
newbigfluxdata(1:990,:) = C(121:end,:);

newbigfluxdata(991:1260,:) = bigflux(3) .* bigfluxdata(991:1260,:) ./ sum(sum(bigfluxdata(991:1260,:)));
newbigfluxdata(1261:1740,:) = bigflux(4) .* bigfluxdata(1261:1740,:) ./ sum(sum(bigfluxdata(1261:1740,:)));
%%%%

newsmallfluxdata(1741:2040,:) = smallflux(1) .* smallfluxdata(1741:2040,:) ./ sum(sum(smallfluxdata(1741:2040,:)));

A=newsmallfluxdata(2041:end,:);
B=newsmallfluxdata(1:990,:);
C=vertcat(A,B);
C = smallflux(2) .* C ./ sum(sum(C));

newsmallfluxdata(2041:end,:) = C(1:120,:);
newsmallfluxdata(1:990,:) = C(121:end,:);

newsmallfluxdata(901:1260,:) = smallflux(3) .* smallfluxdata(901:1260,:) ./ sum(sum(smallfluxdata(901:1260,:)));
newsmallfluxdata(1261:1840,:) = smallflux(4) .* smallfluxdata(1261:1840,:) ./ sum(sum(smallfluxdata(1261:1840,:)));

%%%%%%
% there we go

%newsmallfluxdata = (3.* newsmallfluxdata ./ 4) + (smallfluxdata./4);%%%%%%%%%%%%%%%%%%%%

totalsmallflux = sum(sum(newsmallfluxdata));
newsmallfluxdata = mask .* imgaussfilt(newsmallfluxdata,5);

% rescale
newsmallfluxdata = totalsmallflux .* newsmallfluxdata ./ sum(sum(newsmallfluxdata));


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ flux ] = FluxMap( name, totalvol , mask , facingangles , openingangles , DXC , DYC , XC , YC)
%FluxMap Plots normally decaying flux distribution from each shelf
%   Uses ShelfGen to pick a number of "sources" which share the total
%   volume flux from the shelf. Each source produces its own normal
%   distribution of flux; these are summed to give a distribution which
%   hopefully reflects the geometry of the shelf.

% call ShelfGen to get positions and fluxes of sources for this shelf
[lats,lons,vols] = ShelfGen( name , totalvol ); 

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
    
    V = vols(i);
    
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
    flux = flux + FluxGen(nx,ny,mask,x,y,V,facingangles , openingangles, X,Y);
end


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ freshflux basalflux bigbergflux smallbergflux] = FreshFluxGen( basalflux , bigbergflux , smallbergflux , totalbasalflux , totalbergflux , smallbergfrac )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

totalsmallbergflux = totalbergflux * smallbergfrac;

totalbigbergflux = totalbergflux - totalsmallbergflux;

basalflux = totalbasalflux .* basalflux ./ sum(sum(basalflux));
smallbergflux = totalsmallbergflux .* smallbergflux ./ sum(sum(smallbergflux));
bigbergflux = totalbigbergflux .* bigbergflux ./ sum(sum(bigbergflux));

freshflux = basalflux + smallbergflux + bigbergflux;


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ flux ] = RiverRunoffGen( mask,DXC,DYC,XC,YC,cmap )
%RiverRunoffGen aaa
%   aaa

fileID = fopen('riverrunoff.txt','r');
formatSpec = '%*f %f %f %f %f %f %*f %*f %*f %*f %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s';
A = textscan(fileID,formatSpec);
fclose(fileID);

B=cell2mat(A);

emptyrow = zeros(0,6);

[row col] = find(B(:,3)>-24.7);
B(row,:) = [];


vols = B(:,5);
lons = B(:,2);
lats = B(:,3);


flux = zeros(size(XC,1),size(XC,2));
%length(vols
for i=3:100        % sum fluxes due to each individual shelf

flux = flux + RiverFluxMap( lats(i) , lons(i) , vols(i) , mask ,DXC,DYC,XC,YC);
vols(i);
colormap(cmap);

imagesc(flux);

pause(0.1);
%k = waitforbuttonpress     % can turn this on to wait for click each time

end

flux = 1E9 .* flux ./ (DXC .* DYC);


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ lat , lon , vols ] = ShelfGen( name , TotalVol )
%ShelfGen Build a model of an ice shelf for use in FluxMap
%   Averages lat/lon bedmap2 outline values every n values and divides
%   total flux equally between these points. Passes these back as
%   "mini-shelves" which all produce freshwater, so that geometry of shelf
%   is represented.



n = 200; % average every n values

fristest = strcmp('fris',name); % check if it's the badly behaved filchner-ronne shelf
% (jagged west edge means distribution is skewed)

rosstest = strcmp('ross',name);% check if it's the badly behaved ross shelf
% (rignot says that all basal melting is from ross east)

if(fristest == 1)               % if it is then put in some better data for the north edge
    lat = [-75.46 -76.28 -77.3 -78.24 -79.26];
    lon = [-58.73 -54.48 -50.36 -36.37 -32.93];
else
    if(rosstest == 1)           % if it is then restrict to ross east
        lat = [-77.76 -78 -78.26 -78.5];
        lon = [-176.6 -177.8 -177 -169];
    else
        [h,lat,lon] = outlineashelf(name);  % otherwise call bedmap and average the outline
        lat = arrayfun(@(i) mean(lat(i:i+n-1)),1:n:length(lat)-n+1)';
        lon = arrayfun(@(i) mean(lon(i:i+n-1)),1:n:length(lon)-n+1)';
        
    end
end
    vols = ones(length(lat),1);
    
    vols = vols .* TotalVol ./ length(lat); % scale vols to correct total shelf volume

    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ newmask] = MaskGen( facingangle , openingangle , X , Y , x , y , mask )
%MaskGen Adjusts the mask applied to each shelf source
%   Takes in the original mask and grid, and the position, facing angle and
%   opening angle of the shelf. It excludes an area defined by these angles
%   to prevent e.g. Larsen shelves spreading to the other side of the
%   peninsula.


nx = size(X,1);
ny = size(X,2);

[THETA,RHO] = cart2pol(X-x,Y-y);                    % transform to polars
THETAround = (round(THETA,2) - facingangle); % rotate to facing angle
THETAround =  THETAround - pi / 2;

THETAround = wrapToPi(THETAround);

newmask =  - mask .* THETAround;                    % label mask grid with angles from shelf facing direction

%RHOtest = - mask .* RHO;

for i=1:nx
    for j=1:ny;
        if(abs(newmask(i,j)) > openingangle/2)
            newmask(i,j) = 0;                       % set any excluded areas (outside sector) to zero
        end
    end
end

newmask = newmask ./ newmask;                       % normalise back to 1 (for SOSE, top layer of hFacC so all 1)

newmask(isnan(newmask)) = 0;                        % get rid of any NaN values

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function newnums=roundtowardvec(X,roundvec,type)
%function newnums=roundtowardvec(X,[roundvec],[type])
%
% This function rounds number(s) toward given values. If more than one
% number is given to round, it will return the matrix with each rounded
% value, otherwise it will return the single rounded value. It will ignore
% NaNs and return them back with NaNs.
%
% Inputs: X: the number(s) that you want rounded
%
%         roundvec:(opt) the values to round X to. If none given, it will
%           default to -inf:1:inf (and use the built in functions).
%
%         type:(opt) specifies which kind of rounding you want
%           the function to use.
%
%           Choices are: 'round' - round to nearest value
%                        'floor' - round toward -Inf
%                        'ceil'  - round toward Inf
%                        'fix'   - round toward 0
%                        'away'  - round away from 0 (ceil if positive and floor if negative)
%                     (see help files for more clarity)
%
%           If no type is given, the function will default to rounding to
%           the nearest value.
%
% Outputs: newnums: rounded values, in same shape as X input matrix

% For some reason, loops seem to be faster than vectorizing this code

if nargin==0
    help roundtowardvec; %if nothing given, tell what to give
    return
elseif isempty(X)
    newnums=[]; %if given empty, return empty without going through whole script
    return
end
if ~exist('type','var') || isempty(type)
    type='round';  %%round to nearest value if not specified
end
if ~exist('roundvec','var') || isempty(roundvec)
    if strcmpi(type,'round')
        newnums=round(X);
        %to nearest integer
    elseif strcmpi(type,'away')
        newnums=ceil(abs(X)).*sign(X);
        %nearest integer away from 0
    elseif strcmpi(type,'fix')
        newnums=fix(X);
        %nearest integer toward 0
    elseif strcmpi(type,'floor')
        newnums=floor(X);
        %nearest integer toward -inf
    elseif strcmpi(type,'ceil')
        newnums=ceil(X);
        %nearest integer toward inf
    else
        error(sprintf('Round type not recognized. Options are:\n''round'' - round to nearest value\n''floor'' - round toward -Inf\n''ceil''  - round toward Inf\n''fix''   - round toward 0\n''away''  - round away from 0')) %#ok<SPERR>
    end
else
    %%make newnums size of X
    newnums=X;
    if strcmpi(type,'round') %to nearest value
        roundvec=reshape(unique(roundvec),1,[]);
        for i=numel(X):-1:1
            if ~any(X(i)==roundvec)
                DIFFs=abs(roundvec-X(i));
                if X(i)>=0
                    [~,ind]=min(DIFFs(:,end:-1:1));
                    newnums(i)=roundvec(length(DIFFs)-ind+1);
                elseif X(i)<0
                    [~,ind]=min(DIFFs);
                    newnums(i)=roundvec(ind);
                end
            end
        end
    elseif strcmpi(type,'fix') %to nearest value toward 0
        roundvec=reshape(unique([roundvec 0]),1,[]);
        for i=numel(X):-1:1
            if ~any(X(i)==roundvec)
                if X(i)>0
                    if X(i)>min(roundvec)
                        newnums(i)=roundvec(find(X(i)>roundvec,1,'last'));
                    else
                        newnums(i)=0;
                    end
                elseif X(i)<0
                    if X(i)<max(roundvec)
                        newnums(i)=roundvec(find(X(i)<roundvec,1,'first'));
                    else
                        newnums(i)=0;
                    end
                end
            end
        end
    elseif strcmpi(type,'ceil') %nearest value toward inf
        roundvec=reshape(unique(roundvec),1,[]);
        for i=numel(X):-1:1
            if ~isnan(X(i)) && ~any(X(i)==roundvec)
                if X(i)<max(roundvec)
                    newnums(i)=roundvec(find(X(i)<roundvec,1,'first'));
                else
                    newnums(i)=inf;
                end
            end
        end
    elseif strcmpi(type,'floor') %nearest value toward -inf
        roundvec=reshape(unique(roundvec),1,[]);
        for i=numel(X):-1:1
            if ~isnan(X(i)) && ~any(X(i)==roundvec)
                if X(i)>min(roundvec)
                    newnums(i)=roundvec(find(X(i)>roundvec,1,'last'));
                else
                    newnums(i)=-inf;
                end
            end
        end
    elseif strcmpi(type,'away') %nearest value away from 0
        roundvec=reshape(unique(roundvec),1,[]);
        for i=numel(X):-1:1
            if ~any(X(i)==roundvec)
                if X(i)>0
                    if X(i)<max(roundvec)
                        newnums(i)=roundvec(find(X(i)<roundvec,1,'first'));
                    else
                        newnums(i)=inf;
                    end
                elseif X(i)<0
                    if X(i)>min(roundvec)
                        newnums(i)=roundvec(find(X(i)>roundvec,1,'last'));
                    else
                        newnums(i)=-inf;
                    end
                elseif X(i)==0
                    DIFFs=abs(roundvec-X(i));
                    [~,ind]=min(DIFFs(:,end:-1:1));
                    newnums(i)=roundvec(length(DIFFs)-ind+1);
                end
            end
        end
    else
        error(sprintf('Round type not recognized. Options are:\n''round'' - round to nearest value\n''floor'' - round toward -Inf\n''ceil''  - round toward Inf\n''fix''   - round toward 0\n''away''  - round away from 0')) %#ok<SPERR>
    end
end
end













