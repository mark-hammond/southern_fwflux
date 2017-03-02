% distribute Antarctic runoff more realistically
% 2016 Gt/yr from Jacobs and Helmer 1992
% 2016 Gt/yr = 2016e12 kg/yr = 2016e9 m^3/yr = 63.9e3 m^3/s = 63.9 mSv
% increase linearly from continent to Polar Front,
% then decrease linearly to Subantarctic Front

clear all
cd /skylla/data/runoff

% load Polar and Subantarctic Front locations
load fronts/pf.txt
pf(find(pf(:,1)<0),1)=pf(find(pf(:,1)<0),1)+360;
load fronts/saf.txt
saf(find(saf(:,1)<0),1)=saf(find(saf(:,1)<0),1)+360;

% load coastline location
load tbase1
iy=1:30;
[c,h]=contour(lon,lat(iy),map(iy,:),[0 1e10]);
c=c(:,find(c(1,:)));  

% coarsen to 1 degree
fronts=nan*ones(360,3);
for i=1:360
 ix=find(c(1,:)>=i-1&c(1,:)<i);
 fronts(i,1)=mmax(c(2,ix));
 ix=find(pf(:,1)>=i-1&pf(:,1)<i);
 fronts(i,2)=mmin(pf(ix,2));
 ix=find(saf(:,1)>=i-1&saf(:,1)<i);
 fronts(i,3)=mmax(saf(ix,2)); 
end

% plot coastline and fronts
clf
plot(c(1,:),c(2,:),'g.',pf(:,1),pf(:,2),'b.',saf(:,1),saf(:,2),'r.')
hold on
plot(lon,fronts,'k')

% area of each grid cell
A=111177^2*ones(360,180);
for i=1:180;
 A(:,i)=A(:,i)*cos(deg2rad(i-90.5));
end

% weight starting from zero at continent, 1 at polar front,
% and back to zero at subantarctic front
W=zeros(360,180);
for i=1:360
 for j=1:50
  if lat(j)>fronts(i,1)&lat(j)<=fronts(i,2)
   W(i,j)=(lat(j)-fronts(i,1))/(fronts(i,2)-fronts(i,1));
  end
  if lat(j)>fronts(i,2)&lat(j)<=fronts(i,3)
   W(i,j)=(fronts(i,3)-lat(j))/(fronts(i,3)-fronts(i,2));
  end
 end
end

% distribute 63.9 mSv (2016e9 m^3/yr),
% compute freshwater flux in m/yr,
% and add to runoff file
F=W*2016e9/sum(sum(A.*W));
R=readbin('runoff-360x180x12.bin',[360 180 12]);
R(:,1:30,:)=0;
for m=1:12
 R(:,:,m)=R(:,:,m)+F;
end
writebin('runoff_antarctic_2016GT_360x180x12.bin',R);

% distribute 70 mSv (2208e9 m^3/yr),
% compute freshwater flux in m/yr,
% and add to runoff file
F=W*2208e9/sum(sum(A.*W));
R=readbin('runoff-360x180x12.bin',[360 180 12]);
R(:,1:30,:)=0;
for m=1:12
 R(:,:,m)=R(:,:,m)+F;
end
writebin('runoff_antarctic_2208GT_360x180x12.bin',R);

% distribute 120 mSv (3785e9 m^3/yr),
% compute freshwater flux in m/yr,
% and add to runoff file
F=W*3785e9/sum(sum(A.*W));
R=readbin('runoff-360x180x12.bin',[360 180 12]);
R(:,1:30,:)=0;
for m=1:12
 R(:,:,m)=R(:,:,m)+F;
end
writebin('runoff_antarctic_3785GT_360x180x12.bin',R);

% plot and compare old and new runoff
R1=readbin('runoff-360x180x12.bin',[360 180 12]);
R2=readbin('runoff_antarctic_3785GT_360x180x12.bin',[360 180 12]);
clf, subplot(211), mypcolor(lon,lat,mean(R1,3)');
caxis([0 .1]), thincolorbar, plotland,
title('old Antarctic runoff was 1.2 m/yr')
subplot(212), mypcolor(lon,lat,mean(R2,3)');
caxis([0 .1]), thincolorbar, plotland
title('new Antarctic runoff in m/yr')
print -djpeg AntarcticRuynoff
