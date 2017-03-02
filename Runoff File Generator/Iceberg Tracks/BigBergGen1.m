function [ output_args ] = BigBergGen( XC , YC )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


files = dir('*.dat');
bergs = zeros(0,1);

for file = files'
    
    try
    T=readtable(file.name,'Delimiter',' ','ReadVariableNames',false);
    pos =  table2array(T(:,[2 4]));
    end

    thisberg = pos;

    bergs = vertcat(bergs,thisberg);

    
end
bergs(logical(sum(bergs~=bergs,2)),:)=[];


result = roundtowardvec(bergs(:,1),YC(1,:));
result2 = roundtowardvec(bergs(:,2)+180,XC(:,1));

A=[result result2];

[Auniq,~,IC] = unique(A,'rows');
cnt = accumarray(IC,1);

cnt(:) = 0 + cnt(:);

cnt(:) = log(cnt(:));


[cnt Auniq];

latindex = zeros(length(cnt),1);
lonindex = zeros(length(cnt),1);

iceberggrid = zeros(size(XC,1),size(XC,2));

for i=1:length(cnt)

    latindex(i)=find(round(100*YC(1,:))==round(100*Auniq(i,1)));
    lonindex(i)= find(round(100*XC(:,1))==round(100*Auniq(i,2)));

    iceberggrid(lonindex(i),latindex(i)) = cnt(i);
    
end

iceberggrid =  circshift(iceberggrid,1080,1);


iceberggrid = hFacC(:,:,1) .* imgaussfilt(iceberggrid,20);



totalicebergflux = 500;

iceberggrid = totalicebergflux .* iceberggrid ./ sum(sum(iceberggrid));

imagesc(iceberggrid)

end

