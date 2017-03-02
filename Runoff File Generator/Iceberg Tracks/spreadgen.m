
BW = zeros(2160,320);
BW(:,1:120) = edge(hFacC(:,1:120,1));


    

BW = hFacC(:,:,1) .* imgaussfilt(BW,20);



totalicebergflux = 500;

BW = totalicebergflux .* BW ./ sum(sum(BW));

imagesc(BW)

