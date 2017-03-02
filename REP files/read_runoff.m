% open and read file in correct format
fileID = fopen('runoff_NVND.bin');
A=fread(fileID,'single','ieee-be');

% reshape to SOSE grid size
A = reshape(A,[2160 320]);

% plot (data is in m/yr)
imagesc(A)