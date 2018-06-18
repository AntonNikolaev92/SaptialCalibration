function [ H, I, nFrames ] = loadData( filePath )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

matfiles = dir(filePath);

str = [];
for i = 1:length(matfiles)
    str = [str, matfiles(i).name];
end
nFrames = 0;
imgfnames = regexp(str, 'USIMG\d*_\d*_\d*','match');
for i = 1:length(imgfnames)
    imgfname = imgfnames{i};
    posfname = replace(imgfname,'USIMG','USPOS');
    if ~isempty(posfname)
        nFrames = nFrames+1;
    end
end
if nFrames == 0, return; end;
load( [ filePath, '\', imgfname, '.mat'] );
[nSamples, ~] = size(USIMG);
nElements = 128;
H = zeros(4,4,nFrames);
I = zeros(nSamples, nElements, nFrames);
iFrame = 0;
for i = 1:length(imgfnames)
    imgfname = imgfnames{i};
    posfname = replace(imgfname,'USIMG','USPOS');
    if ~isempty(posfname)
        iFrame = iFrame+1;
        load( [ filePath, '\', imgfname, '.mat'] );
        I(:,:,iFrame) = USIMG(:,65:192); % correct to I(:,:,iFrame) for the normal acquisition
        load( [ filePath, '\', posfname, '.mat'] ); 
        T = eye(4); T(1:3,4) = USPOS.t(1:3,1)*1e-3;
        R = eye(4); R(1:3,1:3) = qGetR(USPOS.q(:,1));
        H(:,:,iFrame) = T*R;
    end
end

end

