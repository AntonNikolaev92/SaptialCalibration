function [  ] = showImageInSpace( inImg, h, w, H )
% [  ] = showImageInSpace( inputImage, h, w, metrix )
% Description:
%   The function displays the gray-scale image inImg 3D space in
%   location , defined by the transformation metrix H. The intitial
%   location of the image is defined image high h and image width w.
%   Initially, the image is located as follow:
%   
%   0________________x
%   |           |w  
%   |           |
%   |   inImg   |
%   |           |
%   |           |
%  h|___________|
%   |
%   |y          
%
% Input Parameters:
%   inImg: input image to show in space
%   h    : high of the image
%   w    : width of the image
%   H    : trasformation metrix
%
% Output parameers: 


szI = size(inImg);
nSamples = szI(1);
nElements = szI(2);

linX = linspace(0, w,nElements)';
linY = linspace(0, h, nSamples)';

linZ = 0;

[ X, Y, Z ] = meshgrid(linX, linY, linZ);
Pf0 = ones(4,3);
Pf0(1,:) = [ min(linX), min(linX), max(linX) ];
Pf0(2,:) = [ min(linY), max(linY), min(linY) ];
Pf0(3,:) = [ 0, 0, 0 ];
Pf = H*Pf0;
P0 = ones(4, length(X(:)));
P0(1,:) = X(:);
P0(2,:) = Y(:);
P0(3,:) = Z(:);
P = H*P0;
[ A, B, C, D ] = getPlaneCoefficients( Pf(:,1), Pf(:,2), Pf(:,3) );
[ X, Y ] = meshgrid(linX, linY);
Z = zeros(size(X));
Idisp = inImg;
coefTreshold = 1e-5;

if abs(C) > coefTreshold
    X = reshape(P(1,1:nElements*nSamples), [ nSamples, nElements]);
    Y = reshape(P(2,1:nElements*nSamples), [ nSamples, nElements]);
    Z = -( A.*X + B.*Y + D )/C;
else
    if abs(A) > coefTreshold
        Z = reshape(P(3,1:nElements*nSamples), [ nSamples, nElements]);
        Y = reshape(P(2,1:nElements*nSamples), [ nSamples, nElements]);
        X = -( C.*Z + B.*Y + D )/A;
    else
        if abs(B) > coefTreshold
            X = reshape(P(1,1:nElements*nSamples), [ nSamples, nElements]);
            Z = reshape(P(3,1:nElements*nSamples), [ nSamples, nElements]);
            Y = -( A.*X + D.*Z + D )/B;
            
        end
    end
end

warp(X, Y, Z, Idisp);

end

