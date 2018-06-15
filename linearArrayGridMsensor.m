function [x, y, z] = linearArrayGridMsensor(pitch,nLines,fs,c,numSamps)
% [x, y, z] = linearArrayGridMsensor(pitch,nLines,fs,c,numSamps)
%   Detailed explanation goes here

[ xhlp, yhlp ] = linearArrayGrid(pitch,nLines,fs,c,0,numSamps);
y = -xhlp;
x = yhlp;
z = zeros(length(x),1);

end

