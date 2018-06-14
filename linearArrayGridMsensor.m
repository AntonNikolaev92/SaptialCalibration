function [x, y, z] = linearArrayGridMsensor(pitch,nLines,fs,c,numSamps)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

[ xhlp, yhlp ] = linearArrayGrid(pitch,nLines,fs,c,0,numSamps);
y = -xhlp;
x =yhlp;
z = zeros(length(x),1);

end

