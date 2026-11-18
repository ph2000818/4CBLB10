function [V] = CylinderVolume(Ca,Cyl)
% This function provides the cylinder volume as function of 
% Ca : Crankangle [degrees]
% Cyl :  a struct containing
%   Cyl.Stroke              : Stroke [m]
%   Cyl.Bore                : Bore [m]
%   Cyl.ConRod              : Connecting Rod length
%   Cyl.CompressionRatio    : Compession Ratio
%   Cyl.TDCangle            : Angle associated with the Top Dead Center

B   = Cyl.Bore;
S   = Cyl.Stroke;
cr  = Cyl.CompressionRatio;
r   = S/2;
l   = Cyl.ConRod;
%-------------------------------------------------------------------------------------------------------
Ap    = pi/4*B^2;
Vd    = Ap*S;
Vc    = Vd/(cr-1);
theta = Ca - Cyl.TDCangle;
V     = Vc + Ap*( r + l - ( r*cosd(theta) + sqrt(l^2 - r^2*sind(theta).^2) ) );




