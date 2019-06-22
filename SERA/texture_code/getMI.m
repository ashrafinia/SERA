function [metrics] = getMI(ROIbox) %,pixelW,sliceS)

% -------------------------------------------------------------------------
% function [metrics] = getMI(ROIonlyPET)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function computes moment invariants of an ROI.
% -------------------------------------------------------------------------
% INPUTS:
% - ROIbox: The smallest box containing the resampled 3D ROI, with the
%           imaging data ready for texture analysis computations. Voxels
%           outside the ROI are set to NaNs.
% -------------------------------------------------------------------------
% OUTPUTS:
% A list of 10 moment invariants features
% -------------------------------------------------------------------------
% AUTHOR(S):    Arman Rahmim
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: 2016
% -------------------------------------------------------------------------
% STATEMENT:
% This file is part of Radiomics Package by Arman Rahmim PhD, Rahmimlab.com
% --> Copyright (C) 2013-2017  Arman Rahmim, PhD; Johns Hopkins University
%
%    This package is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This package is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this package.  If not, see <http://www.gnu.org/licenses/>.
% -------------------------------------------------------------------------

mask = ~isnan(ROIbox); % Find mask covering the ROI
% numberVoxel = sum(mask(:));
%volume = numberVoxel * pixelW * pixelW * sliceS;
%metric=sum(ROIonlyPET(ROIonlyPET>=value))* pixelW * pixelW * sliceS;

[xdim, ydim, zdim]=size(ROIbox);

%sum(mask(:))
%'pause'
%pause

[X,Y,Z]=ndgrid(1:xdim,1:ydim,1:zdim);

min_image=min(ROIbox(mask));
mean_image=mean(ROIbox(mask));

hold_image=(ROIbox-min_image).*mask/((mean_image-min_image));

% Set all values outside of mask to 0 for MI calculation
hold_image(~mask)=0;

% Now, start calculating moment invariants
m000=sum(sum(sum(hold_image)));
om000=sum(sum(sum(mask)));

%pause

u000=m000;
ou000=om000;

% % % m100=sum(sum(sum(X.*hold_image)));
% % % m010=sum(sum(sum(Y.*hold_image)));
% % % m001=sum(sum(sum(Z.*hold_image)));

m100=sum(sum(sum(X.*hold_image)));
m010=sum(sum(sum(Y.*hold_image)));
m001=sum(sum(sum(Z.*hold_image)));

% % % om100=sum(sum(sum(X.*mask)));
% % % om010=sum(sum(sum(Y.*mask)));
% % % om001=sum(sum(sum(Z.*mask)));

om100=sum(sum(sum(X.*mask)));
om010=sum(sum(sum(Y.*mask)));
om001=sum(sum(sum(Z.*mask)));


x_mean=m100./m000;
y_mean=m010./m000;
z_mean=m001./m000;

ox_mean=om100./om000;
oy_mean=om010./om000;
oz_mean=om001./om000;


u200=sum(sum(sum((X-x_mean).^2.*hold_image)));
u020=sum(sum(sum((Y-y_mean).^2.*hold_image)));
u002=sum(sum(sum((Z-z_mean).^2.*hold_image)));

u110=sum(sum(sum((X-x_mean).*(Y-y_mean).*hold_image)));
u101=sum(sum(sum((X-x_mean).*(Z-z_mean).*hold_image)));
u011=sum(sum(sum((Y-y_mean).*(Z-z_mean).*hold_image)));

u300=sum(sum(sum((X-x_mean).^3.*hold_image)));
u030=sum(sum(sum((Y-y_mean).^3.*hold_image)));
u003=sum(sum(sum((Z-z_mean).^3.*hold_image)));

u210=sum(sum(sum((X-x_mean).^2.*(Y-y_mean).*hold_image)));
u201=sum(sum(sum((X-x_mean).^2.*(Z-z_mean).*hold_image)));
u120=sum(sum(sum((X-x_mean).*(Y-y_mean).^2.*hold_image)));
u102=sum(sum(sum((X-x_mean).*(Z-z_mean).^2.*hold_image)));
u021=sum(sum(sum((Y-y_mean).^2.*(Z-z_mean).*hold_image)));
u012=sum(sum(sum((Y-y_mean).*(Z-z_mean).^2.*hold_image)));
u111=sum(sum(sum((X-x_mean).*(Y-y_mean).*(Z-z_mean).*hold_image)));
%%%%
ou200=sum(sum(sum((X-ox_mean).^2.*mask)));
ou020=sum(sum(sum((Y-oy_mean).^2.*mask)));
ou002=sum(sum(sum((Z-oz_mean).^2.*mask)));

ou110=sum(sum(sum((X-ox_mean).*(Y-oy_mean).*mask)));
ou101=sum(sum(sum((X-ox_mean).*(Z-oz_mean).*mask)));
ou011=sum(sum(sum((Y-oy_mean).*(Z-oz_mean).*mask)));

ou300=sum(sum(sum((X-ox_mean).^3.*mask)));
ou030=sum(sum(sum((Y-oy_mean).^3.*mask)));
ou003=sum(sum(sum((Z-oz_mean).^3.*mask)));

ou210=sum(sum(sum((X-ox_mean).^2.*(Y-oy_mean).*mask)));
ou201=sum(sum(sum((X-ox_mean).^2.*(Z-oz_mean).*mask)));
ou120=sum(sum(sum((X-ox_mean).*(Y-oy_mean).^2.*mask)));
ou102=sum(sum(sum((X-ox_mean).*(Z-oz_mean).^2.*mask)));
ou021=sum(sum(sum((Y-oy_mean).^2.*(Z-oz_mean).*mask)));
ou012=sum(sum(sum((Y-oy_mean).*(Z-oz_mean).^2.*mask)));
ou111=sum(sum(sum((X-ox_mean).*(Y-oy_mean).*(Z-oz_mean).*mask)));


n200=u200./u000.^(2/3+1);
n020=u020./u000.^(2/3+1);
n002=u002./u000.^(2/3+1);

n110=u110./u000.^(2/3+1);
n101=u101./u000.^(2/3+1);
n011=u011./u000.^(2/3+1);

n300=u300./u000.^(2);
n030=u030./u000.^(2);
n003=u003./u000.^(2);

n210=u210./u000.^(2);
n201=u201./u000.^(2);
n120=u120./u000.^(2);
n102=u102./u000.^(2);
n021=u021./u000.^(2);
n012=u012./u000.^(2);
n111=u111./u000.^(2);

J1=n200+n020+n002;
Q=n200.^2+n020.^2+n002.^2+2*(n101.^2+n110.^2+n011.^2);
% Note that (as seen by Reiss; J2=1/2*(J1.^2-Q)
J2=n200.*n020+n200.*n002+n020.*n002-n101.^2-n110.^2-n011.^2;
J3=n200.*n020.*n002-n002.*n110.^2+2*n110.*n101.*n011-n020.*n101.^2-n200.*n011.^2;
B3=n300.^2+n030.^2+n003.^2+3*n210.^2+3*n201.^2+3*n120.^2+6*n111.^2+3*n102.^2+3*n021.^2+3*n012.^2;

on200=ou200./ou000.^(2/3+1);
on020=ou020./ou000.^(2/3+1);
on002=ou002./ou000.^(2/3+1);

on110=ou110./ou000.^(2/3+1);
on101=ou101./ou000.^(2/3+1);
on011=ou011./ou000.^(2/3+1);

on300=ou300./ou000.^(2);
on030=ou030./ou000.^(2);
on003=ou003./ou000.^(2);

on210=ou210./ou000.^(2);
on201=ou201./ou000.^(2);
on120=ou120./ou000.^(2);
on102=ou102./ou000.^(2);
on021=ou021./ou000.^(2);
on012=ou012./ou000.^(2);
on111=ou111./ou000.^(2);

oJ1=on200+on020+on002;
oQ=on200.^2+on020.^2+on002.^2+2*(on101.^2+on110.^2+on011.^2);
oJ2=on200.*on020+on200.*on002+on020.*on002-on101.^2-on110.^2-on011.^2;
oJ3=on200.*on020.*on002-on002.*on110.^2+2*on110.*on101.*on011-on020.*on101.^2-on200.*on011.^2;
oB3=on300.^2+on030.^2+on003.^2+3*on210.^2+3*on201.^2+3*on120.^2+6*on111.^2+3*on102.^2+3*on021.^2+3*on012.^2;

metrics=[J1 Q J2 J3 B3 oJ1 oQ oJ2 oJ3 oB3];

end