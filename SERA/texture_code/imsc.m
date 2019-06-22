function imsc(FIG)

% This function plots FIG in a new figure with colorbar.

figure;
imagesc(squeeze(FIG));
colorbar;
axis image
title(inputname(1))
colormap('jet');