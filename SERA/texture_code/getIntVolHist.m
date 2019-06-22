function [textures] = getIntVolHist(IntsROI,ROIBox3D,BinSize,isReSeg,ResegIntrval,IVHconfig)
% -------------------------------------------------------------------------
% function [textures] = getGlobalTextures(ROIonly,Nbins)
% -------------------------------------------------------------------------
% DESCRIPTION:
% This function computes features related to the Intensity Volume Histogram
% ====> Make sure to specify the DataType. If not specified, Binsize will
% be considered as the 1000 bins for histogramming as directed by ISBI.
% -------------------------------------------------------------------------
% INPUTS:
% - ROIonly: Smallest box containing the ROI, with the imaging data ready
%            for texture analysis computations. Voxels outside the ROI are
%            set to NaNs. This should be the intesity ROI.
% - BinSize: scalar indicating the number of discretized bins
%           (or reconstruction levels of quantization).
% - DataType: 'PET': uses 1000 bins by default as suggested by ISBI.
%             'CT': uses the CT HU numbers as bins.
%             Otherwise: uses 1000 bins as suggested by ISBI.
% - minGL: minimum intensity within the ROI
%
% ** 'ROIonly' and 'levels' should be outputs from 'prepareVolume.m' **
% -------------------------------------------------------------------------
% OUTPUTS:
% - textures: Struture specifying the values of different IVH features as
% defined below.
% -------------------------------------------------------------------------
% AUTHOR(S): Saeed Ashrafinia
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: May 2017
% - Modification: July 2017
% -------------------------------------------------------------------------
% STATEMENT:
% This file is part of Radiomics Package by Saeed Ashrafinia, Rahmimlab.com
% --> Copyright (C) 2013-2017  Saeed Ashrafinia, Johns Hopkins University
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


%% Calculate hist and cum. hist
% Unit type = {0: Definite(PET,CT), 1:Arbitrary(MRI,SPECT. This is FNB)}, 
% Disc/Cont = {0: Discrete(for CT), 1:Continuous(for PET/CT, this is FBS)}, 
% Wb = {bin size}. 

% 
if IVHconfig(1) == 0
    %% Definite Units
    ROIarrayValid = squeeze(IntsROI(~isnan(IntsROI)));
    if IVHconfig(2) == 0 
        %% Discrete
        if isReSeg
            minGL = ResegIntrval(1);  maxGL = min(max(ROIarrayValid(:)), ResegIntrval(2));
        else
            minGL = min(ROIarrayValid(:)); maxGL = max(ROIarrayValid(:));
        end
        BinSize = 1;
        G = minGL:BinSize:maxGL;
        gamma = (G - min(G)) / range(G); % fractional grey level
        Hist = hist(ROIarrayValid,G);
        HistAc = cumsum(Hist);
        BinsCenters = G;
        V = 1-(HistAc ./ HistAc(end));
%         gamma = gamma; %((BinsCenters - BinsCenters(1)) / (BinsCenters(end)-BinsCenters(1)));
    elseif IVHconfig(2) == 1 
        %% Continuous (FBS)
        if isReSeg
            minGL = ResegIntrval(1);  maxGL = min(max(ROIarrayValid(:)), ResegIntrval(2));
        else
            minGL = min(ROIarrayValid(:)); maxGL = max(ROIarrayValid(:));
        end
        BinSize = IVHconfig(3);
        G = minGL : BinSize : (ceil(maxGL/BinSize)*BinSize);
        gamma = (G - min(G)) / range(G); % fractional grey level
        try
            Hist = hist(ROIarrayValid,G);
        catch
            disp('--Problem with IVH. Dividing the bin size by 10 to fix.'); 
            BinSize = BinSize/10;
            G = minGL : BinSize : (ceil(maxGL/BinSize)*BinSize);
            gamma = (G - min(G)) / range(G); % fractional grey level
            try
                Hist = hist(ROIarrayValid,G);
            catch
                disp('---It is even worse. Lets add 0.0001 to see if it fixes the histogram.');
                Hist = hist(ROIarrayValid,[G, G+0.0001]);
            end
        end
        HistAc = cumsum(Hist);
        BinsCenters = minGL + BinSize*((1:length(Hist))-0.5);
        V = 1-(HistAc ./ HistAc(end));   
    else
        error('Wrong IVH Config parameter!');  
    end
    
elseif IVHconfig(1) == 1
    %% Arbitrary Units
    ROIarrayValid = squeeze(IntsROI(~isnan(IntsROI)));
    if IVHconfig(2) == 0 
        %% FNB
        if isReSeg
            minGL = ResegIntrval(1);  maxGL = min(max(ROIarrayValid(:)), ResegIntrval(2));
        else
            minGL = min(ROIarrayValid(:)); maxGL = max(ROIarrayValid(:));
        end
        Ng = IVHconfig(3);
        Hist = hist(ROIarrayValid,Ng);
        HistAc = cumsum(Hist);
        BinsCenters = 1:length(Hist);
        BinSize = (maxGL-minGL)/Ng; 
        V = 1-(HistAc ./ HistAc(end));
        gamma = ((BinsCenters - BinsCenters(1)) / (BinsCenters(end)-BinsCenters(1)));
        
    else
        error('Wrong IVH Config parameter!');
    end
    
    % Discritize with 1000 FNB
elseif IVHconfig(1) == 2
    warning('Not a standard IVH Setting!! Make sure you know what you are doing.');
    ROIarrayValid = squeeze(IntsROI(~isnan(IntsROI)));
    BinSize = 1000;
    Hist = hist(ROIarrayValid,BinSize);
    HistAc = cumsum(Hist);
    V = 1-(HistAc ./ HistAc(end));
    BinsCenters = 1:BinSize;
    gamma = (((1:BinSize) - BinsCenters(1)) / (BinsCenters(end)-BinsCenters(1)));
    
    % Use existing intensities as bins (e.g. for CT)
elseif IVHconfig(1) == 3
%     warning('Not a standard IVH Setting!! Make sure you know what you are doing.');
    ROIarrayValid = squeeze(IntsROI(~isnan(IntsROI)));
    if isReSeg
        minGL = max(min(IntsROI(:)), ResegIntrval(1));  maxGL = min(max(IntsROI(:)), ResegIntrval(2));
    else
        minGL = min(ROIarrayValid(:)); maxGL = max(ROIarrayValid(:));
    end
    BinsCenters = minGL:1:ceil(maxGL);
    histraw = histcounts(ROIarrayValid,[BinsCenters, BinsCenters(end)+1]);
    HistAc = cumsum(histraw);
    V = 1 - (HistAc ./ HistAc(end));
    V = [1, V(1:end-1)];
    gamma = ((BinsCenters - BinsCenters(1)) / (BinsCenters(end)-BinsCenters(1)));
    
else
    error('Wrong IVH Config parameter!');
end


% *************************************************************************
%% Calculate Features:
% *************************************************************************

%% Volume at Intensity fraction
% 10%
V_10 = V(find(gamma>=0.1,1,'first'));

% 90%
V_90 = V(find(gamma>=0.9,1,'first'));


%%  Intensity at volume fraction
% 10%
I_10 = BinsCenters(min(length(BinsCenters) , find(V<=0.1,1,'first')));

% 90%
I_90 = BinsCenters(min(length(BinsCenters) , find(V<=0.9,1,'first')));


%% Volume at intensity fraction difference
V_10_90 = V_10 - V_90;


%% Intensity at volume fraction difference
I_10_90 = I_10 - I_90;


%% Area under IVH curve
AUC = trapz(gamma,V);


%% Wrap up
textures = [V_10; V_90; I_10; I_90; V_10_90; I_10_90; AUC ];

end


% % % % 
% % % %         if isReSeg
% % % %             minGL = ResegIntrval(1);  maxGL = min(max(ROIarrayValid(:)), ResegIntrval(2));
% % % %         else
% % % %             minGL = min(ROIarrayValid(:)); maxGL = max(ROIarrayValid(:));
% % % %         end
% % % %         BinSize = IVHconfig(3);
% % % %         G = minGL : BinSize : maxGL;% (ceil(maxGL/BinSize)*BinSize); 
% % % %         gamma = (G - min(G)) / range(G); % fractional grey level
% % % %         Hist = hist(ROIarrayValid,G);
% % % %         HistAcAll = cumsum(Hist);
% % % %         FirstNonZero = find(Hist,1); 
% % % %         HistAc = HistAcAll(FirstNonZero:end);
% % % %         BinsCenters = minGL + BinSize*((1:length(Hist))-0.5);
% % % %         V = 1-(HistAc ./ HistAc(end));   
% % % %         gamma = gamma(FirstNonZero:end);
% % % %         BinsCenters = BinsCenters(FirstNonZero:end);


% % % Using the global discritized ROI with global BinSize
% % if IVHconfig(1) == 0
% %     ROIarrayValid = squeeze(ROIBox3D(~isnan(ROIBox3D)));
% %     BinsCenters = [1, ((2:1:ceil(max(ROIarrayValid(:))/BinSize))-0.5)]*BinSize;
% %     histo = hist(ROIarrayValid,BinsCenters);
% %     histoc = cumsum(histo);
% %     V = 1-(histoc ./ histoc(end));
% %     g = ((BinsCenters - BinsCenters(1)) / (BinsCenters(end)-BinsCenters(1)));
% %     
% %     % using customized settings for IVH
% % elseif IVHconfig(1) == 1
% %     ROIarrayValid = squeeze(IntsROI(~isnan(IntsROI)));
% %     if IVHconfig(2) == 0 %FBS
% %         BinSize = IVHconfig(3);
% %         BinsCenters = [1, ((2:1:ceil(max(ROIarrayValid(:))/BinSize))-0.5)]*BinSize;
% %         histo = hist(ROIarrayValid,BinsCenters);
% %         histoc = cumsum(histo);
% %         V = 1-(histoc ./ histoc(end));
% %         g = ((BinsCenters - BinsCenters(1)) / (BinsCenters(end)-BinsCenters(1)));
% %         
% %     elseif IVHconfig(2) == 1 %FNB
% %         histo = hist(ROIarrayValid,BinSize);
% %         histoc = cumsum(histo);
% %         V = 1-(histoc ./ histoc(end));
% %         BinsCenters = 1:BinSize;
% %         g = (((1:BinSize) - BinsCenters(1)) / (BinsCenters(end)-BinsCenters(1)));
% %     else
% %         error('Wrong IVH Config parameter!');
% %     end
% %     
% %     % Discritize with 1000 FNB
% % elseif IVHconfig(1) == 2
% %     ROIarrayValid = squeeze(IntsROI(~isnan(IntsROI)));
% %     BinSize = 1000;
% %     histo = hist(ROIarrayValid,BinSize);
% %     histoc = cumsum(histo);
% %     V = 1-(histoc ./ histoc(end));
% %     BinsCenters = 1:BinSize;
% %     g = (((1:BinSize) - BinsCenters(1)) / (BinsCenters(end)-BinsCenters(1)));
% %     
% %     % Use existing intensities as bins (e.g. for CT)
% % elseif IVHconfig(1) == 3
% %     ROIarrayValid = squeeze(IntsROI(~isnan(IntsROI)));
% %     if isReSeg
% %         minGL = max(min(IntsROI(:)), ResegIntrval(1));  maxGL = min(max(IntsROI(:)), ResegIntrval(2));
% %     else
% %         minGL = min(ROIarrayValid(:)); maxGL = max(ROIarrayValid(:));
% %     end
% %     BinsCenters = minGL:1:ceil(maxGL);
% %     histraw = histcounts(ROIarrayValid,[BinsCenters, BinsCenters(end)+1]);
% %     histoc = cumsum(histraw);
% %     V = 1 - (histoc ./ histoc(end));
% %     V = [1, V(1:end-1)];
% %     g = ((BinsCenters - BinsCenters(1)) / (BinsCenters(end)-BinsCenters(1)));
% %     
% % else
% %     error('Wrong IVH Config parameter!');
% % end


%
%     if strcmp(DataType, 'PET')
% %         BinSize = 1000;
%         BinsCenters = [1, ((2:1:ceil(max(ROIarrayValid(:))/BinSize))-0.5)]*BinSize;
%         histo = hist(ROIarrayValid,BinsCenters);
%         histoc = cumsum(histo);
%         V = 1-(histoc ./ histoc(end));
%         g = ((BinsCenters - BinsCenters(1)) / (BinsCenters(end)-BinsCenters(1)));
%
%     elseif strcmp(DataType, 'CT') %|| strcmp(DataType, 'MRI')
%         if isReSeg
%             minGL = ResegIntrval(1); maxGL = ResegIntrval(2);
%         else
%             minGL = min(ROIarrayValid(:)); maxGL = max(ROIarrayValid(:));
%         end
%         BinsCenters = minGL:1:ceil(maxGL);
%         histraw = histcounts(ROIarrayValid,[BinsCenters, BinsCenters(end)+1]);
%         histoc = cumsum(histraw);
%         V = 1 - (histoc ./ histoc(end));
%         V = [1, V(1:end-1)];
%         g = ((BinsCenters - BinsCenters(1)) / (BinsCenters(end)-BinsCenters(1)));
%
%     else
%         BinSize = 1000;
%         histo = hist(ROIarrayValid,BinSize);
%         histoc = cumsum(histo);
%         V = 1-(histoc ./ histoc(end));
%         BinsCenters = 1:BinSize;
%         g = (((1:BinSize) - BinsCenters(1)) / (BinsCenters(end)-BinsCenters(1)));
%     end