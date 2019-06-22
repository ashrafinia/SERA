function [coocMat,coocMat_sum,haralick_features]  = cooc3d_fast(varargin)

%Changed by Arman: to allow combining coMat before extracting har
%features (i.e. as opposed to extracting har features for every single
%distance and direction);
% Direction are summed; distances are cumulatively summed
%Authors: Carl Philips and Daniel Li (2008)
%http://facweb.cs.depaul.edu/research/vc/contact.htm
%
%Syntax
%[coocMat,coocMat_sum,features] = cooc3d(data, ROI, normflag, 'distance', [1;2;4;8], ...
%   'direction', [0 1 0;1 1 0; 0 1 -1],'NUMGRAY',16,'combine_all_coMat',1)
%
%Description
%reads in a vector of cubes and outputs a matrix of haralick features and
%the 3D Co-Occurrence matrices.
%
%Input:
%Data: a vector of cubes with the fourth dimension identifying the cube.
%data(:,:,:,1) = rand(20,20,20); %cube 1
%data(:,:,:,2) = rand(20,20,20); %cube 2
%
%Parameters:
%numGray: Integer indicating the number of graylevels to use when
%performing the graylevel resizing (rescaling).
%distance: a nx1 array of distances that will be used when analyzing the
%image. Default is [1,2,4,8];
%direction: a nx3 array of direction offsets in [row, column, vertical]
%format. The vertical value increases from top to bottom
%   the standard 2D directions
%   [0 1 0]    0 degrees
%   [-1 1 0]   45 degrees
%   [-1 0 0]   90 degrees
%   [-1 -1 0]  135 degrees
%
%   The additional 9 directions that make this a 3D Co-Occurrence ...
%       algorithm
%             horizontal, vertical
%   [0 1 -1]   0 degrees, 45 degrees
%   [0 0 -1]   straight up
%   [0 -1 -1]  0 degrees, 135 degrees
%   [-1 0 -1]  90 degrees, 45 degrees
%   [1 0 -1]   90 degrees, 135 degrees
%   [-1 1 -1]  45 degrees, 45 degrees
%   [1 -1 -1]  45 degrees, 135 degrees
%   [-1 -1 -1] 135 degrees, 45 degrees
%   [1 1 -1]   135 degrees, 135 degrees
%   Default is all 13 directions.
%
%Output:
% features:
% Depending on how combine_all_coMat is set, it is either (i) featureVector_overall = haralick values for combined
% coocurrence matrix along all distances and directions
%or is (ii) featureMatrix which is
%reshaped version of featureVector below)
%featureVector = haralick values for each cube (this is what's used for
%                classification. Each row pertains to a different cube.
%coocMat = the Co-Occurrence matrices.
%          coocMat(y,x,direction,distance,cube number)

%featureVector(:,1:13) = the haralick features for distance 1, ...
%   direction 1;
%featureVector(:,14:26) = the haralick features for distance 1, ...
%   direction 2;
%featureVector(:,etc) = the haralick features for distance 2, ...
%   direction 1;
%The haralick features used (in order) are:
%Energy, Entropy, Correlation, Contrast, Variance, SumMean, k-statistic,
%Cluster Shade, Cluster tendendy, Homogeneity,MaxProbability,
%Inverse Variance, Dissimilarity
%
%Designed and tested for cubes with axis 20 voxels long.
%
%normflag = grey level discretization flag.
%           0: takes min and max in each ROI
%           [normflag_min, normflag_max]: set common min and max

haralick_features=NaN;
coocMat= NaN;
coocMat_sum= NaN;

%Default settings (changed by input)
distance = [1;2;4;8]; %more or fewer distances?
numLevels = 16;
numHarFeature=13; %changing this may break the harFeatures function

% If set to 1, sum of all directions, and cumulative sum of cooc matrices will be considered with
% increasing distances (i.e. extending neighborhood); but if set to 0, it
% is better for detection of periodicity (especially using
% contrast/intertia measurement and k-statistic, as it looks at
% individual directions and distances (not sums of them).
combine_all_coMat=0;
accumulate_distance=1;
% Don't touch these above, since below, they are set to 0

offSet = [0 1 0; -1 1 0; -1 0 0; -1 -1 0]; %2D Co-Occurrence directions
dimension3 = [0 1 -1; 0 0 -1; 0 -1 -1; -1 0 -1; 1 0 -1; -1 1 -1; 1 -1 -1;
    -1 -1 -1; 1 1 -1];
offSet = cat(1,offSet,dimension3);
offSet=int32(offSet);

data = varargin{1};
temp = size(data);
if size(temp)<3
    disp('Error: This program is designed for 3 dimensional data')
    return;
end
if nargin>1
    ROI = varargin{2};
    temp = size(ROI);
    if size(temp)<3
        disp('Error: This program is designed for 3 dimensional data (ROI)')
        return;
    end
else
    ROI=ones(size(data));
end
data=data.*(ROI>0);
normflag = double(varargin{3});
numInput = size(varargin,2);
for inputs =4:numInput
    temp = varargin{1,inputs};
    if ~ischar(temp)
        continue;
    end
    temp = upper(temp);
    switch (temp)
        case 'DIRECTION'
            temp2 = int32(varargin{1,inputs+1});
            if size(size(temp2),2) ~=2
                disp('Error: Direction input is formatted poorly')
                return;
            end
            if size(temp2,2) ~=3
                disp(['Error: Incorrect number of columns in ' ...
                    'direction variable'])
                return;
            end
            if max(max(temp2))>1 || min(min(temp2))<-1
                disp('Error: Direction values can only be {-1,0,1}')
                return;
            end
            offSet = temp2;
            
        case 'DISTANCE'
            temp2 = int32(varargin{1,inputs+1});
            if size(size(temp2)) ~= 2
                disp('Error: Incorrect formatting of distance variable')
                return;
            end
            
            if sum(sum(size(temp2))) ~= max(size(temp2)+1)
                disp(['Error: Distance variable is to be a one ' ...
                    'dimensional array'])
                return;
            end
            distance = temp2;
            
        case 'NUMGRAY'
            temp2 = varargin{1,inputs+1};
            if temp2<1
                disp('The number of graylevels must be positive')
                return;
            end
            numLevels = uint16(temp2);
            
        case 'COMBINE'
            temp2 = varargin{1,inputs+1};
            combine_all_coMat=temp2;
            
        case 'ACCUMULATE'
            temp2 = varargin{1,inputs+1};
            accumulate_distance=temp2;
            
    end
end


noDirections = size(offSet,1); %number of directions, currently 13
coocMat = zeros(numLevels, numLevels, noDirections, length(distance), size(data,4));

%featureVector = zeros(numHarFeature*noDirections*length(distance),size(data,4));
featureMatrix = zeros(numHarFeature,noDirections,length(distance),size(data,4));
featureVector_overall = zeros(numHarFeature,size(data,4));
featureMatrix_overall = zeros(numHarFeature,length(distance),size(data,4));

% each iteration is for additional distances
for iteration=1:size(data,4) %each new cube
    tempVec = zeros(1,0);
    for dist =1:length(distance) %distance
        
        coocMat(:,:,:,dist,iteration) = graycooc3d_bis(data(:,:,:,iteration),...
            ROI(:,:,:,iteration),distance(dist),numLevels,offSet,normflag);
        
        %extracting the Haralick features from the Co-Occurrence matrices for each
        %direction separately
        if (~combine_all_coMat)
            harMat = harFeatures(coocMat(:,:,:,dist,iteration),numHarFeature);
        end
        temphar = zeros(1,0);
        if (~combine_all_coMat)
            %organizing the data so each cube's data is on a row
            for clicks =1:size(harMat,1)
                temphar = cat(2,temphar,harMat(clicks,:));
            end
            tempVec = cat(2,tempVec,temphar);
        end
    end
    
    if (combine_all_coMat)
        % Directions are summed; distances are cumulatively summed
        if length(distance)==1
            % 4th elements has a size of 1, anyways
            coocMat_sum=sum(sum(coocMat,3),4);
            harMat = harFeatures(coocMat_sum(:,:,1,1,iteration),numHarFeature);
            %harMat
            featureVector_overall(:,iteration) = harMat;
            haralick_features=featureVector_overall;
        else
            % use cumsum to cumulatively sum distances, while summing
            % directions
            if accumulate_distance
                coocMat_sum=cumsum(sum(coocMat,3),4);
            else
                coocMat_sum=sum(coocMat,3);
            end
            
            %        size(coocMat_sum)
            for dd=1:length(distance)
                harMat = harFeatures(coocMat_sum(:,:,1,dd,iteration),numHarFeature);
                %harMat
                featureMatrix_overall(:,dd,iteration) = harMat;
            end
            haralick_features=featureMatrix_overall;
        end
        
    else % individual distances and directions are considered
        hold=reshape(tempVec,numHarFeature,noDirections,length(distance));
        featureMatrix(:,:,:,iteration)=hold;
        haralick_features=featureMatrix;
        
    end
end
%save('zzzz.mat','data','ROI','coocMat');
end

%**********************************************************************************************************************
%
%I = the 3D image matrix
%distance = a vector of the distances to analyze in
%numLevels = the number of graylevels to be used
%offSet = a matrix of the directions to analyze in
%
%coMat the Co-Occurrence matrices produced

function coMat= graycooc3d_bis(I,ROI,distance,numLevels,offSet,normflag)

%**************Variable initialization/Declaration**********************

noDirections = size(offSet,1);
coMat = zeros(numLevels,numLevels,noDirections);

%************************graylevel resizing*******************************

if size(normflag,2)==2
    minImage=normflag(1);
    maxImage=normflag(2)-minImage;
    if minImage > min(I(find(ROI)))
        error('min too large in graycooc3d_bis.m');
    end
    if normflag(2) < max(I(find(ROI)))
        error('max too small in graycooc3d_bis.m');
    end
    I=I-(minImage);
else
    minImage = min(I(find(ROI)));
    I=I-(minImage);
    maxImage = max(I(find(ROI)));
end
tempShift = double(maxImage)/double(numLevels);
if tempShift==0
    I(find(ROI))=1;
else
    I = ceil(double(I)/double(tempShift));
    I(I==0)=1;
end
I=I.*(ROI>0);

if max(I(ROI(:)>0)) > numLevels
    disp('Error in graylevel resizing.')
    error('graycooc3d_bis.m');
end

[d1,d2,d3] = ind2sub(size(ROI),find(ROI));
% % d1range = min(d1)-1:max(d1)+1;
% % d2range = min(d2)-1:max(d2)+1;
% % d3range = min(d3)-1:max(d3)+1;
% Saeed: rewrite the following:
% d1range = max(1,min(d1)-1):max(d1)+1;
% d2range = max(1,min(d2)-1):max(d2)+1;
% d3range = max(1,min(d3)-1):max(d3)+1;
d1range = min(d1):max(d1);
d2range = min(d2):max(d2);
d3range = min(d3):max(d3);
Icropped = I(d1range, d2range, d3range);

for direction =1:noDirections
    coMat(:,:,direction) = grayLevelMatrixSingle(Icropped, numLevels, distance*offSet(direction,:), false);
end

end

%**********************************************************************************************************************
%
%coMat = Co-occurrence matrices 2D stack upon each other (i j k) k is number
%of directions analyzed. For 3d that's 13. Created in cooc3d.m
%numHarFeature is the number of variables you will be extracting.
%Haralick order
%Energy, Entropy, Correlation, Contrast, Homogeneity, Variance, SumMean, k-statistic,
%Cluster Shade, Cluster tendendy, MaxProbability,
%Inverse Variance.
%harMat = matrix of the haralick features in the format harMat(direction,
%feature)

function [harMat]= harFeatures(coMat, numHarFeature)

%numHarFeature=13;
%numPosFeature=13; %If you add any more features bump this up.
numLevels = size(coMat,1); %number of graylevels
harMat = zeros(numHarFeature,size(coMat,3));
%%%%%%tempHarMat = zeros(numPosFeature,1);  %continue working here....
%tempCoMat=zeros(size(coMat,1),size(coMat,2));


for iteration = 1:size(coMat,3) %directions
    
    
    %%%%%%%%%%%%%%%%%%%%Preparation
    
    %%%%%%%determining various p values
    % Normalization
    pij = sum(sum(coMat(:,:,iteration)));
    coMat(:,:,iteration)=coMat(:,:,iteration)./pij;
    
    tempmux=0;
    tempmuy=0;
    for j=1:numLevels
        for i=1:numLevels
            tempmux =  tempmux+(i*(coMat(j,i,iteration)));
            tempmuy =  tempmuy+(j*(coMat(j,i,iteration)));
        end
    end
    mux=tempmux; %mux
    muy=tempmuy;
    
    tempx=0;
    tempy=0;
    for j=1:numLevels
        for i=1:numLevels
            tempx = tempx+ (i-mux)^2*coMat(j,i,iteration);
            tempy = tempy+ (j-muy)^2*coMat(j,i,iteration);
        end
    end
    sigx=tempx; %sigx
    sigy=tempy;
    
    
    
    
    %Calculations
    tempEnergy =0;
    tempEntropy=0;
    tempCorr=0;
    tempCont=0;
    tempGen=0;
    tempVar=0;
    tempMean=0;
    %tempInert=0;
    tempShade=0;
    tempTen=0;
    tempInVar=0;
    tempDis=0;
    
    for j=1:numLevels
        for i=1:numLevels
            value = coMat(j,i,iteration);
            
            tempEnergy = tempEnergy+ value^2;
            if(value~=0)
                % Negative sign is below, later.
                tempEntropy = tempEntropy + (value * log10(value));
            end
            tempCorr = tempCorr+ ((i-mux)*(j-muy)*(value/(sigy*sigx)));
            n=(abs(i-j))^2;
            tempCont = tempCont+ value*n;
            tempDis = tempDis+ value*(abs(i-j));
            tempGen = tempGen+ value/(1+abs(i-j)); % some people
                                                   % do:  tempGen =
                                                   % tempGen+
                                                   % value/(1+abs(i-j)^2);
                                                   % Two types of
                                                   % homogenuity
                                                   % definition (1
                                                   % vs. 2)!!
            tempVar = tempVar + ((i - mux)^2)*value+((j-muy)^2)*value;
            tempMean = tempMean + (i+j)*(value);
            % tempInert = tempInert+ (i-j)^2*(value); same as tempCont  Contrast=Inertia
            tempShade=tempShade+ ((i+j-mux-muy)^3)*(value);
            tempTen = tempTen+ (((i + j - mux - muy)^4) .* (value));
            if i~=j
                tempInVar=tempInVar+ value/(i-j)^2;
            end
        end
    end
    
    % See Parkkinen 1990 (k-statistic; Cohen's kappa; is preferred method of measuring
    % periodicity than chi-squared; measures how diagonal co-occurrence
    % matrix is
    
    hold_matrix=coMat(:,:,iteration);
    Ci=sum(hold_matrix,1);
    Cj=sum(hold_matrix,2);
    %  size(Ci)
    %  size(Cj)
    
    P0=trace(hold_matrix);
    Pc=sum(Ci(:).*Cj(:)); % as one is row vector; one is column vector
    k_statistic=(P0-Pc)/(1-Pc);
    
    
    
    harMat(1,iteration)=tempEnergy;         %Energy (angular second moment)
    harMat(2,iteration) = -tempEntropy;     %Entropy
    harMat(3,iteration)=tempCorr;           %Correlation
    harMat(4,iteration)=tempCont;           %Contrast (is the same as inertia)
    harMat(5,iteration) = tempGen;          %(Local) Homogeneity
    harMat(6,iteration) = tempVar/2;        %Variance
    harMat(7,iteration)=tempMean/2;         %Sum Mean
    
    %%%harMat(8,iteration)=tempInert;       %Inertia (is the same as Contrast)
    harMat(8,iteration)=k_statistic;          % Cohen's kappa;
    % diagonality; measure of periodicity
    % Also called agreement
    
    harMat(9,iteration)=tempShade;          %Cluster Shade
    harMat(10,iteration) = tempTen;         %Cluster Tendency (or Prominence)
    harMat(11,iteration) = max(max(coMat(:,:,iteration))); %Max Probability
    harMat(12,iteration) = tempInVar;       %Inverse Variance
    harMat(13,iteration) = tempDis;         %Dissimilarity
    
    clear 'tempEnergy' 'tempEntropy' 'tempCorr' 'tempCont' 'tempGen';
    clear 'tempVar' 'tempMean' 'tempInert' 'tempShade';
    clear 'tempTen' 'tempInVar' 'tempDis';
    
end
%makes it so that rows are cases
harMat = harMat';

end

%**********************************************************************************************************************

function glcm = grayLevelMatrix(img, numLevels, offsets, symmetric)

% compute the glcm for all specified offsets
% Ivan Klyuzhin
%
% img is prepared before calling grayLevelMatrix
% rescaled into numLevels and zero outside ROI.
% Adapted by Stephan Blinder

numOffsets = size(offsets,1);
glcm = zeros(numLevels, numLevels, numOffsets);
for n = 1:numOffsets
    offset = offsets(n,:);
    glcm_tmp = grayLevelMatrixSingle(img, numLevels, offset, symmetric);
    glcm(:,:,n) = glcm_tmp;
end
glcm = mean(glcm,3);

end

%**********************************************************************************************************************

function glcmSingle = grayLevelMatrixSingle(img, numLevels, offset, symmetric)
% compute the gray level matrix for single offset
glcmSingle = zeros(numLevels, numLevels);

d1N = size(img, 1);
d2N = size(img, 2);
d3N = size(img, 3);

for d1 = 1:d1N
    for d2 = 1:d2N
        for d3 = 1:d3N
            v1 = img(d1,d2,d3);
            if ~v1, continue, end;
            t1 = d1 + offset(1);
            t2 = d2 + offset(2);
            t3 = d3 + offset(3);
            if (t1>0)&&(t1<=d1N)&&(t2>0)&&(t2<=d2N)&&(t3>0)&&(t3<=d3N)
                v2 = img(t1, t2, t3);
                if ~v2, continue, end;
                glcmSingle(v1,v2) = glcmSingle(v1,v2) + 1;
                if (symmetric)
                    glcmSingle(v2,v1) = glcmSingle(v2,v1) + 1;
                end
            end
        end
    end
end

end






% function metrics = getHaralikFeatures(img, roi, normflag)
%

%%% Arman changes (Sep 2015)
%%%[M0 M haralick_features_hold]=cooc3d_fast(matrix3D,ROI,0,'distance',distance,'direction',direction,'numgray',numgray,'COMBINE',combine_all_coMat,'ACCUMULATE',accumulate_distance);
%%%where accumulate_distance is either 0 or 1


% % Compute Haralik Features.
% % Based on Ivan Klyuzhin's GLCM and Arman Rahmim's Haralik features computation.
% % Struct has fields that correspond to metric names
% % ***Haralick feats (GLSD)***
% % metrics.harfeats
% %
% % Stephan blinder
% % May 8th, 2015
%
% %2D Co-Occurrence directions
% direction = [0 1 0; -1 1 0; -1 0 0; -1 -1 0];
% %direction = [1 0 0; 0 1 0; 1 1 0; 1 -1 0];
% %the additional 9 directions that make 3D Co-Occurrence from 2D
% dimension3 = [0 1 -1; 0 0 -1; 0 -1 -1; -1 0 -1; 1 0 -1; -1 1 -1; 1 -1 -1;-1 -1 -1; 1 1 -1];
% direction = cat(1,direction,dimension3);
% distance=[1;2;3;4;8;16];
% distance=[4];
% distance=[4;8];
%


% [M0 M haralick_features_hold]=cooc3d_bis(img,roi,normflag,'distance',distance,'direction',direction,'numgray',16,'combine_all_coMat',1);
% metrics.harfeats = haralick_features_hold;
% % 13x4 matrix is returned
% % 13 features - by - 6 distances
% % Feature order: (look into the end of "cooc3d.m")
% % 1 - Energy
% % 2 - Entropy
% % 3 - Correlation
% % 4 - Contrast
% % 5 - Homogeneity
% % 6 - Variance
% % 7 - Sum Mean
% % 8 - Cohen's kappa (agreement)
% % 9 - cluster shade
% % 10 - cluster tendency
% % 11 - max probability
% % 12 - inverse variance
% % 13 - Dissimilarity
%
%
% end
%**********************************************************************************************************************
