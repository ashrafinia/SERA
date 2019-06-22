%cooc3d
%Authors: Carl Philips and Daniel Li (2008)
%http://facweb.cs.depaul.edu/research/vc/contact.htm
%
%Syntax
%[featureVector, coocMat] = cooc3d (data, 'distance', [1;2;4;8], ...
%   'direction', [0 1 0;1 1 0; 0 1 -1])
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
%featureVector = haralick values for each cube (this is what's used for
%                classification. Each row pertains to a different cube.
%coocMat = the Co-Occurrence matrices.
%          coocMat(y,x,direction,distance,cube number)

%featureVector(:,1:12) = the haralick features for distance 1, ...
%   direction 1;
%featureVector(:,13:24) = the haralick features for distance 1, ...
%   direction 2;
%featureVector(:,127:167) = the haralick features for distance 2, ...
%   direction 1;
%The haralick features used (in order) are:
%Energy, Entropy, Correlation, Contrast, Variance, SumMean, Inertia, 
%Cluster Shade, Cluster tendendy, Homogeneity,MaxProbability, 
%Inverse Variance.


%Designed and tested for cubes with axis 20 voxels long.


%function [featureVector,coocMat] = cooc3d (data,distance, directions)
function [featureVector,coocMat] = cooc3d (varargin)
featureVector=NaN;
coocMat= NaN;
%inputStr = {'Distance','Direction'};

%Default settings
distance = [1;2;4;8]; %more or fewer distances?
numLevels = 16;
numHarFeature=12; %changing this may break the harFeatures function

offSet = [0 1 0; -1 1 0; -1 0 0; -1 -1 0]; %2D Co-Occurrence directions
%o,45,90,135 degrees

    %the additional 9 directions that make 3D Co-Occurrence from 2D
dimension3 = [0 1 -1; 0 0 -1; 0 -1 -1; -1 0 -1; 1 0 -1; -1 1 -1; 1 -1 -1;
           -1 -1 -1; 1 1 -1];
offSet = cat(1,offSet,dimension3);

%checking inputs

data = varargin{1};
temp = size(data);
if size(temp)<3
    disp('Error: This program is designed for 3 dimensional data')
    return;
end
numInput = size(varargin,2);
for inputs =2:numInput
    temp = varargin{1,inputs};
    if ~ischar(temp)
        continue;
    end
    temp = upper(temp);
    switch (temp)
     
         case 'DIRECTION'
             temp2 = int8(varargin{1,inputs+1});
             if size(size(temp2),2) ~=2
                 disp('Error: Direction input is formatted poorly')
                 return;
             end
             if size(temp2,2) ~=3
                 disp(['Error: Incorrect number of columns in ' ... 
                     'direction variable'])
                 return;
             end
             if max(max(temp2))>1 | min(min(temp2))<-1
                 disp('Error: Direction values can only be {-1,0,1}')
                 return;
             end
             offSet = temp2;
        
        case 'DISTANCE'
            temp2 = int8(varargin{1,inputs+1});
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
    end
end

noDirections = size(offSet,1); %number of directions, currently 13
coocMat = zeros(numLevels, numLevels, noDirections, size(distance,2), ...
    size(data,4));
featureVector = zeros(size(data,4),noDirections*size(distance,1)*...
    numHarFeature);





for iteration=1:size(data,4) %each new cube
    tempVec = zeros(1,0);
    for dist =1:size(distance,1) %distance
        [harMat, coocMat(:,:,:,dist,iteration)] = graycooc3d(...
            data(:,:,:,iteration),distance(dist),numLevels,...
            numHarFeature,offSet); 
        temphar = zeros(1,0);
        
        %organizing the data so each cube's data is on a row
        for clicks =1:size(harMat,1)
            temphar = cat(2,temphar,harMat(clicks,:));
        end
        tempVec = cat(2,tempVec,temphar);
              
    end

    %produces a larger space to separate cubes
    %haralickMat = cat(1,haralickMat,space3);
    featureVector(iteration,:) = tempVec;

    
    disp(['completed cube number' num2str(iteration)])
end
return


%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%I = the 3D image matrix
%distance = a vector of the distances to analyze in
%numLevels = the number of graylevels to be used
%numHarFeature = the number of haralick features to compute
%offSet = a matrix of the directions to analyze in
%
%harMat = a matrix of the haralick features in the format harMat(direction,
%feature)
%coMat the Co-Occurrence matrices produced
function [harMat,coMat]= graycooc3d(I,distance,numLevels,numHarFeature,...
    offSet)

%**************Variable initialization/Declaration**********************
harMat =0;


noDirections = size(offSet,1); %number of directions, currently 13
coMat = zeros(numLevels,numLevels,noDirections);


%************************graylevel resizing*******************************
numLevels = numLevels-1; %don't touch. Logical adding issue.
minImage = min(min(min(I)));
I=I-(minImage);
min(min(min(I)));
maxImage = max(max(max(I)));
tempShift = double(maxImage)/double(numLevels);
I = floor(double(I)/double(tempShift));
I=I+1;
numLevels = numLevels+1; %don't touch. Logical adding issue.
if max(max(max(I))) > numLevels
    disp('Error is graylevel resizing.')
    disp('cooc3d.m');
    return
end


%**************************Beginning analysis*************************
%Order of loops: Direction, slice, graylevel, graylevel locations
for direction =1:noDirections %currently 13 (for the 3d image)

    tempMat = zeros(numLevels,numLevels,size(I,3));
    for slicej =1:size(I,3)
         for j=1:numLevels %graylevel
             
             %finds all the instances of that graylevel
            [rowj,colj] = find(I(:,:,slicej)==j);  

            %populating the Cooc matrix.
            for tempCount = 1:size(rowj,1) 
                rowT = rowj(tempCount) + distance*offSet(direction,1);
                colT = colj(tempCount) + distance*offSet(direction,2);
                sliceT = slicej+distance*offSet(direction,3);
                [I1, I2, I3] = size(I);
                if rowT <= I1 && colT <= I2 && sliceT <= I3
                    if rowT > 0 && colT > 0 && sliceT > 0
                        
                        %Error checking for NANs and Infinite numbers
                        IIntensity = I(rowT,colT,sliceT);
                        if ~isnan(IIntensity)
                            if ~isinf(IIntensity)
                                %Matlab doesn't have a ++ operator.
                                tempMat(j,IIntensity,slicej)= tempMat...
                                    (j,IIntensity,slicej)+1;
                            end
                        end
                    end
                end
            end
            
        end

    end
    for slicej =1:size(I,3)
        coMat(:,:,direction)= coMat(:,:,direction)+tempMat(:,:,slicej);
    end
end

%extracting the Haralick features from the Co-Occurrence matrices
harMat = harFeatures(coMat,numHarFeature);
return


%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%coMat = Co-occurrence matrices 2D stack upon eachother (i j k) k is number
%of directions analyzed. For 3d that's 13. Created in cooc3d.m
%numHarFeature is the number of variables you will be extracting.
%Haralick order
%Energy, Entropy, Correlation, Contrast, Variance, SumMean, Inertia, 
%Cluster Shade, Cluster tendendy, Homogeneity,MaxProbability, 
%Inverse Variance.
%harMat = matrix of the haralick features in the format harMat(direction,
%feature)
function [harMat]= harFeatures(coMat, numHarFeature)

%numHarFeature=12;
%numPosFeature=12; %If you add any more features bump this up.
numLevels = size(coMat,1); %number of graylevels
harMat = zeros(numHarFeature,size(coMat,3));
%%%%%%tempHarMat = zeros(numPosFeature,1);  %continue working here....
%tempCoMat=zeros(size(coMat,1),size(coMat,2));


for iteration = 1:size(coMat,3) %directions

    
%%%%%%%%%%%%%%%%%%%%Preparation

%%%%%%%determining various p values

    pij = sum(sum(coMat(:,:,iteration))); %already normalized
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
    tempInert=0;
    tempShade=0;
    tempTen=0;
    tempInVar=0;
    for j=1:numLevels
        for i=1:numLevels
            value = coMat(j,i,iteration);
            
            tempEnergy = tempEnergy+ value^2;
            if(value~=0) 
                tempEntropy = tempEntropy + (value * log10(value));
            end
            tempCorr = tempCorr+ ((i-mux)*(j-muy)*(value/(sigy*sigx)));
            n=(abs(i-j))^2;
            tempCont = tempCont+ value*n;
            tempGen = tempGen+ value/(1+abs(1-j));
            tempVar = tempVar + ((i - mux)^2)*value+((j-muy)^2)*value;
            tempMean = tempMean + (i+j)*(value);
            tempInert = tempInert+ (i-j)^2*(value);
            tempShade=tempShade+ ((i+j-mux-muy)^3)*(value);
            tempTen = tempTen+ (((i + j - mux - muy)^4) .* (value));
            if i~=j
                tempInVar=tempInVar+ value/(i-j)^2;
            end
        end
    end
    harMat(1,iteration)=tempEnergy;         %Energy
    harMat(2,iteration) = -tempEntropy;     %Entropy
    harMat(3,iteration)=tempCorr;           %Correlation
    harMat(4,iteration)=tempCont;           %Contrast
    harMat(5,iteration) = tempGen;          %Homogeneity
    harMat(6,iteration) = tempVar/2;        %Variance
    harMat(7,iteration)=tempMean/2;         %Sum Mean
    harMat(8,iteration)=tempInert;          %Inertia
    harMat(9,iteration)=tempShade;          %Cluster Shade
    harMat(10,iteration) = tempTen;         %Cluster Tendency
    harMat(11,iteration) = max(max(coMat(:,:,iteration))); %Max Probability
    harMat(12,iteration) = tempInVar;       %Inverse Variance
    
    clear 'tempEnergy' 'tempEntropy' 'tempCorr' 'tempCont' 'tempGen';
    clear 'tempVar' 'tempMean' 'tempInert' 'tempShade';
    clear 'tempTen' 'tempInVar';

end
%makes it so that rows are cases
harMat = harMat';
return