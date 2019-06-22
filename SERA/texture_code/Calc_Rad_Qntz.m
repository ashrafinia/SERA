clear
tic

%% Choose a dataset:
% Link to datasets. Each dataset should include folders containing two
% files: contour.mat and PETimg.mat. It's highly recommended to move
% dataset to a directory that only includes the dataset folders including
% the two files.
% PET_SPECT_CT : PET:0, SPECT:1, CT:2

% Renal Sestamibi Dataset:
% DataDir = 'C:\Users\ashra\OneDrive - Johns Hopkins University\Projects\Renal_Sestamibi\Renal_SPECT_CT\ExportToMatlab\'; PET_SPECT_CT = 0; disp('Selecting the RENAL SESTAMIBI dataset.');

% Head and Neck FDG Dataset
% DataDir = 'C:\Users\ashra\OneDrive - Johns Hopkins University\Projects\HeadnNeck_MIM\PatientData\'; PET_SPECT_CT = 0; 4disp('Selecting HEAD AND NECK PET dataset.');

% Cardiac SPECT Dataset:
% DataDir = 'C:\Users\ashra\OneDrive - Johns Hopkins University\Projects\CardiacSPECT\DataExportedtoMatlab\'; PET_SPECT_CT = 0;  disp('Selecting CARDIAC SPECT dataset.');

% Prostate PyL Dataset:
DataDir = 'C:\Users\Saeed\OneDrive - Johns Hopkins University\Projects\PyL_Prostate\ExportedImgContours_V2\'; PET_SPECT_CT = 0; dbs = 'PyL'; disp('Selecting Prostate PyL dataset.');
% DataDir = 'C:\Users\ashra\OneDrive - Johns Hopkins University\Projects\PyL_Prostate\ExportedImgContours_V2\'; PET_SPECT_CT = 0; dbs = 'PyL'; disp('Selecting Prostate PyL dataset.');
% DataDir = 'E:\OneDrive - Johns Hopkins University\Projects\PyL_Prostate\ExportedImgContours_V2\'; PET_SPECT_CT = 0; ROIperTmr = 8; disp('Selecting Prostate PyL dataset.');

% Radiomics features calcualtion codes folders
% addpath 'NIfTI_20140122'
addpath 'texture_code'

isotVoxSize = 4;
qntz = 'Uniform';

%% Retreiving folders inside the directory

ListDirName = dir(DataDir);
NumCases = size(ListDirName , 1) - 2;                % How many patients are there in the DIR

% Load DICOM info
% This file includes a 2-column vector involving voxel size in x and z
% dimensions. MAKE SURE THESE VALUES ARE IN MM.
DICOMinfo = xlsread([DataDir , '..\DICOMinfo\VoxelSizeInfo.xlsx']);

% Create a variable containing skipped cases:
% Flags: 1: less than 2 files inside dir, 2:no PETimg, 3:no contours,
% 4:conversion to uint16 error 5:empty ROI, 6:1D ROI, 7:image contains NaNs
% 8:contours.m is Empty 
SkippedCases = [];                           

% Check the number of files inside each directory
NumOKfolders = 0; % This counter only counts the folders that have exactly two files.
for casenum = 1:NumCases
    PatDir = dir([DataDir , char(cellstr(ListDirName(casenum + 2).name))]);         % get the info of the patient dir
    if (size(PatDir , 1) - 2) >= 2
        NumOKfolders = NumOKfolders +1;
        PatDirName(NumOKfolders) = cellstr(ListDirName(casenum + 2).name); %#ok<SAGROW> % save the patient directory name
    else
        warning(['"', char(cellstr(ListDirName(casenum + 2).name)) , '" has ',int2str((size(PatDir , 1) - 2)), ' files instead of 2.']);
        SkippedCases = cat(1,SkippedCases , { char(cellstr(ListDirName(casenum + 2).name)) , 0 , 1});
    end
end
% disp('*******************');
disp(['Total number of folders with at least 2 files: ',int2str(NumOKfolders), ' out of ',int2str(NumCases)]);

% % Check if we have DICOM data for all patients
% if size(DICOMinfo , 1) ~= NumOKfolders
%     error('Make sure number of rows in "DICOMinfo" is the same as OK patients.');
% end

%% Loop over every patient
for casenum = 1:NumOKfolders
    CurrentCaseName = PatDirName(casenum); %cellstr(ListDirName(casenum + 2).name); % save the patient directory name
    CurrentCaseDir = [DataDir , char(CurrentCaseName) , '\' ];
    
    %% Loading image and contours
    disp(['Loading image and contours for patient "', char(CurrentCaseName) , '", case # ',int2str(casenum),' at ' , num2str(toc,5)]);
    
    % Load image
    try
        load([CurrentCaseDir , 'PETimg']);
    catch
        warning(['No "PETimg" was found. Skipping case # ',char(CurrentCaseName)]);
        SkippedCases = cat(1,SkippedCases , { char(cellstr(ListDirName(casenum).name)) , 0 , 2});
        continue;
    end
    % Load contours
    try
        load([CurrentCaseDir , 'contours']);
    catch
        warning(['No "contours" was found. Skipping case # ',char(CurrentCaseName)]);
        SkippedCases = cat(1,SkippedCases , { char(cellstr(ListDirName(casenum).name)) , 0 , 3});
        continue;
    end
    if isempty(total)
        warning('"contours.m" is empty. Skipping the case.');
        SkippedCases = cat(1,SkippedCases , { char(cellstr(ListDirName(casenum).name)) , 0 , 8});
        continue;
    end
    
    % convert ROIs to uint16 for space conserving
    encryptedROIs = uint16(total);
    if sum(encryptedROIs(:)) ~= sum(total(:))
        encryptedROIs = uint32(total);
        if sum(encryptedROIs(:)) ~= sum(total(:))
            warning(['Error in converting ROI matrix to uint32. Skipping case # ',char(CurrentCaseName)]);
            SkippedCases = cat(1,SkippedCases , { char(cellstr(ListDirName(casenum).name)) , 0 , 4});
            continue;
        end
        disp('Error in converting ROI matrix to uint16. Convert to uint32 instead.');
    end
    clear total
    

    %% Decrypting contours
    [dxROI, dyROI, dzROI] = size(encryptedROIs);
    [dxImg, dyImg, dzImg] = size(vol_vals);
    R  = numel(de2bi(max(encryptedROIs(:)))); % ROIsEncrypt number of ROIs + 1
    
    ContourDecrypt  = zeros(dxROI,dyROI,dzROI,R-1,'uint8'); % store decrypted contours
    ROIsResized4D   = zeros(dxImg,dyImg,dzImg,R-1,'uint8');
    for k=1:dzROI
        roi = de2bi(encryptedROIs(:,:,k),R);
        for r = 1:R-1     % The first bit should be ignored!
            ContourDecrypt(:,:,k,r) = reshape(roi(:,r+1), dxROI,dyROI);
            ROIsResized4D(:,:,k,r) = imresize(ContourDecrypt(:,:,k,r) , dxImg/dxROI);
        end
    end
    clear ContourDecrypt
    
%     %% Resizing images and contours
%     if 4.69 ~= DICOMinfo(casenum ,1)
%         disp('.....Different voxel size; resizing the image......');
%         for k = 1:size(vol_vals,3)
%             imgresized(:,:,k) = imresize(vol_vals(:,:,k) , 4.69 / DICOMinfo(casenum ,1));
%             for r = 1:R-1
%                 roiresized(:,:,k,r) = ceil(imresize(ROIsResized4D(:,:,k,r) , 4.69 / DICOMinfo(casenum ,1)));
%             end
%         end
%         
%         % Display changing voxel size
%         if sum(size(ROIsResized4D)) ~= sum(size(roiresized))
%             disp('***3D image and ROIs resized****');
%         end
%         vol_vals = imgresized;
%         ROIsResized4D = roiresized;
%         %     [dxROI, dyROI, dzROI] = size(ROIsResized4D);
%         [dxImg, dyImg, dzImg] = size(imgresized);
%         
%     end
    
    
    
%     %% Plot the decrypted contours on top of CT
%     fig=figure; set(fig, 'Units', 'normalized', 'Position', [0,0,1,1]);
%     colormap('Jet');
%     MINC = min(vol_vals(:));
%     MAXC = max(vol_vals(:));
%     
%     seg2start = 2;  % 2 for all contours, R-1 for only epi and endo
%     if seg2start ~= 2, spi=[1 2]; else spi=[2 4]; end;
%     
%     for k=1:dzROI
%         TXT = ['Slice # ',int2str(k),', ROI # '];
%         for r = seg2start:R     % The first bit should be ignored!
%             subplot(spi(1),spi(2),r+1-seg2start)%-(R-3))
% %             imagesc(double(ROIsResized4D(:,:,3*(k-1)+1,r-1))*5000 + vol_vals(:,:,k));
%             imagesc(double(ROIsResized4D(:,:,k,r-1))*400 + vol_vals(:,:,k));
%             axis image
%             title([TXT,int2str(r-1)]);
%             caxis([MINC MAXC]);
%         end
%         pause(0.5);
%     end
    
    
    %% Start processing
    % We change data types from integer (even in DaTscan! and also ROI) to single
    % (for texture calculation code which requires singles or doubles)
    
    for roi = 1:R-1
        %% Prepare image and ROI
        imgvol=single(vol_vals);
        
        % Check if ROI is not empty
        try
            ROIbox = computeBoundingBox(ROIsResized4D(:,:,:,roi));
        catch
            if isempty(find(ROIsResized4D(:,:,:,roi)~=0)) %#ok<*EFIND>
                warning('ROI has no voxels, Skipping this ROI');
                SkippedCases = cat(1,SkippedCases , { char(cellstr(ListDirName(casenum).name)) , roi , 5});
                continue
            end
        end
        % Check that ROI is not one-dimensional
%         if sum((ROIbox(:,2)-ROIbox(:,1))==0) >= 2
%             warning('ROI is one-dimensional; skipping this ROI')
%             SkippedCases = cat(1,SkippedCases , { char(cellstr(ListDirName(casenum).name)) , roi , 6});
%             continue
%         end
        
        % Correct for SPECT images
        if PET_SPECT_CT == 1
            try
                ROI_ref = single(squeeze(ROIsResized4D(:,:,:,correction_roi)));
            catch
                error('An SPECT image requires a reference ROI for normalization');
            end
            mean_ROI_ref=sum(sum(sum(imgvol .* ROI_ref))) / sum(ROI_ref(:));
            imgvol = imgvol / mean_ROI_ref;
        end
        
        % Check image for NaNs
        if (isnan(sum(imgvol(:))))
            warning(['Vol image has NaNs. Case # ',int2str(casenum),', ROI # ',int2str(roi)]);
            SkippedCases = cat(1,SkippedCases , { char(cellstr(ListDirName(casenum).name)) , roi , 7});
            continue
        end
        loadedROI = single(squeeze(ROIsResized4D(:,:,:,roi))) ;
        
        %% Calculating features
        disp(['Calculating features for ROI # ',int2str(roi)]);
        
        [x_all, x_select]=main_texture_compute_qntz(imgvol,loadedROI , DICOMinfo(casenum ,1) , DICOMinfo(casenum ,2) , isotVoxSize , qntz); % assuming x and z dim are in mm, convert it to cm
        features_all(casenum,roi,:)=x_all; %#ok<SAGROW>
        features_select(casenum,roi,:)=x_select; %#ok<SAGROW>
        
    end
end


%%
size(features_all)
size(features_select)

%% save
time    = now; 
str     = datestr(time,0); str(str==':')='_'; str(str==' ')=',';
save([DataDir,'\..\CalculatedFeatures\Radiomics_',dbs,'_Isotrop_',num2str(isotVoxSize),'mm_Nbins_Qntz',qntz,'_',str,'.mat'] ,...
    'features_all' , 'features_select' , 'SkippedCases');

toc

