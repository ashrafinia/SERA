function D = bwdistNew(ROIOnly)
%BWDIST Distance transform of binary image.
%   D = BWDISTROIOnlyBW) computes the Euclidean distance transform of the
%   binary image BW. For each pixel in BW, the distance transform assigns
%   a number that is the distance between that pixel and the nearest
%   nonzero pixel of BW. BWDIST uses the chessboard distance with 4
%   connectedness for 2D and 6 connectedness for 3D matrices. 
%
% -------------------------------------------------------------------------
% INPUTS:
% - ROIonly: A 1D, 2D or 3D matrix of 0 and 1. 
% -------------------------------------------------------------------------
% OUTPUTS:
% - D: Distance matrix with the same size as the input.
% -------------------------------------------------------------------------
% AUTHOR(S): 
% - Saeed Ashrafinia 
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: Dec 2017
% -------------------------------------------------------------------------
% --> Copyright (c) 2007-2017, Saeed Ashrafinia
%     All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%      
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.
% -------------------------------------------------------------------------

%Initialize
tmpimg = ROIOnly; 
nD = ndims(ROIOnly);
flg = 0; 

%% Take care of 1D ROIs
if isrow(ROIOnly) || iscolumn(ROIOnly) 
    tmpimg=repmat(ROIOnly(:),1,sum(ROIOnly(:)));
    flg = 1;
end

%% Set the connectedness
if nD == 2 
    con = 4; 
elseif nD == 3
    con = 6; 
else
    error('What is going on with the dimenstions?!!! It should be 1, 2 or 3.');
end  

%% Calculate distance to edges
D = zeros(size(tmpimg),'single'); 
for d = 1: max(size(tmpimg))
    edges = bwperim(tmpimg,con);
    if sum(edges(:)) == 0
        continue
    end
    D = D + d*edges; 
    tmpimg = tmpimg - edges;
    
end

%% Fix the 1D
if flg
    D = D(: , ceil(size(D,2)/2));
    if isrow(ROIOnly) 
        %% make sure to return the same array size as the input
        D = D';
    end
end


end
