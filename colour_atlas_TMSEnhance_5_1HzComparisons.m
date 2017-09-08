function Vo = colour_atlas(P,x,outname)
% Colour in a brain region atlas using a vector of values for each region
% FORMAT:
% Q = colour_atlas(P,x,fname)
%
% INPUT:
% P          - filename of atlas image
% x          - vector of values to assign to atlas regions
% outname    - filename of output image
%
% OUTPUT:
% Vo         - Header information for output image (see 'spm_vol')
%
% Vector x must have the same length as the number of unique values in the atlas
% image, excluding zero.
% Values in x are assigned to atlas values in numerical order, so if atlas
% values are continuous:
%     region 1 = x(1), region 2 = x(2) etc.
% For an atlas with regions numbered 10, 15 & 5:
%     region 10 = x(2), region 15 = x(3), region 5 = x(1)
% This function assumes that all atlases have integer-numbered regions, and
% rounds the values loaded from the atlas image, because sometimes SPM
% incorrectly loads integer-valued atlases with non-integer values (perhaps
% a rounding error). If your atlas has non-integer values, modify the code.
%

% load data file
homedir = '/Users/simonwdavis/Dropbox/Duke/TMSEnhance/Scripts/';
cd(homedir)
load StimCompar.mat
penis = penis';
% define the background image--should have the same number of ROIs as your matrix
if nargin<1 || isempty(P)
    [P,sts] = spm_select(1,'image','Select atlas image ...');
    if sts == 0
        error('User quit')
    end
end
Vatlas = spm_vol(P);
Yatlas = spm_read_vols(Vatlas);
atvals = unique(Yatlas);
atvals = atvals(atvals~=0);

% change 'correls' to the name of your n (kinds of maps you want to make) x m (# of ROIs) matrix
for behav = 1:size(penis, 1)
    x = penis(behav, :);
    bname = sprintf('%d', behav);
    
    % change what's in '' to something meaningful
    [p,n,e] = spm_fileparts(P);
    if nargin<3 || isempty(outname)
        outname = fullfile(p,[n,'_ESA_', bname,e]);
    end
    
    Yo = zeros(size(Yatlas));
    for region = 1:length(atvals)
        % Use round in case integers in atlas are loaded incorrectly
        Yo(round(Yatlas)==atvals(region)) = x(region);
    end
    
    Vo = Vatlas;
    Vo = rmfield(Vo,{'pinfo','private'});
    Vo.fname = outname;
    Vo.dt = [spm_type('int16'),0]; % Should be signed data type
    spm_write_vol(Vo,Yo);
end
disp('Done');


