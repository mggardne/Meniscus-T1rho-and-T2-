function bbox = fnd_bboxc(mask,nrsls,iszs,ipad)
%FND_BBOXC Given a cell array containing logical masks for MRI image
%          slices, an array with the number of slices in each
%          compartment, and the size of the images (number of rows and
%          number of columns), finds the bounding box in pixels that
%          will contain all of the mask pixels.
%
%          BBOX = FND_BBOXC(MASK,NRSLS,ISZS) Given two-dimensional
%          logical masks with the first dimension being the image, and
%          the second dimension being slices in a cell array of masks
%          with the first index to lateral and medial and the second
%          index to the anterior and posterior, MASK, the number of
%          slices in each compartment, NRSLS, the size of the images
%          in pixels (number of rows and number of columns), ISZS,
%          finds the four-element array with the rectangular bounding
%          box in pixels that will contain all of the mask pixels,
%          BBOX.
%
%          BBOX = FND_BBOXC(MASK,NRSLS,ISZS,IPAD) Given the positive
%          integer, IPAD, pads all the edges of the bounding box by
%          IPAD pixels.  The default IPAD value is 1 to insure all of
%          pixels are within the bounding box.
%
%          NOTES:  None.
%
%          10-Aug-2022 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if nargin<3
  error([' *** ERROR in fnd_bboxc:  Three input variables are', ...
         ' required!']);
end
%
if nargin<4
  ipad = 1;
end
%
% Check For Valid IPAD
%
if isempty(ipad)||ipad<1||any(ipad>iszs)
  ipad = 1;
end
%
ipad = round(ipad);
%
% Get Starting Mask and Number of Top Level Masks
%
% Cell Levels in Variable "mask":
%   First Cell Level - Compartments:  1 - Lateral and 2 - Medial
%   Second Cell Level - Positions:  1 - Anterior and 2 - Posterior
%
msk = mask{1}{1}(:,1);  % First mask
%
nc = size(mask,1);      % Number of compartments (usually two [2])
%
% Combine All the Masks
%
for kl = 1:nc
   np = size(mask{kl},1);              % Number of positions (usually two [2])
   for ka = 1:np
      for ks = 1:nrsls(kl);
         msk = msk|mask{kl}{ka}(:,ks);
      end
   end
end
%
% Get Range of Pixels
%
msk2 = reshape(msk,iszs);
[i,j] = find(msk2);
imn = min(i);
imx = max(i);
jmn = min(j);
jmx = max(j);
%
% Get Bounding Box with Padding
%
bbox = [jmn-ipad jmx+ipad imn-ipad imx+ipad];
%
return