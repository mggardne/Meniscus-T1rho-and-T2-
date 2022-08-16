function bbox = fnd_bbox(masks,iszs,ipad)
%FND_BBOX  Given a two-dimensional logical array of masks for MRI
%          image slices, and the size of the images (number of rows and
%          number of columns), finds the bounding box in pixels that
%          will contain all of the masks' pixels.
%
%          BBOX = FND_BBOX(MASKS,ISZS) Given a two-dimensional logical
%          array of masks for MRI image slices, MASKS, with the first
%          dimension being the image and the second dimension being
%          slices, and the number of rows and number of columns in the
%          original image in a two element array, ISZS, finds the
%          four-element array with the rectangular bounding box in
%          pixels that will contain all of the mask pixels, BBOX.
%
%          BBOX = FND_BBOX(MASKS,ISZS,IPAD) Given the positive
%          integer, IPAD, pads all the edges of the bounding box by
%          IPAD pixels.  The default IPAD value is 1 to insure all of
%          pixels are within the bounding box.
%
%          NOTES:  None.
%
%          12-Aug-2022 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if nargin<2
  error([' *** ERROR in fnd_bbox:  Two input variables are', ...
         ' required!']);
end
%
if nargin<3
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
% Number of Slices
%
msk = masks(:,1);
ns = size(masks,2);      % Number of slices
%
% Combine All the Masks
%
for k = 2:ns
   msk = msk|masks(:,k);
end
%
% Get Range of Pixels
%
msk = reshape(msk,iszs);
[i,j] = find(msk);
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