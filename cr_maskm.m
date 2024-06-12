function mask = cr_maskm(xy,npx,iplt)
%CR_MASKM  Creates an image mask based on the two-dimensional (2-D)
%          coordinates of a line.  The line coordinates are assumed to
%          enclose an area.  This function is for making a mask for the
%          knee meniscus.
%
%          MASK = CR_MASKM(XY,NPX) given the two-dimensional
%          coordinates of a line enclosing an area in a two column
%          array (one coordinate in each column), XY, and the size of
%          the image for the mask in array, NPX, creates a logical
%          array mask, MASK.  Note:  If NPX has only one element, the
%          image is assumed to be square (symmetrical) (NPX by NPX).
%
%          MASK = CR_MASKM(XY,NPX,IPLT) if IPLT is true (or nonzero),
%          the triangular mesh of the segmentation is plotted in two
%          new figures.  These figures are closed after a pause.
%
%          NOTES:  1.  M-files in_tri2d.m, and mk_tri_m.m must be in
%                  the current directory or path.
%
%                  2.  For generating a mask of one meniscus
%                  segmentation on an individual slice.
%
%                  3.  The coordinates should enclose an area, but
%                  should not be closed (i.e., the first and last
%                  points should be close, but not the same
%                  coordinates).
%
%          21-Jul-2022 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<2)
  error(' *** ERROR in cr_maskm:  Two inputs are required!');
end
%
[npts,ncol] = size(xy);
%
if ncol~=2
  error([' *** ERROR in cr_maskm:  Two-dimensional input', ...
         ' coordinates must be in a two-column array!']);
end
%
ndim = size(npx(:),1);
if ndim>2||ndim<1
  error([' *** ERROR in cr_maskm:  Incorrect number of image ', ...
         'dimensions!']);
end
%
if ndim==1
  npx = [npx; npx];
end
%
if nargin<3
  iplt = false;
end
%
% Create Triangle Meshes for the Region of Interest
%
tri = mk_tri_m(xy,iplt);
%
% Find Pixels within the Region of Interest
%
mask = false(npx(1)*npx(2),1);
%
minr = floor(min(xy));
if any(minr(:)==0)      % Trap for zeros
  minr(minr(:)==0) = 1;
end
maxr = ceil(max(xy));
idx = minr(:,1):maxr(:,1);
idy = minr(:,2):maxr(:,2);
[xg,yg] = meshgrid(idx,idy);
xym = [xg(:) yg(:)];
%
idr = sub2ind([npx(1) npx(2)],xym(:,2),xym(:,1));
%
in_roi = in_tri2d(tri,xy,xym);
idr = idr(in_roi);
mask(idr) = true;
%
return