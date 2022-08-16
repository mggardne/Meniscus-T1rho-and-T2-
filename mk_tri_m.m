function tri = mk_tri_m(dat,iplt)
%MK_TRI_M  Makes a triangular mesh by using boundary line data from a
%          a two-dimensional digitized MRI slice.  The line data is
%          assumed to be an enclosed area.  This function is for making
%          a triangular mesh for the knee meniscus.
%
%          TRI = MK_TRI_M(DAT) given an array with two (2) columns with
%          the boundary line coordinate point data, DAT, returns a
%          three (3) column triangle connectivity matrix, TRI.
%
%          TRI = MK_TRI_M(DAT,IPLT) if IPLT is true (or nonzero), the
%          triangular mesh of the segmentation is plotted in two
%          new figures.  These figures are closed after a pause.
%
%          NOTES:  1.  For generating a triangular mesh of one meniscus
%                  segmentation on an individual slice.
%
%                  2.  The coordinates should form an enclosed area,
%                  but should not be closed (i.e., the first and last
%                  points should be close, but not the same
%                  coordinates).
%
%          21-Jul-2022 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<1)
  error(' *** ERROR in mk_tri_m:  No input coordinates!');
end
%
if nargin<2
  iplt = false;         % No plotting
end
%
% Check Input Coordinates
%
[npts,ncol] = size(dat);
%
if ncol~=2
  error([' *** ERROR in mk_tri_m:  Two-dimensional input', ...
         ' coordinates must be in a two-column array!']);
end
%
% Delaunay Triangulation
%
c = (1:npts)';
c = [c [2:npts 1]'];    % Constraints for Delaunay triangulation
%
dt = delaunayTriangulation(dat,c);
%
idin = isInterior(dt);  % Only interior triangles
tri = dt(idin,:);
%
% Plot Triangulations?
%
if iplt
%
  h1 = figure;
  orient landscape;
%
  nt = size(tri,1);
  xt = dat(:,1);
  yt = dat(:,2);
%
  plot(xt,yt,'k.','LineWidth',1,'MarkerSize',7);
  hold on;
  text(xt,yt,int2str((1:npts)'),'Color','k','FontSize',10);
%
  trimesh(tri,xt,yt);
  text(mean(xt(tri),2),mean(yt(tri),2),int2str((1:nt)'), ...
       'Color','r','FontSize',10);
%
  set(gca,'YDir','reverse');
  axis equal;
%
  h2 = figure;
  orient landscape;
  plot(xt,yt,'k.','LineWidth',1,'MarkerSize',7);
  hold on;
  text(xt,yt,int2str((1:npts)'),'Color','k','FontSize',10);
%
  xp = reshape(dat(tri,1),nt,3)';
  yp = reshape(dat(tri,2),nt,3)';
  xp = repmat(mean(xp),3,1)+0.75*(xp-repmat(mean(xp),3,1));
  yp = repmat(mean(yp),3,1)+0.75*(yp-repmat(mean(yp),3,1));
  patch(xp,yp,[1 0.7 0.7]);
  text(mean(xt(tri),2),mean(yt(tri),2),int2str((1:nt)'), ...
       'Color','r','FontSize',10);
%
  set(gca,'YDir','reverse');
  axis equal;
  pause;
  close(h1,h2);
%
end
%
return