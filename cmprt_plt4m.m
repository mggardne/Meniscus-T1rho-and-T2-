function cmprt_plt4m(v,mask,rsls,nrsls,idt,tcp,nps,mxtc,cmap,bx, ...
                     txt1,psnam)
%CMPRT_PLT4M  Plots T1rho/T2* values with the underlying images by slice 
%          within regions of interest (ROIs).
%
%          CMPRT_PLT4M(V,MASK,RSLS,NRSLS,IDT,TCP,NPS) Given a four-
%          dimensional matrix of T1/T2 intensities from a MRI image
%          volume, V, where the first two dimensions are an image, the
%          third dimension are the slices, and the fourth dimension are
%          the spin lock/echo times, three dimensional logical masks
%          with the first dimension being the image, the second
%          dimension being the superficial layer in the first column and
%          deep layer in the second column and the third dimension being
%          slices in a cell array of masks with the first index to the
%          lateral and medial compartments and the second index to the
%          anterior and posterior, MASK, a cell array with the slices
%          within each compartment, RSLS, the number of slices in each
%          compartment, NRSLS, index to the  spin lock/echo time to use
%          for plotting, IDT, cell array of T1rho/T2* values, TCP, and
%          the number of fitted pixels in each slice within the
%          compartments in a cell array, NPS, plots the T1rho/T2*
%          values with the underlying images by slice within the
%          regions of interest defined by the MASK.
%
%          CMPRT_PLT4M(V,MASK,RSLS,NRSLS,IDT,TCP,NPS,MXTC,CMAP,BX,TXT1)
%          Given the maximum plotting value for the color scale, MXTC,
%          a three color map, CMAP, an array defining a plot bounding
%          box in pixels, BX, and a text string for the first line of
%          the plot title, TXT1, plots the T1rho/T2* values with a
%          color maximum of MXTC using the color map, CMAP, within the
%          plot window BX and using TXT1 for the first line of the plot
%          title.  The default maximum value is 70.  The default color
%          map is gray for the image and jet for the T1rho/T2* values.
%          The default window is the whole image.  The default first
%          line title text is "Results Plot".
%
%          CMPRT_PLT4M(V,MASK,RSLS,NRSLS,IDT,TCP,NPS,MXTC,CMAP,BX,TXT1,
%          PSNAM) Given the name for a PS file, PSNAM, prints the plots
%          to the PS file.  By default, the plots are not printed.
%
%          NOTES:  1.  Plots the pixel output of cmprt_ana4m.m.  See
%                  cmprt_ana4m.m and mri_m_fit.m.
%
%                  2.  Based on cmprt_plt.m and cmprt_plt4.m, but
%                  modified by removing the combining of deep and
%                  superficial cartilage layers.  Adapted for use with
%                  plotting knee meniscus results and added a bounding
%                  box to zoom in on the meniscus.
%
%          10-Aug-2022 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if nargin<7
  error([' *** ERROR in cmprt_plt:  Seven input variables are', ...
         ' required!']);
end
%
if nargin<8||isempty(mxtc)
  mxtc = 70;
end
%
if nargin<9||isempty(cmap)
%
% Default Color Map
%
  gmap = gray(128);     % Gray color map for not cartilage
  jmap = jet(128);      % Jet color map for cartilage measures
  cmap = [gmap; jmap];
end
%
if nargin<10||isempty(bx)
  bx = [];
end
%
if nargin<11||isempty(txt1)
  txt1 = 'Results Plot';
end
%
if nargin<12||isempty(psnam)
  isave = false;
else
  isave = true;
end
%
% Compartment Labels
%
tcmprts = {'Lateral'; 'Medial'};
%
% Initialize Arrays
%
idxs = [2 2];           % Maximum indices for compartment and position
%
% Loop through Compartments
%
for kr = 1:2
%
   mskr = mask{kr};     % Mask for this compartment
   rsl = rsls{kr};      % Slices for this compartment
   nrsl = nrsls(kr);    % Number of slices in this compartment
   tcmprt = [tcmprts{kr} ' Compartment'];
%
% Loop through Slices
%
   for ks = 1:nrsl
%
      slk = rsl(ks);    % Slice
%
% Get Slice Image
%
      rimg = squeeze(v(:,:,slk,idt));   % T1/T2 data for slice and plot spin lock/echo time
      rimgr = rimg;     % Image for plotting results
%
% Scale T1rho/T2* Image to -mxtc to Zero
%      
      rimgr = rimgr-min(rimgr(:));
      imgmx = max(rimgr);
      rimgr = mxtc*rimgr./imgmx;
      rimgr = rimgr-(mxtc+0.01);
%
% Loop through Position
%
      for kb = 1:2
%
         mskb = mskr{kb};              % Mask for this position
         msk = mskb(:,ks);             % Mask for this slice
%          msk = squeeze(msk(:,1)|msk(:,2));  % Combine cartilage layers
%
         idx = sub2ind(idxs,kb,kr); % Index to T1rho/T2* results
         npsk = nps{idx};
         npsks = sum(npsk(1:ks));
         npsks = (npsks-npsk(ks)+1:npsks)';
%
         tcpk = tcp{idx};
%
         rimgr(msk) = tcpk(npsks);  % T1rho/T2* values
%
      end               % End of kb loop - positions loop
%
% Plot Slice
%
      figure;
      orient landscape;
%
      imagesc(rimgr,[-mxtc mxtc]);
      colormap(cmap);
      axis image;
      if ~isempty(bx)
        axis(bx);
      end
      axis off;
      title({txt1; ['Slice ' int2str(slk)]; tcmprt},'FontSize',16, ...
            'FontWeight','bold');
%
      hb = colorbar;
      set(hb,'Limits',[0 mxtc]);
%
      if isave          % Print plots
        if kr==1&&ks==1
          print('-dpsc2','-r600','-fillpage',psnam);
        else
          print('-dpsc2','-r600','-fillpage','-append',psnam);
        end
      end
%
   end                  % End of ks loop - slices loop
%
%    close all;
%
end                     % End of kr loop - compartments loop
%
return