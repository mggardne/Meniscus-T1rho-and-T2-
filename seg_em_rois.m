%#######################################################################
%
% * SEGmentation to Eroded Meniscus Regions of Interest (ROIs) Program *
%
%          M-File which reads the registered MRI data and meniscus
%     segmentation CSV files to create masks for the meniscus regions
%     of interest in the lateral and medial compartments.  The masks
%     are saved in MAT files with the series number and ending in
%     "_emrois.mat."
%
%          The regions of interests are eroded by one pixel inferiorly
%     and superiorly to avoid including femoral or tibial cartilage and
%     any synovial fluid.
%
%     NOTES:  1.  The registered MRI MAT files must be in subject
%             directories starting with "PTOA.
%
%             2.  T1rho MAT files must start with "T1rho_S" and T2* MAT
%             files must start with "T2star_S".  See rd_m_dicom.m.
%
%             3.  M-files circ_plt.m, cr_mask2.m, cr_mask2f.m,
%             fnd_bbox.m, lsect2.m, lsect2a.m, lsect3.m, lsect4.m,
%             lsect5.m, midline.m,  mk2_tri_2d.m, mk2_tri_2df.m,
%             rd_roi6.m and tri2d.m must be in the current directory or
%             path.
%
%     12-Aug-2022 * Mack Gardner-Morse
%

%#######################################################################
%
% Setup Color Maps, Line Types and Color, and Legend Strings
%
gmap = gray(128);       % Gray color map for not cartilage
jmap = jet(128);        % Jet color map for cartilage
cmap = [gmap; jmap];
rmap = [[0.85 0 0]; [0 0 0.85]; gmap]; % Red/blue color map for ROIs
%
lt = ['g.-'; 'b.-'; 'r.-'; 'c.-'; 'y.-'; 'm.-']; % Line color and type
legds = {'ANT_LAT'; 'ANT_MED'; 'POST_LAT'; 'POST_MED'};
%
ipad = 25;              % Bounding box padding in pixels
%
% Compartment Labels
%
cmprt = {'Lateral'; 'Medial'};
%
% Structural Element to Erode Masks Inferiorly and Superiorly by a Pixel 
%
se = strel('line',3,90);               % Vertical line
%
% Directory with MAT Files as a String
%
dirstr = split(pwd,filesep);
dirstr = dirstr{end};
dirstr = split(dirstr,' ');
dirstr = dirstr{end};
dirstr = dirstr(1:3);   % Subject number as text
%
% T1rho
%
id5m = 1;               % Use first spin lock time for plots
rdir = '..\WB';         % Directory for T1rho segmentations
%
% Get T1rho MAT Files in Directory
%
d = dir('T1rho_S*.mat');
mnams = {d.name}';
idr = contains(mnams,'roi','IgnoreCase',true);
mnams = mnams(~idr);
nmat = size(mnams,1);
%
% Loop through T1rho MAT Files
%
for m = 1:nmat
%
   mnam = mnams{m};
   load(mnam);
   fs = ['S' snt];      % Series number prefaced with a 'S'
%
% Parse Series Text for Leg
%
   if strcmpi(st(1),'L')
     leg = 'L';
     lstr = 'left';
   else
     leg = 'R';
     lstr = 'right';
   end
%
% Get Segmentation File Name and Read ROIs
%
   rnam = dir(fullfile(rdir,[dirstr '_' leg '_sag_menisci_RHO_*.csv']));
   rnam = char({rnam.name}');
   if isempty(rnam)
     continue;
   end
   if size(rnam,1)>1
     error([' *** ERROR in seg_m_rois:  More than one segmentation', ...
            ' file for subject ' dirstr ' on the ' lstr ' leg!']);
   end
%
   rois = rd_roi6(fullfile(rdir,rnam),true);
   nrois = size(rois,1);
   idl = 1:nrois;
%
% Get MRI Slices
%
   rsl = [rois.imageno]';
   rsl = unique(rsl);        % Ensure unique slices in sorted order
   nrsl = size(rsl,1);
   rslmn = rsl(1);
   rslmx = rsl(nrsl);
%
% Loop through ROI Slices, Plot ROIs and Generate Masks for Each Slice
%
   pnam1 = [fs '_ROIs1e.ps'];          % ROI lines print file name
   pnam2 = [fs '_ROIs2e.ps'];          % ROI areas print file name
%
   mcoords = cell(nrsl,nrois);        % Meniscus coordinates
%
   iap = zeros(nrsl,nrois);            % -1 - anterior, 1 - posterior
   icmprt = zeros(nrsl,nrois);         % -1 - lateral, 1 - medial
   npxt = prod(iszs);                  % Total number of pixels in a slice
   maskm = false(npxt,nrsl,nrois);     % Masks for all meniscus segmentations
   masks = false(npxt,nrsl);           % Masks for meniscus slices
%
% Loop through Slices
%
   for k = 1:nrsl       % Loop through slices
%
% Get T1 Slice Image
%
      slk = rsl(k);     % Slice number
      img = squeeze(v(:,:,slk,id5m));
%
      figure;
      orient landscape;
      imagesc(img);
      colormap gray;
      axis image;
      axis off;
      hold on;
%
      lh = gobjects(nrois,1);          % Line graphic handles
      idxs = false(nrois,1);
%
% Get ROI Data for this Slice and Plot ROIs
%
      for l = 1:nrois   % Loop through ROIs
         idxr = rois(l).imageno==slk;
         if any(idxr)
           mcoords(k,l) = rois(l).data(idxr);   % Meniscus coordinates
%
           snam = rois(l).name;        % Get ROI name
%
           if strcmpi(snam(1),'A')
             iap(k,l) = -1;            % -1 - anterior
           else
             iap(k,l) = 1;             % 1 - posterior
           end
           if strcmpi(snam(end-1),'L')
             icmprt(k,l) = -1;         % -1 - lateral
           else
             icmprt(k,l) = 1;          % 1 - medial
           end
%
% Create Logical Masks for this Slice
%
           msk = cr_maskm(mcoords{k,l},iszs);    % Masks for all segmentations
           msk = reshape(msk,iszs);
           msk = imerode(msk,se);      % Erode mask by a pixel inferiorly/superiorly
           maskm(:,k,l) = msk(:);
           masks(:,k) = masks(:,k)|maskm(:,k,l);           % Masks for slices
%
           idx = iap(k,l)+icmprt(k,l)/2+2.5;
           lh(idx) = plot(mcoords{k,l}(:,1),mcoords{k,l}(:,2), ...
                          lt(idx,:));
           idxs(idx) = true;
         end            % End of if - any rois for this slice
      end               % End of rois loop
%
% Add Legend and Title, and Print Slice Plots
%
      bbx = fnd_bbox(masks(:,k),iszs,ipad);
      axis(bbx);
      legend(lh(idl(idxs)),legds(idl(idxs),:),'Interpreter','none');
      title({dirstr; [fs ' Slice ' int2str(slk)]; ...
            [cmprt{any(icmprt(k,:))+1} ' Compartment']}, ...
            'FontSize',16,'FontWeight','bold');
      if k==1
        print('-dpsc2','-r600','-fillpage',pnam1);
      else
        print('-dpsc2','-r600','-fillpage','-append',pnam1);
      end
%
% Plot ROIs
%
      cmx = max(img(:));
      img1 = img-cmx-1;
      dcmx = 16*cmx/128;
      maskms = squeeze(maskm(:,k,idxs));
      img1(maskms(:,iap(k,idxs)==-1)) = dcmx;         % Blue - Anterior
      img1(maskms(:,iap(k,idxs)==1)) = cmx-dcmx;      % Red - Posterior
%
      figure;
      orient landscape;
      imagesc(img1);
      colormap(cmap);
      caxis([-cmx cmx]);
      axis image;
      axis(bbx);
      axis off;
      title({dirstr; [fs ' Slice ' int2str(slk)]; ...
            [cmprt{any(icmprt(k,:))+1} ' Compartment']}, ...
            'FontSize',16,'FontWeight','bold');
%
      if k==1
        print('-dpsc2','-r600','-fillpage',pnam2);
      else
        print('-dpsc2','-r600','-fillpage','-append',pnam2);
      end
%
   end
%
% Get Bounding Box for All Slices
%
   bbox = fnd_bbox(masks,iszs,ipad);
%
% Get Compartment Identifier for Slices
%
   idl = sum(icmprt,2)<0;              % Index to lateral compartment slices
   rsll = rsl(idl);                    % Lateral compartment slices
   rslm = rsl(~idl);                   % Medial compartment slices
%
% Get Compartment and AP Specific Masks
%
   maskml = maskm(:,idl,:);            % Masks for lateral slices
   maskmm = maskm(:,~idl,:);           % Masks for medial slices
%
   idapl = sum(iap(idl,:));
   idal = idapl<0;                     % Column index to anterior lateral masks
   idpl = idapl>0;                     % Column index to posterior lateral masks
   maskmal = maskml(:,:,idal);         % Anterior lateral mask
   maskmpl = maskml(:,:,idpl);         % Posterior lateral mask
%
   idapm = sum(iap(~idl,:));
   idam = idapm<0;                     % Column index to anterior medial masks
   idpm = idapm>0;                     % Column index to posterior medial masks
   maskmam = maskmm(:,:,idam);         % Anterior medial mask
   maskmpm = maskmm(:,:,idpm);         % Posterior medial mask
%
% Get Compartment and AP Specific Meniscus Coordinates
%
   mcoordl = mcoords(idl,:);           % Meniscus coordinates for lateral slices
   mcoordm = mcoords(~idl,:);          % Meniscus coordinates for medial slices
%
   mcoordal = mcoordl(:,idal);         % Anterior lateral meniscus coordinates
   mcoordpl = mcoordl(:,idpl);         % Posterior lateral meniscus coordinates
%
   mcoordam = mcoordm(:,idam);         % Anterior medial meniscus coordinates
   mcoordpm = mcoordm(:,idpm);         % Posterior medial meniscus coordinates
%
% Save Masks, ROIS and Slice Information into MAT File
%
   savnam = [mnam(1:end-4) '_emrois.mat'];
   save(savnam,'bbox','cmprt','maskmal','maskmam', ...
        'maskmpl','maskmpm','masks','mcoordal','mcoordam', ...
        'mcoordpl','mcoordpm','rsl','rsll','rslm');
%
   close all;
%
end                     % End of m loop - MAT file loop
%
% T2star
%
rdir = '..\WB';         % Directory for T2* segmentations (same as T1rho)
id5 = repmat(4,6,1);    % Default echo time (5 ms) for segmentations
%
% Get T2* MAT Files in Directory
%
d = dir('T2star_S*.mat');
mnams = {d.name}';
idr = contains(mnams,'roi','IgnoreCase',true);
mnams = mnams(~idr);
nmat = size(mnams,1);
%
% Loop through T2* MAT Files
%
for m = 1:nmat
%
   mnam = mnams{m};
   load(mnam);
   fs = ['S' snt];      % Series number prefaced with a 'S'
%    if length(id5)>=m
%      id5m = id5(m);
%    else
     id5m = id5(1);
%    end
%
% Parse Series Text for Leg
%
   if strcmpi(st(1),'L')
     leg = 'L';
     lstr = 'left';
   else
     leg = 'R';
     lstr = 'right';
   end
%
% Get Segmentation File Name and Read ROIs
%
   rnam = dir(fullfile(rdir,[dirstr '_' leg '_sag_menisci_T2S_*.csv']));
   rnam = char({rnam.name}');
   if isempty(rnam)
     continue;
   end
   if size(rnam,1)>1
     error([' *** ERROR in seg_m_rois:  More than one segmentation', ...
            ' file for subject ' dirstr ' on the ' lstr ' leg!']);
   end
%
   rois = rd_roi6(fullfile(rdir,rnam),true);
   nrois = size(rois,1);
   idl = 1:nrois;
%
% Get MRI Slices
%
   rsl = [rois.imageno]';
   rsl = unique(rsl);        % Ensure unique slices in sorted order
   nrsl = size(rsl,1);
   rslmn = rsl(1);
   rslmx = rsl(nrsl);
%
% Loop through ROI Slices, Plot ROIs and Generate Masks for Each Slice
%
   pnam1 = [fs '_ROIs1e.ps'];          % ROI lines print file name
   pnam2 = [fs '_ROIs2e.ps'];          % ROI areas print file name
%
   mcoords = cell(nrsl,nrois);         % Meniscus coordinates
%
   iap = zeros(nrsl,nrois);            % -1 - anterior, 1 - posterior
   icmprt = zeros(nrsl,nrois);         % -1 - lateral, 1 - medial
   npxt = prod(iszs);                  % Total number of pixels in a slice
   maskm = false(npxt,nrsl,nrois);     % Masks for all meniscus segmentations
   masks = false(npxt,nrsl);           % Masks for meniscus slices
%
% Loop through Slices
%
   for k = 1:nrsl       % Loop through slices
%
% Get T2 Slice Image
%
      slk = rsl(k);     % Slice number
      img = squeeze(v(:,:,slk,id5m));
%
      figure;
      orient landscape;
      imagesc(img);
      colormap gray;
      axis image;
      axis off;
      hold on;
%
      lh = gobjects(nrois,1);          % Line graphic handles
      idxs = false(nrois,1);
%
% Get ROI Data for this Slice and Plot ROIs
%
      for l = 1:nrois   % Loop through ROIs
         idxr = rois(l).imageno==slk;
         if any(idxr)
           mcoords(k,l) = rois(l).data(idxr);   % Meniscus coordinates
%
           snam = rois(l).name;        % Get ROI name
%
           if strcmpi(snam(1),'A')
             iap(k,l) = -1;            % -1 - anterior
           else
             iap(k,l) = 1;             % 1 - posterior
           end
           if strcmpi(snam(end-1),'L')
             icmprt(k,l) = -1;         % -1 - lateral
           else
             icmprt(k,l) = 1;          % 1 - medial
           end
%
% Create Logical Masks for this Slice
%
           msk = cr_maskm(mcoords{k,l},iszs);    % Masks for all segmentations
           msk = reshape(msk,iszs);
           msk = imerode(msk,se);      % Erode mask by a pixel inferiorly/superiorly
           maskm(:,k,l) = msk(:);
           masks(:,k) = masks(:,k)|maskm(:,k,l);           % Masks for slices
%
           idx = iap(k,l)+icmprt(k,l)/2+2.5;
           lh(idx) = plot(mcoords{k,l}(:,1),mcoords{k,l}(:,2), ...
                          lt(idx,:));
           idxs(idx) = true;
         end            % End of if - any rois for this slice
      end               % End of rois loop
%
% Add Legend and Title, and Print Slice Plots
%
      bbx = fnd_bbox(masks(:,k),iszs,ipad);
      axis(bbx);
      legend(lh(idl(idxs)),legds(idl(idxs),:),'Interpreter','none');
      title({dirstr; [fs ' Slice ' int2str(slk)]; ...
            [cmprt{any(icmprt(k,:))+1} ' Compartment']}, ...
            'FontSize',16,'FontWeight','bold');
      if k==1
        print('-dpsc2','-r600','-fillpage',pnam1);
      else
        print('-dpsc2','-r600','-fillpage','-append',pnam1);
      end
%
% Plot ROIs
%
      cmx = max(img(:));
      img1 = img-cmx-1;
      dcmx = 16*cmx/128;
      maskms = squeeze(maskm(:,k,idxs));
      img1(maskms(:,iap(k,idxs)==-1)) = dcmx;         % Blue - Anterior
      img1(maskms(:,iap(k,idxs)==1)) = cmx-dcmx;      % Red - Posterior
%
      figure;
      orient landscape;
      imagesc(img1);
      colormap(cmap);
      caxis([-cmx cmx]);
      axis image;
      axis(bbx);
      axis off;
      title({dirstr; [fs ' Slice ' int2str(slk)]; ...
            [cmprt{any(icmprt(k,:))+1} ' Compartment']}, ...
            'FontSize',16,'FontWeight','bold');
%
      if k==1
        print('-dpsc2','-r600','-fillpage',pnam2);
      else
        print('-dpsc2','-r600','-fillpage','-append',pnam2);
      end
%
   end
%
% Get Bounding Box for All Slices
%
   bbox = fnd_bbox(masks,iszs,ipad);
%
% Get Compartment Identifier for Slices
%
   idl = sum(icmprt,2)<0;              % Index to lateral compartment slices
   rsll = rsl(idl);                    % Lateral compartment slices
   rslm = rsl(~idl);                   % Medial compartment slices
%
% Get Compartment and AP Specific Masks
%
   maskml = maskm(:,idl,:);            % Masks for lateral slices
   maskmm = maskm(:,~idl,:);           % Masks for medial slices
%
   idapl = sum(iap(idl,:));
   idal = idapl<0;                     % Column index to anterior lateral masks
   idpl = idapl>0;                     % Column index to posterior lateral masks
   maskmal = maskml(:,:,idal);         % Anterior lateral mask
   maskmpl = maskml(:,:,idpl);         % Posterior lateral mask
%
   idapm = sum(iap(~idl,:));
   idam = idapm<0;                     % Column index to anterior medial masks
   idpm = idapm>0;                     % Column index to posterior medial masks
   maskmam = maskmm(:,:,idam);         % Anterior medial mask
   maskmpm = maskmm(:,:,idpm);         % Posterior medial mask
%
% Get Compartment and AP Specific Meniscus Coordinates
%
   mcoordl = mcoords(idl,:);           % Meniscus coordinates for lateral slices
   mcoordm = mcoords(~idl,:);          % Meniscus coordinates for medial slices
%
   mcoordal = mcoordl(:,idal);         % Anterior lateral meniscus coordinates
   mcoordpl = mcoordl(:,idpl);         % Posterior lateral meniscus coordinates
%
   mcoordam = mcoordm(:,idam);         % Anterior medial meniscus coordinates
   mcoordpm = mcoordm(:,idpm);         % Posterior medial meniscus coordinates
%
% Save Masks, ROIS and Slice Information into MAT File
%
   savnam = [mnam(1:end-4) '_emrois.mat'];
   save(savnam,'bbox','cmprt','maskmal','maskmam', ...
        'maskmpl','maskmpm','masks','mcoordal','mcoordam', ...
        'mcoordpl','mcoordpm','rsl','rsll','rslm');
%
   close all;
%
end                     % End of m loop - MAT file loop
%
return