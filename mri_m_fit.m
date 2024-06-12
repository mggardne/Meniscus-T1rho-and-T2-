%#######################################################################
%
%                     * MRI Meniscus FIT Program *
%
%          M-File which reads the registered MRI data and segmentation 
%     MAT files and fits a monoexponential to the MRI data as a function
%     of spin lock or echo times where T1rho or T2* are the time
%     constants of the fits.  Resulting T1rho and T2* values and summary
%     statistics are written to the MS-Excel spreadsheet,
%     mri_m_fit.xlsx, in the "Results" directory.
%
%     NOTES:  1.  Data MAT files must be in subject directories 
%             "..\PTOA\0*" where "*" is the subject number.
%
%             2.  T1rho MAT files must start with "T1rho_S" and T2* MAT
%             files must start with "T2star_S".  Segmentation MAT file
%             names must contain "mroi".  See rd_m_dicom.m and
%             seg_m_rois.m.
%
%             3.  M-file exp_fun1.m, cmprt_ana4m.m and cmprt_plt4m.m
%             must be in the current directory or path.
%
%     22-Jul-2022 * Mack Gardner-Morse
%

%#######################################################################
%
% Setup Color Map
%
gmap = gray(128);       % Gray color map for not cartilage
jmap = jet(128);        % Jet color map for cartilage
cmap = [gmap; jmap];
%
% Bounding Box for Plot Zoom Window
%
% boxr = [150 345 185 320];    % T1rho bounding box - not used - see bbox
% boxs = [90 250 150 250];     % T2* bounding box - not used - see bbox
%
% Set Curvefit Optimization Parameters
%
opt = optimset('Display','off','TolFun',1e-8,'TolX',1e-8,'MaxIter', ...
               2e+3,'Algorithm','levenberg-marquardt','Jacobian', ...
               'on','UseParallel',true);
%
fun = @exp_fun1;        % Exponential function
%
% Initialize Parameters
%
% init = -1;              % Use weighted least squares for starting parameters
% init = 0;               % Use linear least squares for starting parameters
init = 1;               % Use fixed starting parameters
tr0 = 20;               % 21.45 ms - mean w/ 50 ms threshold
% tr0 = 65;               % Initial T1rho estimate in ms
% tr0 = 80;               % Initial T1rho estimate in ms
trmx = 100;             % Maximum valid T1rho result
trmn = 0;               % Minimum valid T1rho result
ts0 = 12;               % 12.38 ms - mean w/ 50 ms threshold
% ts0 = 35;               % Initial T2* estimate in ms
tsmx = 100;             % Maximum valid T2* result
tsmn = 0;               % Minimum valid T2* result
%
mxtr = 60;              % Maximum scale on T1rho plots
mxts = 55;              % Maximum scale on T2* plots
%
% Output Directory, Output Files and Output Labels
%
resdir = fullfile('Results');          % Results directory
%
ifirst = true;          % First write to file
xlsnam = 'mri_m_fit.xlsx';             % Results spreadsheet
xlsnam = fullfile(resdir,xlsnam);      % Include output directory
hdrs1 = {'Subject' 'Result' 'Leg' 'Comprt' 'AP'};
hdrs2 = {'Pixels' 'T1R/T2S' 'RSS' 'ValidPix' 'Mean' 'Min' 'Max' ...
         'SD' 'COV'};
%
psnam = fullfile(resdir,'mri_m_fit_'); % Start of PS file name
pstyp = '.ps';          % PS file type
%
% Get Subject Directories
%
spath = fullfile('..','PTOA');         % Path to series MAT files
sdirs = dir('..\PTOA\0*');
sdirs = {sdirs([sdirs.isdir]').name}'; % Subject directories
nsubj = size(sdirs,1);
%
o3 = ones(1,3);         % Column index for variable "id"
o5 = ones(1,5);         % Column index for variable "ids"
%
% Initialize Results Variables
%
% Indices key:
%   Index 1 - Subject
%   Index 2 - Leg - 1 = left and 2 = right
%   Index 3 - Compartment - 1 = lateral and 2 = medial
%   Index 4 - AP - 1 = anterior and 2 = posterior
%
t1r_res = zeros(nsubj,2,2,2);
t1r_npx = zeros(nsubj,2,2,2);
t1r_rss = zeros(nsubj,2,2,2);
%
t1r_respx = cell(nsubj,2,2,2);
t1r_rsspx = cell(nsubj,2,2,2);
t1r_nps = cell(nsubj,2,2,2);
%
t2s_res = zeros(nsubj,2,2,2);
t2s_npx = zeros(nsubj,2,2,2);
t2s_rss = zeros(nsubj,2,2,2);
%
t2s_respx = cell(nsubj,2,2,2);
t2s_rsspx = cell(nsubj,2,2,2);
t2s_nps = cell(nsubj,2,2,2);
%
% Loop through Subjects
%
for ks = 1:nsubj
%
% Get Subject Directory, Name and Number
%
   sdir = sdirs{ks};    % Current subject directory
   subj = eval(sdir);   % Subject number
%
   spdir = fullfile(spath,sdir);       % Include series path
%
   psnams = [psnam sdir];              % Add subject to PS file name
%
% Get T1rho MAT Files in Directory
%
   d = dir(fullfile(spdir,'T1rho_S*.mat'));
   rhonams = {d.name}';
   idr = contains(rhonams,'roi','IgnoreCase',true);   % ROI files
   rhonams = rhonams(~idr);
   idr = contains(rhonams,'chk','IgnoreCase',true);   % Check files
   rhonams = rhonams(~idr);
   nrho = size(rhonams,1);
%
   d = dir(fullfile(sdir,'T1rho_S*mrois.mat'));       % ROI files
   roinams = {d.name}';
   nroi = size(roinams,1);
%
   if nrho~=nroi
     error([' *** ERROR in mri_m_fit:  Number of T1rho MAT', ...
            ' files not equal to number of ROI MAT files!']);
   end
%
% T1rho Identifier
%
   ires = 0;            % ires = 0 - T1rho, ires = 1 - T2*
   idt = 1;             % Spin lock/echo time for plots - 1 = 0 ms spin lock time
%
   psnamr = [psnams '_T1R_'];          % Add result type to PS file name
%
% Loop through T1rho MAT Files
%
   for km = 1:nroi
%
% Load Data
%
      roinam = roinams{km};
      load(fullfile(sdir,roinam),'bbox','cmprt','maskmal','maskmam', ...
           'maskmpl','maskmpm','rsll','rslm');
%
      idm = contains(rhonams,roinam(1:end-10));  % Get matching file
      if ~any(idm)
        error([' *** ERROR in mri_m_fit:  Matching T1rho MAT', ...
               ' file not found for ROI MAT file:  ' roinam '!']);
      end
      rhonam = rhonams{idm};
      load(fullfile(spdir,rhonam),'iszs','nslt','scmx','sns', ...
           'snt','splt','st','v');
%
      npix = prod(iszs);     % Number of pixels in an image
      fs = ['S' snt];        % Series number prefaced with a 'S'
%
% Parse Series Text for Leg
%
      if strcmpi(st(1),'L')
        leg = 'L';
        ileg = 0;       % Coding for leg
        ltxt = 'Left';
      else
        leg = 'R';
        ileg = 1;
        ltxt = 'Right';
      end
%
% Add Leg to PS File Name
%
      psnamf = [psnamr leg pstyp];     % Add leg to PS file name
%
% Get Compartment Identifier for Slices
%
      rsls = {rsll'; rslm'};           % Lateral - row 1, medial - row 2
      nrsls = [size(rsll,1); size(rslm,1)];
%
% Combine Masks into a Cell Array
%
      maskmapl = {maskmal; maskmpl};   % Combine anterior and posterior masks
      maskmapm = {maskmam; maskmpm};   % Combine anterior and posterior masks
      mask = {maskmapl; maskmapm};     % Combine compartment masks
%
% Do Compartmental Analysis
%
      [tc,~,rss,npx,id,tcp,ampp,rssp,nps] = cmprt_ana4m(v,mask, ...
                                 rsls,nrsls,splt,nslt,fun,init,tr0,opt);
      na = size(tc,1);                 % Number of results
%
% Save Results
%
% Indices key:
%   Index 1 - Subject
%   Index 2 - Leg - 1 = left and 2 = right
%   Index 3 - Compartment - 1 = lateral and 2 = medial
%   Index 4 - AP - 1 = anterior and 2 = posterior
%
      for ka = 1:na
         t1r_res(ks,ileg+1,id(ka,1)+1,id(ka,2)+1) = tc(ka);
         t1r_npx(ks,ileg+1,id(ka,1)+1,id(ka,2)+1) = npx(ka);
         t1r_rss(ks,ileg+1,id(ka,1)+1,id(ka,2)+1) = rss(ka);
         t1r_respx{ks,ileg+1,id(ka,1)+1,id(ka,2)+1} = tcp{ka};
         t1r_rsspx{ks,ileg+1,id(ka,1)+1,id(ka,2)+1} = rssp{ka};
         t1r_nps{ks,ileg+1,id(ka,1)+1,id(ka,2)+1} = nps{ka};
      end
%
% Plot Results
%
      sid = ['Subject ' sdir ', ' ltxt ' Leg, T1\rho'];
      cmprt_plt4m(v,mask,rsls,nrsls,idt,tcp,nps,mxtr,cmap,bbox,sid, ...
                  psnamf);
%
% Get Statistics on Pixel Results
%
      npxv = zeros(na,1);              % Number of valid results
      tcpm = zeros(na,1);              % Mean
      tcpmn = zeros(na,1);             % Minimum
      tcpmx = zeros(na,1);             % Maximum
      tcpsd = zeros(na,1);             % SD
%
      for ka = 1:na
         idv = tcp{ka}>=trmn&tcp{ka}<=trmx;
         npxv(ka) = sum(idv);          % Number of valid results
         if npxv(ka)>0
           tcpv = tcp{ka}(idv);        % Valid T1rho values
           tcpm(ka) = mean(tcpv);      % Mean
           tcpmn(ka) = min(tcpv);      % Minimum
           tcpmx(ka) = max(tcpv);      % Maximum
           tcpsd(ka) = std(tcpv);      % SD
         end
      end
%
      tcpcov = 100*tcpsd./tcpm;        % Coefficient of variation
      tcpcov(isnan(tcpcov)) = 0;       % Catch any NaNs
%
% Combine Identifiers
%
      ids = [subj ires ileg];          % MAT file identifiers
      ids = repmat(ids,na,1);
      ids = [ids id];                  % All identifiers
%
% Create and Write Table of Results
%
      t1 = array2table(ids,'VariableNames',hdrs1);
      t2 = table(npx,tc,rss,npxv,tcpm,tcpmn,tcpmx,tcpsd,tcpcov, ...
                 'VariableNames',hdrs2);
      t = [t1 t2];
%
      if ifirst
        writetable(t,xlsnam,'WriteMode','replacefile');
        ifirst = false;
      else
        writetable(t,xlsnam,'WriteMode','append', ...
                   'WriteVariableNames',false);
      end
%
   end                  % End of km loop - T1rho MAT file loop
%
   close all;           % Close all plot windows
%
% Get T2* MAT Files in Directory
%
   d = dir(fullfile(spdir,'T2star_S*.mat'));
   starnams = {d.name}';
   idr = contains(starnams,'roi','IgnoreCase',true);   % ROI files
   starnams = starnams(~idr);
   idr = contains(starnams,'chk','IgnoreCase',true);   % Check files
   starnams = starnams(~idr);
   nstar = size(starnams,1);
%
   d = dir(fullfile(sdir,'T2star_S*mrois.mat'));
   roinams = {d.name}';
   nroi = size(roinams,1);
%
   if nstar~=nroi
     error([' *** ERROR in mri_m_fit:  Number of T2* MAT', ...
            ' files not equal to number of ROI MAT files!']);
   end
%
% T2* Identifier
%
   ires = 1;            % ires = 0 - T1rho, ires = 1 - T2*
   idt = 3;             % Spin lock/echo time for plots - 3 = 5 ms echo time
%
   psnamr = [psnams '_T2S_'];          % Add result type to PS file name
%
% Loop through T2* MAT Files
%
   for km = 1:nroi
%
% Load Data
%
      roinam = roinams{km};
      load(fullfile(sdir,roinam),'bbox','cmprt','maskmal','maskmam', ...
           'maskmpl','maskmpm','rsll','rslm');
%
      idm = contains(starnams,roinam(1:end-10)); % Get matching file
      if ~any(idm)
        error([' *** ERROR in mri_m_fit:  Matching T2* MAT', ...
               ' file not found for ROI MAT file:  ' roinam '!']);
      end
      starnam = starnams{idm};
      load(fullfile(spdir,starnam),'etns','iszs','netn','scmx', ...
           'sns','snt','st','v');
%
      npix = prod(iszs);     % Number of pixels in an image
      fs = ['S' snt];        % Series number prefaced with a 'S'
%
% Parse Series Text for Leg
%
      if strcmpi(st(1),'L')
        leg = 'L';
        ileg = 0;       % Coding for leg
        ltxt = 'Left';
      else
        leg = 'R';
        ileg = 1;
        ltxt = 'Right';
      end
%
% Add Leg to PS File Name
%
      psnamf = [psnamr leg pstyp];     % Add leg to PS file name
%
% Get Compartment Identifier for Slices
%
      rsls = {rsll'; rslm'};           % Lateral - row 1, medial - row 2
      nrsls = [size(rsll,1); size(rslm,1)];
%
% Combine Masks into a Cell Array
%
      maskmapl = {maskmal; maskmpl};   % Combine anterior and posterior masks
      maskmapm = {maskmam; maskmpm};   % Combine anterior and posterior masks
      mask = {maskmapl; maskmapm};     % Combine compartment masks
%
% Do Compartmental Analysis
%
      [tc,~,rss,npx,id,tcp,ampp,rssp,nps] = cmprt_ana4m(v,mask, ...
                                 rsls,nrsls,etns,netn,fun,init,ts0,opt);
      na = size(tc,1);              % Number of results
%
% Save Results
%
% Indices key:
%   Index 1 - Subject
%   Index 2 - Leg - 1 = left and 2 = right
%   Index 3 - Compartment - 1 = lateral and 2 = medial
%   Index 4 - AP - 1 = anterior and 2 = posterior
%
      for ka = 1:na
         t2s_res(ks,ileg+1,id(ka,1)+1,id(ka,2)+1) = tc(ka);
         t2s_npx(ks,ileg+1,id(ka,1)+1,id(ka,2)+1) = npx(ka);
         t2s_rss(ks,ileg+1,id(ka,1)+1,id(ka,2)+1) = rss(ka);
         t2s_respx{ks,ileg+1,id(ka,1)+1,id(ka,2)+1} = tcp{ka};
         t2s_rsspx{ks,ileg+1,id(ka,1)+1,id(ka,2)+1} = rssp{ka};
         t2s_nps{ks,ileg+1,id(ka,1)+1,id(ka,2)+1} = nps{ka};
      end
%
% Plot Results
%
      sid = ['Subject ' sdir ', ' ltxt ' Leg, T2*'];
      cmprt_plt4m(v,mask,rsls,nrsls,idt,tcp,nps,mxts,cmap,bbox,sid, ...
                  psnamf);
%
% Get Statistics on Pixel Results
%
      npxv = zeros(na,1);              % Number of valid results
      tcpm = zeros(na,1);              % Mean
      tcpmn = zeros(na,1);             % Minimum
      tcpmx = zeros(na,1);             % Maximum
      tcpsd = zeros(na,1);             % SD
%
      for ka = 1:na
         idv = tcp{ka}>=tsmn&tcp{ka}<=tsmx;
         npxv(ka) = sum(idv);          % Number of valid results
         if npxv(ka)>0
           tcpv = tcp{ka}(idv);        % Valid T2* values
           tcpm(ka) = mean(tcpv);      % Mean
           tcpmn(ka) = min(tcpv);      % Minimum
           tcpmx(ka) = max(tcpv);      % Maximum
           tcpsd(ka) = std(tcpv);      % SD
         end
      end
%
      tcpcov = 100*tcpsd./tcpm;        % Coefficient of variation
      tcpcov(isnan(tcpcov)) = 0;       % Catch any NaNs
%
% Combine Identifiers
%
      ids = [subj ires ileg];          % MAT file identifiers
      ids = repmat(ids,na,1);
      ids = [ids id];                  % All identifiers
%
% Create and Write Table of Results
%
      t1 = array2table(ids,'VariableNames',hdrs1);
      t2 = table(npx,tc,rss,npxv,tcpm,tcpmn,tcpmx,tcpsd,tcpcov, ...
                 'VariableNames',hdrs2);
      t = [t1 t2];
%
      writetable(t,xlsnam,'WriteMode','append', ...
                 'WriteVariableNames',false);
%
   end                  % End of km loop - T2* MAT file loop
%
   close all;           % Close all plot windows
%
end                     % End of ks loop - subjects loop
%
% Save to MAT File
%
save(fullfile(resdir,'mri_m_fit.mat'),'t1r_res','t1r_npx','t1r_rss', ...
     't1r_respx','t1r_rsspx','t1r_nps','t2s_res','t2s_npx', ...
     't2s_rss','t2s_respx','t2s_rsspx','t2s_nps');
%
return