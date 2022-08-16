function [tc,amp,rss,npx,id,tcpx,amppx,rsspx,nps] = cmprt_ana4m(v, ...
                             mask,rsls,nrsls,time,ntime,fun,init,t0,opt)
%CMPRT_ANA4M  Given a four-dimensional matrix of T1/T2 intensities from
%          a MRI image volume, a mask, and a vector of spin lock/echo
%          times, calculates the T1rho/T2* map for the volume.
%
%          TC = CMPRT_ANA4M(V,MASK,RSLS,NRSLS,TIME,NTIME,FUN,INIT) Given
%          a four-dimensional matrix of T1/T2 intensities from a
%          MRI image volume, V, where the first two dimensions are an
%          image, the third dimension are the slices and the fourth
%          dimension are the spin lock/echo times, three dimensional
%          logical masks with the first dimension being the image, the
%          second dimension being the superficial layer in the first
%          column and deep layer in the second column and the third
%          dimension being slices in a cell array of masks with the
%          first index to lateral and medial and the second index to the
%          anterior and posterior, MASK, a cell array with the slices
%          within each compartment, RSLS, the number of slices in each
%          compartment, NRSLS, vector of spin lock/ echo times, TIME,
%          the number of spin lock/echo times, NTIME, a function handle,
%          FUN, to evaluate the monoexponential and its derivatives, an
%          initialization flag, INIT, to determine the initial starting
%          parameters, calculates T1rho/T2* for the two compartments
%          (lateral and medial) and two positions (anterior and
%          posterior) in array TC.
%
%          TC = CMPRT_ANA4M(V,MASK,RSLS,NRSLS,TIME,NTIME,FUN,INIT,T0) If
%          INIT is greater than zero, uses T0 as the initial default
%          time constant for the nonlinear exponential fit.  The
%          default value is 50.  The mean value of the maximum
%          intensities are used as the initial amplitude.  If INIT is
%          equal to zero, the exponential equation is linearized using
%          logarithms and least squares are used to solve for the
%          initial parameters for the nonlinear least squares.  If INIT
%          is less than zero, weighted least squares using the
%          intensities is used to solve for the starting parameters.
%
%          TC = CMPRT_ANA4M(V,MASK,RSLS,NRSLS,TIME,NTIME,FUN,INIT,T0,
%          OPT) a structure, OPT, with options for the curvefit solver.
%          See Matlab command "optimset" for a list of parameters.
%
%          [TC,AMP,RSS,NPX,ID] = CMPRT_ANA4M(...) returns the amplitudes
%          of the exponential fits, AMP, the sum of squared residuals,
%          RSS, the number of pixels in the curvefits, NPX, and a two
%          column array with the first column being the compartment (0 -
%          lateral and 1 - medial), and the second column being position
%          (0 - anterior and 1 - posterior) that identifies the
%          analyzes, ID. 
%
%          [TC,AMP,RSS,NPX,ID,TCPX,AMPPX,RSSPX,NPS] = CMPRT_ANA4M(...)
%          returns the T1rho/T2* for each pixel in the two compartments
%          (lateral and medial) and two positions (anterior and
%          posterior) in a cell array, TCPX, the exponential amplitudes
%          for each pixel in a cell array, AMPPX, the sum of squared
%          residuals for each pixel, RSSPX, and the number of pixels in
%          each slice within the regions in a cell array, NPS.
%
%          NOTES:  1.  The T1rho/T2* map is calculated as a mono-
%                  exponential as a function of spin lock/echo times.
%
%                  2.  The images volume data and spin lock/echo times
%                  are assumed to be in ascending order with increasing
%                  index into the fourth dimension of the image data V.
%
%                  3.  The function handle, FUN, is to the M-file
%                  exp_fun1.m which must be in the current directory or
%                  path.
%
%                  4.  T0 must be between 1 and 100.
%
%                  5.  The Matlab Optimization toolbox is required.
%                  The nonlinear least squares is performed by the
%                  Matlab function lsqcurvefit in the optimization
%                  toolbox.
%
%                  6.  The nonlinear least squares uses the Levenberg-
%                  Marquardt algorithm in the Matlab function
%                  lsqcurvefit.
%
%                  7.  The Matlab Parallel Computing toolbox is 
%                  required.  The Matlab parallel construct parfor is
%                  used to calculate the T1rho/T2* values in parallel.
%                  Use Matlab command parpool to control the number of 
%                  workers.
%
%                  8.  Based on the M-files cmprt_ana.m and
%                  cmprt_ana4.m, but modified by removing the combining
%                  of deep and superficial cartilage layers.  Adapted
%                  for use with calculating T1rho/T2* for the knee
%                  meniscus.
%
%          10-Aug-2022 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if nargin<8
  error([' *** ERROR in cmprt_ana4:  Eight input variables are', ...
         ' required!']);
end
%
if nargin<9
  t0 = 50;              % Initial estimate of time constant
end
%
% Check Initial Time Constant
%
if isempty(t0)||t0<1||t0>100
  t0 = 50;
end
%
% Default Optimization Parameters
%
if nargin<10
  opt = optimset('Display','off','TolFun',1e-8,'TolX',1e-8, ...
                 'MaxIter',2e+3,'Algorithm','levenberg-marquardt', ...
                 'Jacobian','on','UseParallel',true);
end
%
% Initialize Variables
%
% Columns in variable "id":
%   1 - Compartments:  0 - Lateral and 1 - Medial
%   2 - Positions:  0 - Anterior and 1 - Posterior
%
tc = zeros(4,1);        % Time constant - T1rho/T2*
amp = zeros(4,1);       % Amplitude of exponential
rss = zeros(4,1);       % Sum of squared residuals
npx = zeros(4,1);       % Number of pixels in regional fit
id = zeros(4,2);        % Columns: 1 - Compartment and 2 - Bone
%
nps = cell(4,1);        % Number of pixels in each slice
tcpx = cell(4,1);       % Time constant - T1rho/T2* for each pixel
amppx = cell(4,1);      % Amplitude of exponential for each pixel
rsspx = cell(4,1);      % Sum of squared residuals for each pixel
%
% Loop through Compartments
%
for kr = 1:2
%
   ikr = kr-1;          % Compartment code (0 - lateral, 1 - medial)
   mskr = mask{kr};     % Mask for this compartment
   rsl = rsls{kr};      % Slices for this compartment
   nrsl = nrsls(kr);    % Number of slices in this compartment
%
% Loop through Bones
%
   for kb = 1:2
%
      ikb = kb-1;       % Position code identifier (0 - anterior, 1 - posterior)
      msk = mskr{kb};   % Mask for this bone
%       msk = squeeze(msk(:,1,:)|msk(:,2,:)); % Combine cartilage layers
      rimgs = cell(nrsl,1);            % # of slices
%
% Loop through Slices
%
      np = zeros(nrsl,1);
%
      for ks = 1:nrsl
%
         slk = rsl(ks);                % Slice number
%
         rimgt = cell(ntime,1);        % # of spin lock/echo times
%
% Loop through Spin Lock/Echo Times
%
         for kt = 1:ntime
%
            rimg = squeeze(v(:,:,slk,kt));  % T1/T2 data for slice and spin lock/echo time
            rimgt{kt} = rimg(msk(:,ks))';   % Apply mask
%
         end            % End of kt loop - Spin lock/echo times
%
         rimgs{ks} = cat(1,rimgt{:});  % Combine spin lock/echo times
         np(ks) = sum(msk(:,ks));
%
      end               % End of ks loop - slices loop
%
      clear rimg rimgt;
%
% Calculate T1rho/T2*
%
      idx = 2*kr+kb-2;  % Index to variables
      nps{idx} = np;    % Number of pixels in slice
%
      rimgl = cat(2,rimgs{:});
      npx(idx) = size(rimgl,2);
      r = rimgl(:);
%
      amp0 = mean(rimgl(1,:));
      rp0c = [amp0; t0];               % Constant initial parameters
%
      tl = repmat(time,npx(idx),1);
%
      if init==0     % Use linear least squares for initial constants
        xdat = [tl ones(ntime*npx(idx),1)];
        lr = log(r);
        rpl = xdat\lr;
        rpl = abs(rpl);   % Get magnitude of complex numbers
        t0l = -1/rpl(1);
        if t0l>1&&t0l<100
          rp0 = [exp(rpl(2)); t0l];
        else
          rp0 = rp0c;
        end
      elseif init<0  % Use weighted least squares for initial constants
        xdat = [tl.*r r];
        lr = log(r);
        rpl = xdat\(lr.*r);
        rpl = abs(rpl);   % Get magnitude of complex numbers
        t0l = -1/rpl(1);
        if t0l>1&&t0l<100
          rp0 = [exp(rpl(2)); t0l];
        else
          rp0 = rp0c;
        end
      else           % Use constant initial parameters
        rp0 = rp0c;
      end
%
% Nonlinear Least Squares Exponential Fit to Get T1rho/T2* Value
%
      [rp,err] = lsqcurvefit(fun,rp0,tl,r,[],[],opt);
%
      tc(idx) = rp(2);                 % Time constant - T1rho/T2* 
      amp(idx) = rp(1);                % Amplitude of exponential
      rss(idx) = err;                  % Residual sum of squares
      id(idx,:) = [ikr ikb];           % Compartment, bone and layer IDs
%
% Calculate T1rho/T2* for Each Pixel in Each Layer
%
      if nargout>5
        tcp = zeros(npx(idx),1);
        ampp = tcp;
        rssp = tcp;
%
% Get Initial Parameters
%
        if init==0      % Use linear least squares for initial constants
          xdat = [time ones(ntime,1)];
          rpa = xdat\log(rimgl);
          rpa = abs(rpa);              % Get magnitude of complex numbers
          t0a = -1./rpa(1,:);
          ampa = exp(rpa(2,:));
          rpa = [ampa; t0a];
          idv = t0a<1|t0a>100;
          if any(idv)
            rpa(:,idv) = repmat(rp0c,1,sum(idv));
          end
        elseif init<0   % Use weighted least squares for initial constants
          rpa = repmat(rp0c,1,npx(idx));
          for p = 1:npx(idx)
             rimgp = rimgl(:,p);
             xdat = [time.*rimgp rimgp];
             lr = log(rimgp);
             rpl = xdat\(lr.*rimgp);
             rpl = abs(rpl);           % Get magnitude of complex numbers
             t0l = -1/rpl(1);
             if t0l>=1&&t0l<=100
               rpa(:,p) = [exp(rpl(2)); t0l];
             end
          end
        else            % Use constant initial parameters
          rpa = repmat(rp0c,1,npx(idx));
        end
%
% Nonlinear Least Squares Exponential Fit to Get T1rho/T2* Pixel Values
%
        parfor p = 1:npx(idx)
%      
           rimgp = rimgl(:,p);
%      
           [rp,err] = lsqcurvefit(fun,rpa(:,p),time,rimgp,[],[],opt);
           tcp(p) = rp(2);
           ampp(p) = rp(1);
           rssp(p) = err;
%
        end
%
% Save Results to Outputs
%
        tcpx{idx} = tcp;
        amppx{idx} = ampp;
        rsspx{idx} = rssp;
%
      end
%
   end                  % End of kb loop - bones loop
end                     % End of kr loop - compartments loop
%
return