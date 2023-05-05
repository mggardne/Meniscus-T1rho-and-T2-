%
% Get Subject Directory
%
sdir = '012-JH';        % Subject directory
%
% Load Masks
%
load(fullfile(sdir,'T1rho_S1401_mrois.mat'),'masks');
%
% Get Slice Mask as an Image
%
msk = masks(:,5);       % Slice 5
msk = reshape(msk,512,512);
%
% Get Different Structured Elements for Eroding the Image by One Pixel
% (Increase 3 to 5 to erode by two (2) pixels.)
%
se = strel('line',3,0);                % Horizontal line
se90 = strel('line',3,90);             % Vertical line
sec = strel(logical([0 1 0; 1 1 1; 0 1 0]));     % Crossed lines
ses = strel('square',3);
%
% Erode Images
%
emsk = imerode(msk,se);                % Horizontal line
emsk90 = imerode(msk,se90);            % Vertical line
emskc = imerode(msk,sec);              % Crossed lines
emsks = imerode(msk,ses);              % Square
%
% Show Differences
%
bbox = fnd_bbox(masks,[512,512]);
%
figure;
orient tall;
subplot(4,1,1);
imshowpair(msk,emsk);
axis(bbox);
title('Horizontal Line Erosion','FontSize',16,'FontWeight','bold');
%
subplot(4,1,2);
imshowpair(msk,emsk90);
axis(bbox);
title('Vertical Line Erosion','FontSize',16,'FontWeight','bold');
%
subplot(4,1,3);
imshowpair(msk,emskc);
axis(bbox);
title('Cross Erosion','FontSize',16,'FontWeight','bold');
%
subplot(4,1,4);
imshowpair(msk,emsks);
axis(bbox);
title('Square Erosion','FontSize',16,'FontWeight','bold');
%
print -dpdf -r600 img_erode_diff.pdf;
%
return