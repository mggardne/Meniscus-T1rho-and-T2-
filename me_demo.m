%#######################################################################
%
%               * Meniscus Erosion DEMOgraphics Program *
%
%          M-File which reads the PTOA demographics data from the CSV 
%     file in the ..\PTOA\Results\ subdirectory:
%     PTOAStudyTheAdaptive-DataPull3May2023_DATA_2023-05-03_1359.csv.
%
%          The meniscus T1rho and T2* results are read and the
%     demographics are incorporated into the results.  The expanded
%     results are written to the MS-Excel spreadsheets:
%     mri_me_fit_no_stat.xlsx, mri_me_fit_ec_stat.xlsx, and
%     mri_me_fit_es_stat.xlsx in the Meniscus\Results\ subdirectory.
%
%     NOTES:  1.  Assumes Matlab is starting in the Meniscus directory.
%
%     05-Jun-2023 * Mack Gardner-Morse
%

%#######################################################################
%
% Read Demographics Data for All PTOA Subjects from CSV File
%
data = readtable(fullfile('..','PTOA','Results', ...
     'PTOAStudyTheAdaptive-DataPull3May2023_DATA_2023-05-03_1359.csv'));
subjID = data.ptoa_subj;
subjID = strrep(subjID,'PTOA','');
subjID = strrep(subjID,'-1','');
%
% Get PTOA Subject Demographics
%
subjn = str2double(subjID);
sex = 2-(data(:,2).dem01);
age = data(:,3).dem02;
lmenis = table2array(data(:,4));
bmi = data(:,5).bmi2;
%
% Get Results Spreadsheets
%
xlsnams = {'mri_me_fit_no.xlsx'
           'mri_me_fit_ec.xlsx'
           'mri_me_fit_es.xlsx'};
xlsnams = fullfile('Results',xlsnams);
nr = size(xlsnams,1);
%
% Read Results Spreadsheet
%
for k = 1:nr            % Loop through results spreadsheets
%
   xlsnam = xlsnams{k};
   res = readtable(xlsnam);
%
% Get Variable Names (Column Headers)
%
   hdrs = res.Properties.VariableDescriptions;
   hdrs = [hdrs(1) {'Sex' 'Age' 'BMI' 'LatMeniscus'} hdrs(2:end)];
%
% Expand Subject Demographics to All Trials
%
   subjs = res(:,1).Subject;
%
   idx = double(subjs==subjn');
   sexk = idx*sex;
   agek = idx*age;
   bmik = idx*bmi;
   lmenisk = idx*lmenis;
%
% Expand Table with Demographics
%
   res = addvars(res,sexk,agek,bmik,lmenisk,'After','Subject');
   res.Properties.VariableNames = hdrs;
%
% Write Table to MS-Excel Spreadsheet
%
   idot = strfind(xlsnam,'.');
   idot = idot(end);
   xlsnamo = [xlsnam(1:idot-1) '_stat' xlsnam(idot:end)];
   writetable(res,xlsnamo,'WriteMode','replacefile');
%
end
%
return