%#######################################################################
%
%       * Meniscus Erosion DAYS between Injury and MRI Program *
%
%          M-File which reads the PTOA days between injury and MRI data
%     from the MS-Excel file in the ..\PTOA\Results\ subdirectory:
%     "Copy of mri_fitps_stat_MK_19May23.xlsx".
%
%          The meniscus T1rho and T2* results are read and the
%     days are incorporated into the results.  The expanded results
%     are written to the MS-Excel spreadsheets:
%     mri_me_fit_no_stat.xlsx, mri_me_fit_ec_stat.xlsx, and
%     mri_me_fit_es_stat.xlsx in the Meniscus\Results\ subdirectory.
%
%     NOTES:  1.  Assumes Matlab is starting in the Meniscus directory.
%
%     09-Jun-2023 * Mack Gardner-Morse
%

%#######################################################################
%
% Read Days between Injury and MRI Data for All PTOA Subjects from
% MS-Excel Spreadsheet File
%
data = readtable(fullfile('..','PTOA','Results', ...
                 'Copy of mri_fitps_stat_MK_19May23.xlsx'));
subjn = data.Subject;
days = data.DaysBtwInj_MRI;
[subjn,idu] = unique(subjn);
days = days(idu);
%
% Get Results Spreadsheets
%
xlsnams = {'mri_me_fit_no_stat.xlsx'
           'mri_me_fit_ec_stat.xlsx'
           'mri_me_fit_es_stat.xlsx'};
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
   hdrs = res.Properties.VariableNames;
   hdrs = [hdrs {'DaysBtwInj_MRI'}];
%
% Expand Subject Demographics to All Trials
%
   subjs = res(:,1).Subject;
%
   idx = double(subjs==subjn');
   daysk = idx*days;
%
% Expand Table with Demographics
%
   res = addvars(res,daysk);
   res.Properties.VariableNames = hdrs;
%
% Write Table to MS-Excel Spreadsheet
%
   writetable(res,xlsnam,'WriteMode','replacefile');
%
end
%
return