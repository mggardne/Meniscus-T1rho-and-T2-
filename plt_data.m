%
% Run Matlab program "ForceClamp_UO_Aurora_GM.m" first to load the data.
%
t = data(:,1);
id2 = 14001:size(t,1);
%
txt = {'L_{in}'; 'L_{out}'; 'F_{in}'; 'F_{out}'};
%
% Plot Force Data
%
figure;
orient landscape;
%
plot(t(id2)/1000-20,data(id2,4:5));
xlabel('Time (s)','FontSize',12,'FontWeight','bold');
ylabel('Force (mN)','FontSize',12,'FontWeight','bold');
legend(txt(3:4),'Location','West','FontSize',12,'FontWeight','bold');
print -dpsc2 -r600 -fillpage forceclamp.ps
%
% Plot Length Data
%
figure;
orient landscape;
%
plot(t(id2)/1000-20,data(id2,2:3));
xlabel('Time (s)','FontSize',12,'FontWeight','bold');
ylabel('Length (mm)','FontSize',12,'FontWeight','bold');
legend(txt(1:2),'Location','West','FontSize',12,'FontWeight','bold');
print -dpsc2 -r600 -fillpage -append forceclamp.ps
%
return