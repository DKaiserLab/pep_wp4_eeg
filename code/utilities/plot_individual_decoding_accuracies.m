%% Housekeeping
clear 
clc
close all

%% load data
filepath = fullfile(pwd, '..', '..', 'derivatives', 'PEP_WP4_EEG_decoding_accuracy6_TP_RA.mat');
load(filepath);

%% 
ss=[117 104 103 106 109 111 113 110 122 112];
an=1;
sn=size(res.dec_acc,1);
analysis_title="individual decoding accuracies" ;
analysis_color=[104, 138, 189];
chance_level=1/2;

hf=figure('position',[1,1,1000,600], 'unit','centimeters');
dec_acc_min=min(min(res.dec_acc));
dec_acc_max=max(max(res.dec_acc));
   
for s=1:sn
subplot(ceil((sn)/2),2,s)
xline(0); 

hold on
h=line([min(res.time),max(res.time)],[chance_level,chance_level]);
set(h,'color','k');
set(h,'linewidth',1);
set(h,'linestyle','--');
hold on
subj_id=title(['subject ' num2str(ss(s))]);
p=plot(res.time,squeeze(res.dec_acc(s,:)));
set(p,'linewidth',1.25);
set(p,'Color',analysis_color./255);
yline(max(res.dec_acc(s,:)), 'k--', ['max ', num2str(max(res.dec_acc(s,:)))])
ylim([dec_acc_min, dec_acc_max+0.02])
end
subplot_title=sgtitle(char(analysis_title));
set(subplot_title,'FontWeight','bold');

output_dir=fullfile('..','..', 'Plots');
filename = fullfile(output_dir, 'decoding_accuracy_individual.jpg');
saveas(gcf, filename);
