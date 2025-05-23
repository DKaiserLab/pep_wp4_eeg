%% Housekeeping
clear 
clc
close all

%% define parameter
slist=[6]; % smoothing window

%% plot decoding accuracy
for s=slist

    % load data
    filepath = fullfile(pwd, '..', '..', 'derivatives', ['PEP_WP4_EEG_decoding_accuracy', num2str(s), '_TP_RA.mat']);
    load(filepath);


    % plot decoding results

    dec_acc_mean=squeeze(mean(res.dec_acc,1));

    figure();
    decoding_plot=plot(res.time, dec_acc_mean(:));
    chance=1/2;
    grid("on")
    ylim([0.48,0.55])
    yline(max(dec_acc_mean), 'k--', ['max ', num2str(max(dec_acc_mean))])
    yline(chance,'k--');
    xline(0,'k');
    title(['decoding accuracy']);

    % save plot
    output_dir=fullfile('..','..', 'Plots');
    filename = fullfile(output_dir, ['decoding_accuracy', num2str(s), '.jpg']);

    saveas(gcf, filename); 

end