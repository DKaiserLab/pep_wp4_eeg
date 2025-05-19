%% Housekeeping
clear 
clc
close all

slist=[6];

for s=slist

    %% load data
    filepath = fullfile(pwd, '..', '..', 'derivatives', ['PEP_WP4_EEG_decoding_accuracy', num2str(s), '_TP_RA.mat']);
    load(filepath);


    %% plot decoding results

    dec_acc_mean=squeeze(mean(res.dec_acc,1));

    figure();
    decoding_plot=plot(res.time, dec_acc_mean(:));
    chance=1/2;
    grid("on")
    ylim([0.48,0.55])
    yline(max(dec_acc_mean), 'k--', ['max ', num2str(max(dec_acc_mean))])
    yline(chance,'k--');
    xline(0,'k');
    title(['eeg decoding']);

    %%

    output_dir=fullfile('..','..','..', 'Plots');
    filename = fullfile(output_dir, ['decoding_accuracy_without110', num2str(s), '.jpg']);

    % Speichern
    saveas(gcf, filename);  % oder: print(gcf, filename, '-djpeg', '-r300');

end