function individual_decoding_accuracies(cfg)

if ~isfield(cfg, 'smoothing_window'); cfg.smoothing_window = 6;end

%% load data
filepath = fullfile(cfg.outputPath, ['PEP_WP4_EEG_decoding_accuracy', num2str(cfg.smoothing_window),'_TP_RA_reref_w.mat']);
load(filepath);

%% plot

% define parameters
analysis_title="individual decoding accuracies" ;
analysis_color=[104, 138, 189];
chance_level=1/2;

%figure('position',[1,1,1000,600], 'unit','centimeters');
figure;
set(gcf, 'color', [1 1 1]); % white background

% graphic parameters
set(0, 'DefaultTextFontSize', cfg.FontSize);
set(0, 'DefaultAxesFontSize', cfg.FontSize);
set(0, 'DefaultTextFontName', cfg.FontName);
set(0, 'DefaultAxesFontName', cfg.FontName);

dec_acc_min=min(min(res.dec_acc));
dec_acc_max=max(max(res.dec_acc));

% plot individual decoding results
for s=1:cfg.n
    subplot(ceil((cfg.n)/2),2,s)
    xline(0);

    hold on
    h=line([min(res.time),max(res.time)],[chance_level,chance_level]);
    set(h,'color','k');
    set(h,'linewidth',1);
    set(h,'linestyle','--');
    hold on;
    title(['subject ' num2str(res.order(s))]);
    p=plot(res.time,squeeze(res.dec_acc(s,:)));
    set(p,'linewidth',1.25);
    set(p,'Color',analysis_color./255);
    yline(max(res.dec_acc(s,:)), 'k--', ['max ', num2str(max(res.dec_acc(s,:)))])
    ylim([dec_acc_min, dec_acc_max+0.02])
end

subplot_title=sgtitle(char(analysis_title));
set(subplot_title,'FontWeight','bold');

% save
if cfg.saving
    output_dir=fullfile('..', 'Plots');
    filename = fullfile(output_dir, 'decoding_accuracy_individual_reref_w.jpg');
    saveas(gcf, filename);
end
end