function plot_corr_sig(cfg, d, pval, save_name)

% evaluate input
if ~isfield(cfg, 'smoothing_window'); cfg.smoothing_window = 10; end
if ~isfield(cfg, 'ylim'); cfg.ylim = [-0.1, 0.1]; end
if ~isfield(cfg, 'xlim'); cfg.xlim = [-200, 500]; end

% create figure and layout
fig = figure;
set(gcf, 'Units','centimeters', 'Position',[2 2 80 20]);
%fig.Position = [100 100 1400 500];
tiledlayout(1,2);

% graphic parameters
set(0, 'DefaultTextFontSize', cfg.FontSize);
set(0, 'DefaultAxesFontSize', cfg.FontSize);
set(0, 'DefaultTextFontName', cfg.FontName);
set(0, 'DefaultAxesFontName', cfg.FontName);
set(gcf, 'color', [1 1 1]);

% prep time points and variables
timepoints = d.pairRep.included_time;
timeseries = d.pairRep.all_time;
x = timeseries(timepoints)*1000;
ntp_before_stim = length(d.pairRep.included_time) - length(d.stats.ISC_RSA.testedTime);

for icat = 1:length(cfg.categories)
    nexttile;
    hold on;
    category = cfg.categories{icat};
    plot_counter = 1;
    g = gobjects(1, length(cfg.RDM_to_partial_out));

    for var = 1:length(cfg.RDM_to_partial_out)

        % smooth data
        y = d.resMat.partial_cor.(category)(var, :);
        y = smoothdata(y, 2, 'movmean', cfg.smoothing_window+5);

        sigMat = false(1,length(timepoints));
        % compute significant time points
        if isfield(pval, category)
            sigMat(ntp_before_stim+1:end) = pval.(category)(var,:) <= 0.05;
        elseif isfield(pval, 'FDR')
            sigMat(ntp_before_stim+1:end) = pval.FDR.(category)(var,:);
        elseif isfield(pval, 'cluster')
           if ~isempty(pval.cluster.(category){var}.obs_cluster_stats)
                for p=1:length(pval.cluster.(category){var}.pvals_cluster)
                    if pval.cluster.(category){var}.pvals_cluster(p) <= cfg.alpha
                        idx = pval.cluster.(category){var}.obs_clusters.PixelIdxList{p} + ntp_before_stim ;
                        sigMat(idx)=1;
                    end
                end
            end
        else
            error('Not defined.')
        end


        % get color
        if startsWith(cfg.RDM_to_partial_out{var}, 'typical')
            clr = [1, 0, 1];
        elseif startsWith(cfg.RDM_to_partial_out{var}, 'control')
            clr = [.7, .7, .7];
        elseif startsWith(cfg.RDM_to_partial_out{var}, 'photos')
            clr = [.4, .9, 1];
        end

        % plot line
        g(plot_counter) = plot(x, y, 'color', clr, 'LineWidth', 5, 'DisplayName', [regexprep(cfg.RDM_to_partial_out{var}, "_.*", ""), ' ', category]);

        % mark signifikant time points
        pos = 0.1 - var*0.01;
        plot(x(sigMat), repmat(pos, 1, sum(sigMat)), ...
            'color', clr, 'marker' ,'O', 'MarkerFaceColor', clr ,'MarkerSize', 7, 'LineStyle','none');

        plot_counter = plot_counter+1;
    end


    xlim(cfg.xlim);
    ylim(cfg.ylim);
    yline(0, '--', 'LineWidth', 2);
    xline(0, '--', 'LineWidth', 2);
    legend(g, 'Location', 'southeast');
    legend('boxoff');
    xlabel('Time (ms)');
    if icat ==1; ylabel('Partial correlation'); end
    set(gca, 'box', 'off');
    set(gca, 'FontWeight', 'bold', 'LineWidth', 2);

end % category

sgtitle('Time-resolved IS-RSA', 'FontSize', cfg.FontSize+10, 'FontWeight', 'bold');
hold off;
save_plot(fig, ['corr_sig_', save_name], cfg.figPath)

%% plot average across categories
% create figure
fig = figure;
set(gcf, 'Units','centimeters', 'Position',[2 2 40 25]);
%fig.Position = [100 100 1000 500]; 
set(gcf, 'color', [1 1 1]);
hold on;

% prep time points and variables
plot_counter = 1;
g = gobjects(1, length(cfg.RDM_to_partial_out));

for var = 1:numel(cfg.RDM_to_partial_out)
    % smooth data
    % plot mean
    y1 = d.resMat.partial_cor.bathroom(var, :);
    y2 = d.resMat.partial_cor.kitchen(var, :);
    y = mean([y1; y2]);
    y = smoothdata(y, 2, 'movmean', cfg.smoothing_window);

    % compute signifikant time points
    sigMat = false(1, length(timepoints));
    if isfield(pval, category)
        sigMat(ntp_before_stim+1:end) = pval.all(var,:) <= 0.05;
    elseif isfield(pval, 'FDR')
        sigMat(ntp_before_stim+1:end) = pval.FDR.all(var,:);
    elseif isfield(pval, 'cluster')
        if ~isempty(pval.cluster.all{var}.obs_cluster_stats)
            for p=1:length(pval.cluster.all{var}.pvals_cluster)
                if pval.cluster.all{var}.pvals_cluster(p) <= cfg.alpha
                    idx = pval.cluster.all{var}.obs_clusters.PixelIdxList{p} + ntp_before_stim ;
                    sigMat(idx)=1;
                end
            end
        end
    else
        error('Not defined.')
    end

    % get color
    if startsWith(cfg.RDM_to_partial_out{var}, 'typical')
        clr = [1, 0, 1];
    elseif startsWith(cfg.RDM_to_partial_out{var}, 'control')
        clr = [.7, .7, .7];
    elseif startsWith(cfg.RDM_to_partial_out{var}, 'photos')
        clr = [.4, .9, 1];
    end

    % plot mean
    g(plot_counter) = plot(x, y, 'color', clr, 'LineStyle', '-', 'LineWidth', 5, ...
        'DisplayName', regexprep(cfg.RDM_to_partial_out{var}, "_.*", ""));

    % mark signifikant time points
    pos = 0.1 - var*0.01;
    plot(x(sigMat), repmat(pos, 1, sum(sigMat)), ...
        'color', clr, 'marker' ,'O', 'MarkerFaceColor', clr ,'MarkerSize', 7, 'LineStyle','none');

    plot_counter = plot_counter + 1;
end

if cfg.partial_cor
    ylabel('Partial correlation');
else
    ylabel([cfg.correlation_type, ' correlation [r]', newline]);
end


ylim(cfg.ylim);
xlim(cfg.xlim);
yline(0, '--', 'LineWidth', 2);
xline(0, '--', 'LineWidth', 2);
xlabel('Time (ms)');
legend(g, 'Location', 'northeast');
legend('boxoff');
set(gca, 'box', 'off');
%sgtitle('Time-resolved IS-RSA', 'FontSize', cfg.FontSize + 10, 'FontWeight', 'bold');
set(gca, 'FontWeight', 'bold', 'LineWidth', 2);
hold off;

save_plot(fig, ['mean_corr_sig', save_name], cfg.figPath)

end