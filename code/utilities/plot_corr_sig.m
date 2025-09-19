function plot_corr_sig(cfg, d, pval, save_name)

% evaluate input
if ~isfield(cfg, 'smoothing_window'); cfg.smoothing_window = 10; end
if ~isfield(cfg, 'ylim'); cfg.ylim = [-0.1, 0.1]; end

% create figure and layout
fig = figure;
fig.Position = [100 100 1400 500]; 
tiledlayout(1,2);

% prep time points and variables
timepoints = d.pairRep.included_time;
timeseries = d.pairRep.all_time;
x = timeseries(timepoints);
x_lim = ([min(x)-0.01, max(x)+0.01]);
ntp_before_stim = length(d.pairRep.included_time) - length(d.stats.testedTime);

for category = cfg.categories
    nexttile;
    hold on;
    category = char(category);
    legend_labels = cell(1,3);
    plot_counter = 1;
    g = gobjects(1, length(cfg.RDM_to_partial_out)*length(cfg.categories));

    for var = 1:length(cfg.RDM_to_partial_out)

        % smooth data
        y = d.resMat.partial_cor.(category)(var, :);
        y = smoothdata(y, 2, 'movmean', cfg.smoothing_window);

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
        g(plot_counter) = plot(x, y, 'color', clr, 'LineWidth', 3);

        % mark signifikant time points
        pos = 0.1 - var*0.01;
        plot(x(sigMat), repmat(pos, 1, sum(sigMat)), ...
            'color', clr, 'marker' ,'O', 'MarkerFaceColor', clr ,'MarkerSize', 5, 'LineStyle','none');

        legend_labels{plot_counter} = [strrep(cfg.RDM_to_partial_out{var}, '_', '-'), ' ', category];
        plot_counter = plot_counter+1;
    end


    xlim(x_lim);
    ylim(cfg.ylim);
    yline(0, '--')
    xline(0, '--')
    xlabel('time');
    legend([g(1:plot_counter-1)], legend_labels, 'Location', 'best');
    set(gca, 'box', 'off');

end % category

sgtitle('Time-resolved correlation ISC-RDM with predictors', 'FontSize', 18);
hold off;
save_plot(fig, ['corr_sig_', save_name], cfg.figPath)

% create figure
fig = figure;
fig.Position = [100 100 1000 500]; 
hold on;

% prep time points and variables
legend_labels = cell(1,3);
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
    g(plot_counter) = plot(x, y, 'color', clr, 'LineStyle', '-', 'LineWidth', 3);

    % mark signifikant time points
    pos = 0.1 - var*0.01;
    plot(x(sigMat), repmat(pos, 1, sum(sigMat)), ...
        'color', clr, 'marker' ,'O', 'MarkerFaceColor', clr ,'MarkerSize', 5, 'LineStyle','none');

    legend_labels{plot_counter} = [strrep(cfg.RDM_to_partial_out{var}, '_', '-'),' mean'];
    plot_counter = plot_counter + 1;
end

if cfg.partial_cor
    ylabel(['Partial correlation [r]', newline]);
else
    ylabel([cfg.correlation_type, ' correlation [r]', newline]);
end

set(gca, 'LineWidth', 1, 'FontName', cfg.FontName, 'FontSize', cfg.FontSize, 'FontWeight', 'bold');
ylim(cfg.ylim);
xlim(x_lim);
yline(0, '--');
xline(0, '--');
xlabel('time');
legend([g(1:plot_counter-1)], legend_labels, 'Location', 'best');
set(gca, 'box', 'off');

sgtitle('Time-resolved correlation ISC-RDM with predictors', 'FontSize', 18);
hold off;

save_plot(fig, ['mean_corr_sig', save_name], cfg.figPath)

end