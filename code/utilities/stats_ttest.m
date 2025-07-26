function stats_ttest(cfg, d)

if ~isfield(cfg, 'chance_level'); cfg.chance_level = 0.5; end
if ~isfield(cfg, 'plotting'); cfg.plotting = true; end

% Initialize data storage
n_timepoints = size(d.pairRep.included_time, 2);
p_val = struct();
t_stat = struct();
p_val_fdr = struct();

for i = 1:length(cfg.categories)
    cat = cfg.categories{i};
    
    % Initialize data storage
    p_vals_cat = zeros(1, n_timepoints);
    t_vals_cat = zeros(1, n_timepoints);

    % compute t-test for each time point
    for t = 1:n_timepoints
        [~, p, ~, stats] = ttest(d.meanAcc.(cat).pairRep(t,:), cfg.chance_level, 'tail', 'right');
        p_vals_cat(t) = p;
        t_vals_cat(t) = stats.tstat;
    end

    % FDR
    p_vals_cat_fdr = fdr_bh(p_vals_cat);

    % save
    p_val.(cat) = p_vals_cat;
    t_stat.(cat) = t_vals_cat;
    p_val_fdr.(cat) = p_vals_cat_fdr;

end

% stats.p_val = p_val;
% stats.t_stat = t_stat;
% stats.p_val_fdr = p_val_fdr;

if cfg.plotting

    fig = figure;
    hold on;

    % graphic parameters
    set(0, 'DefaultTextFontSize', cfg.FontSize);
    set(0, 'DefaultAxesFontSize', cfg.FontSize);
    set(0, 'DefaultTextFontName', cfg.FontName);
    set(0, 'DefaultAxesFontName', cfg.FontName);
    set(gcf, 'color', [1 1 1]);

    colors = lines(length(cfg.categories));
    x = d.pairRep.all_time(d.pairRep.included_time);
    
    for i=1:length(cfg.categories)
        
        cat = cfg.categories{i};
        sem = std(d.meanAcc.(cat).pairRep, 0, 2)' ./ sqrt(cfg.n);
        mean_acc = mean(d.meanAcc.(cat).pairRep,2)';
        upper = mean_acc + sem;
        lower = mean_acc - sem;

        fill([x, fliplr(x)], ... 
            [upper, fliplr(lower)], colors(i, :), ...
            'EdgeColor', 'none', 'FaceAlpha', 0.3);
        h(i) = plot(x, mean_acc, 'Color', colors(i, :), 'LineWidth', 2);

    end

    % mark significant time points
    y_min = min([mean(d.meanAcc.kitchen.pairRep,2); mean(d.meanAcc.bathroom.pairRep,2)]) - 0.002; % etwas unter Minimalwert
    plot(x(p_val_fdr.kitchen), repmat(y_min,1,sum(p_val_fdr.kitchen)), 'b*', 'MarkerSize', 3);
    plot(x(p_val_fdr.bathroom), repmat(y_min-0.002,1,sum(p_val_fdr.bathroom)), 'r*', 'MarkerSize', 3);

    xlim([min(x)-min(x)-0.01, max(x)+0.01])
    xline(0, '--');
    yline(0.5, '--');

    % add labels and legend
    xlabel('Time (s)');
    ylabel('Mean decoding accuracy');
    legend(h, cfg.categories, 'Location', 'best');
    title(sprintf('Mean pairwise decoding accuracy\n with FDR-corrected significant time points'));

    % save figure
    if cfg.saving
        save_plot(fig, 'pairwise-decoding-with-sig-tp', cfg.figPath)
    end
end

end
