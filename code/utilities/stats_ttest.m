function is_significant = stats_ttest(cfg, d)

x = d.pairRep.all_time(d.pairRep.included_time)*1000;

if ~isfield(cfg, 'chance_level'); cfg.chance_level = 0.5; end
if ~isfield(cfg, 'plotting'); cfg.plotting = true; end
if ~isfield(cfg, 'xlim'); cfg.xlim = [-200, 500]; end
if ~isfield(cfg, 'alpha'); cfg.alpha = 0.05; end

% Initialize data storage
n_timepoints = length(d.pairRep.included_time);
p_val = struct();
t_stat = struct();
is_significant = struct();
p_val_fdr = struct();

for i = 1:length(cfg.categories)
    cat = cfg.categories{i};
    
    % Initialize data storage
    p_vals_cat = zeros(1, n_timepoints);
    t_vals_cat = zeros(1, n_timepoints);

    % compute t-test for each time point
    for t = 1:n_timepoints
        [~, p, ~, stats] = ttest(d.meanAcc.(cat)(t,:), cfg.chance_level, 'tail', 'right', 'alpha', cfg.alpha);
        p_vals_cat(t) = p;
        t_vals_cat(t) = stats.tstat;
    end

    % FDR
    [is_significant.(cat), ~, ~, adj_p] = fdr_bh(p_vals_cat, cfg.alpha);

    % save
    p_val.(cat) = p_vals_cat;
    t_stat.(cat) = t_vals_cat;
    p_val_fdr.(cat) = adj_p;

end

% stats.p_val = p_val;
% stats.t_stat = t_stat;
% stats.p_val_fdr = p_val_fdr;

if cfg.plotting

    fig = figure;
    set(gcf, 'Units','centimeters', 'Position',[2 2 40 25]);
    hold on;

    % graphic parameters
    set(0, 'DefaultTextFontSize', cfg.FontSize);
    set(0, 'DefaultAxesFontSize', cfg.FontSize);
    set(0, 'DefaultTextFontName', cfg.FontName);
    set(0, 'DefaultAxesFontName', cfg.FontName);
    set(gcf, 'color', [1 1 1]);

    % preallocate
    h = gobjects(1, length(cfg.categories));

    colors = lines(length(cfg.categories));
    pos = (max([mean(d.meanAcc.kitchen,2); mean(d.meanAcc.bathroom,2)]) + 0.016) * 100; % position to mark significant time points

    for i=1:length(cfg.categories)
        
        category = cfg.categories{i};
        sem = (std(d.meanAcc.(category), 0, 2)' ./ sqrt(cfg.n))* 100;
        mean_acc = mean(d.meanAcc.(category),2)' * 100;
        upper = mean_acc + sem;
        lower = mean_acc - sem;

        fill([x, fliplr(x)], ... 
            [upper, fliplr(lower)], colors(i, :), ...
            'EdgeColor', 'none', 'FaceAlpha', 0.3);
        h(i) = plot(x, mean_acc, 'Color', colors(i, :), 'LineWidth', 5, 'DisplayName', category);

        % mark significant time points
        plot(x(is_significant.(category)), repmat(pos-0.2*i,1,sum(is_significant.(category))),...
            'marker' ,'O', 'MarkerFaceColor', colors(i, :) , 'Color', colors(i, :), 'MarkerSize', 7);
    end

    xlim(cfg.xlim)
    ylim([49, pos+0.6])
    xline(0, '--', 'LineWidth', 2);
    yline(50, '--', 'LineWidth', 2);

    % add labels and legend
    xlabel('Time (ms)');
    ylabel('Decoding accuracy (%)');
    lgd = legend(h, 'Location', 'none');
    set(lgd, 'Position', [0.7 0.7 0.3 0.1]);
    legend('boxoff');
    %title(sprintf('Mean pairwise decoding accuracy'), "FontSize", cfg.FontSize+10);

    set(gca, 'FontWeight', 'bold', 'LineWidth', 2);
    set(gca, 'box', 'off');


    % save figure
    if cfg.saving
        save_plot(fig, 'pairwise-decoding-with-sig-tp', cfg.figPath)
    end
end

end
