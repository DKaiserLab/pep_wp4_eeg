function stats = stats_ttest_stimuli(cfg, d)

if ~isfield(cfg, 'plotting'); cfg.plotting = true; end
if ~isfield(cfg, 'alpha'); cfg.alpha = 0.05; end
if ~isfield(cfg, 'smoothing_window'); cfg.smoothing_window = 6; end
x=d.pairRep.all_time(d.pairRep.included_time)*1000;
if ~isfield(cfg, 'xlim'); cfg.xlim = [-200, 500]; end
if ~isfield(cfg, 'ylim'); cfg.ylim = [-0.05, 0.1]; end


% Initialize data storage
n_timepoints = length(d.pairRep.included_time);
n_pos_timepoints = sum(d.pairRep.all_time(d.pairRep.included_time)>=0);
pos_tp_idx = find((d.pairRep.all_time(d.pairRep.included_time)>=0), 1, 'first');
p_val = struct();
%t_stat = struct();
is_significant = struct();
p_val_fdr = struct();

for l=1:length(cfg.RDM_to_partial_out)
    layer = cfg.RDM_to_partial_out{l};
    p_val.(layer) = struct;
    %t_stat.(layer) = struct;
    p_val_fdr.(layer) = struct;
    is_significant.(layer) = struct;
end

for l=1:length(cfg.RDM_to_partial_out)
    layer = cfg.RDM_to_partial_out{l};
    for icat = 1:length(cfg.categories)
        category = cfg.categories{icat};

        % Initialize data storage
        p_vals_cat = struct;
%         t_vals_cat = struct;
        p_vals_cat.(layer) = zeros(1, n_pos_timepoints);
%         t_vals_cat.(layer) = zeros(1, n_pos_timepoints);

        % compute t-test for each time point
        for iTp = 1:n_pos_timepoints
            idx = iTp + pos_tp_idx - 1;
            %p = signrank(d.stimuli.all.(layer).(category)(:,idx), 0, 'tail', 'right', 'alpha', cfg.alpha);
            [~, p, ~, ~] = ttest(d.stimuli.all.(layer).(category)(:,idx), 0, 'tail', 'right', 'alpha', cfg.alpha);
            p_vals_cat.(layer)(iTp) = p;
            %t_vals_cat.(layer)(iTp) = stats.tstat;
        end

        % FDR
        [is_significant.(layer).(category), ~, ~, adj_p] = fdr_bh(p_vals_cat.(layer), cfg.alpha);

        % save
        p_val.(layer).(category) = p_vals_cat.(layer);
        %t_stat.(layer).(category) = t_vals_cat.(layer);
        p_val_fdr.(layer).(category) = adj_p;

    end % cat
end % layer

% mean
for l=1:length(cfg.RDM_to_partial_out)
    layer = cfg.RDM_to_partial_out{l};
    y_all = (d.stimuli.all.(layer).bathroom + d.stimuli.all.(layer).kitchen)/2;
    p_mean = zeros(1, n_pos_timepoints);
%     t_mean = zeros(1, n_pos_timepoints);

    for iTp = 1:n_pos_timepoints
        idx = iTp + pos_tp_idx - 1; 
        [~, p_mean(iTp), ~, ~] = ttest(y_all(:, idx), 0, 'tail', 'right', 'alpha', cfg.alpha);
        %p_mean(iTp) = signrank(y_all(:, idx), 0, 'tail', 'right', 'alpha', cfg.alpha);
        %t_mean(iTp) = stats.tstat;
    end
    
    % fdr
    [is_significant.(layer).all,~,~, adj_p] = fdr_bh(p_mean, cfg.alpha);

    p_val.(layer).all = p_mean;
    p_val_fdr.(layer).all = adj_p;

end % RDM

stats = struct;
stats.pval = p_val;
stats.pval_FDR = p_val_fdr;


%% plotting
if cfg.plotting
    % evtl. noch sem hinzufügen

    % graphic parameters
    set(0, 'DefaultTextFontSize', cfg.FontSize);
    set(0, 'DefaultAxesFontSize', cfg.FontSize);
    set(0, 'DefaultTextFontName', cfg.FontName);
    set(0, 'DefaultAxesFontName', cfg.FontName);
    set(gcf, 'color', [1 1 1]);

    fig = figure;
    set(gcf, 'Units','centimeters', 'Position',[2 2 40 25]);
%     fig.Position = [100 100 500 500]; 
    tiledlayout(length(cfg.RDM_to_partial_out),1);

    x=d.pairRep.all_time(d.pairRep.included_time)*1000;
    colors = lines(length(cfg.categories));
    
    for i=1:length(cfg.RDM_to_partial_out)
        nexttile; hold on;
        
        %preallocate
        lower = nan(length(cfg.categories),n_timepoints);
        upper = nan(length(cfg.categories),n_timepoints);
        p = gobjects(1,length(cfg.categories));

        layer = cfg.RDM_to_partial_out{i};
        for cat=1:length(cfg.categories)
            category = cfg.categories{cat};

            % compute and plot mean
            mean_cat = mean(d.stimuli.all.(layer).(category),1);
            p(cat)=plot(x, mean_cat, 'LineWidth', 5, 'Color', colors(cat,:));

            % plot standard error of the mean
            sem = std(d.stimuli.all.(layer).(category),0, 1)/sqrt(cfg.n);
            upper(cat,:) = mean_cat + sem;
            lower(cat,:) = mean_cat - sem;
            fill([x, fliplr(x)], ...
                [upper(cat,:), fliplr(lower(cat,:))], colors(cat,:),...
                'EdgeColor', 'none', 'FaceAlpha', 0.3, 'HandleVisibility', 'off');
        end
       
        % mark significant time points
        for cat=1:length(cfg.categories)
            category = cfg.categories{cat};
            h_idx = logical([false(1, pos_tp_idx-1), is_significant.(layer).(category)]);
            plot(x(h_idx), repmat(max(max(upper))+0.005*cat, 1, sum(h_idx)), ...
                 'marker' ,'O', 'LineStyle', 'none', 'MarkerFaceColor', colors(cat,:), 'Color', colors(cat,:), 'MarkerSize', 7, 'HandleVisibility','off');
        end

        title(strrep(layer, '_', '-'));
        ylim(cfg.ylim);
        xlim(cfg.xlim);
        yline(0,'--', 'LineWidth',2);
        xline(0,'--', 'LineWidth',2);
        xlabel('Time (ms)');
        ylabel({'Partial correlation', '1 - Spearman''s R'});
        set(gca, 'box', 'off');

    end % RDM
    sgtitle('Time-resolved RSA');
    lg = legend(p, cfg.categories);
    legend('boxoff');
    lg.Layout.Tile = 'east';
    % save plot
    if cfg.saving
        fig_path = fullfile(cfg.figPath, ['exp_', num2str(cfg.exp_num)], 'compare_tp_RDMs_to_stimuli_RDMs');
        save_plot(fig, cfg.save_name, fig_path);
    end

    %%
    fig = figure;
    fig.Position = [100 100 1000 500]; 
    hold on;
    p = gobjects(1, length(cfg.RDM_to_partial_out));
    for i = 1:length(cfg.RDM_to_partial_out)

        layer = cfg.RDM_to_partial_out{i};
        % smooth data
        % plot mean
        y1 = mean(d.stimuli.all.(layer).bathroom,1);
        y2 = mean(d.stimuli.all.(layer).kitchen,1);
        y = mean([y1; y2]);
        y = smoothdata(y, 2, 'movmean', cfg.smoothing_window);

        % get color
        if endsWith(cfg.RDM_to_partial_out{i}, 'early') 
            clr = [0.6, 0.8, 0.2];
        elseif endsWith(cfg.RDM_to_partial_out{i}, 'intermediate') 
            clr = [0.0, 0.55, 0.27];
        elseif endsWith(cfg.RDM_to_partial_out{i}, 'late')
            clr = [0.0, 0.27, 0.13];
        end
        % plot mean
        p(i) = plot(x, y, 'color', clr, 'LineStyle', '-', 'LineWidth', 5, 'DisplayName', regexprep(cfg.RDM_to_partial_out{i}, '.*_', ''));

        % mark significant timepoints
        h_idx = logical([false(1, pos_tp_idx-1), is_significant.(layer).all]);
        plot(x(h_idx), repmat(0.06+0.01*i, 1, sum(h_idx)), ...
            'marker' ,'O', 'MarkerFaceColor', clr, 'Color', clr, 'MarkerSize', 7, 'LineStyle', 'none', 'HandleVisibility','off');
        
    end

    legend(p, 'Location','northeastoutside');
    legend('boxoff');
    ylim(cfg.ylim);
    xlim(cfg.xlim);
    yline(0,'--', 'LineWidth', 2, 'HandleVisibility', 'off');
    xline(0,'--', 'LineWidth', 2, 'HandleVisibility', 'off');
    xlabel('Time (ms)');
    ylabel('Spearman''s R');

%     title('Time-resolved RSA');
    set(gca, 'box', 'off');
    set(gca, 'FontWeight', 'bold', 'LineWidth', 2);

    % saving
    if cfg.saving
        fig_path = fullfile(cfg.figPath, ['exp_', num2str(cfg.exp_num)], 'compare_tp_RDMs_to_stimuli_RDMs');
        save_plot(fig, [cfg.save_name, '_mean'], fig_path)
    end
end

end
