function beh = beh_task(cfg)
%% define parameters
if ~isfield(cfg, 'chance_level'); cfg.chance_level = 0.5; end

%% array preparation
beh = struct;
beh.mean_acc.(cfg.categories{1}) = zeros(1, length(cfg.subNums));
beh.mean_acc.(cfg.categories{2}) = zeros(1, length(cfg.subNums));

%% compute accuracy
idx = 1;

for s=cfg.subNums
    % load data
    filepath=fullfile(cfg.sourcedataPath, ['sub-', num2str(s)], 'beh', ['sub-', num2str(s), '_task-main_events.tsv']);
    T = readtable(filepath, 'FileType', 'text', 'Delimiter', '\t');

    % Initialize data storage
    beh.(['sub', num2str(s)]).(cfg.categories{1}) = zeros(1,cfg.nBlocks/2);
    beh.(['sub', num2str(s)]).(cfg.categories{2}) = zeros(1,cfg.nBlocks/2);
    
    % Initialize indices and counter
    idx_k = 1;
    idx_b = 1;
    correct = 0;
    all_trials_with_question = 0;

    % compute block accuracy
    for i=1:20
        T_block = T(T.block == i, :);

        trials_with_question = height(T_block)-sum(isnan(T_block.accuracy));
        all_trials_with_question = all_trials_with_question + trials_with_question;

        if mod(i,2) == 0 % kitchen
            beh.(['sub', num2str(s)]).kitchen(idx_k) = sum(T_block.accuracy, "omitnan")/trials_with_question;
            idx_k = idx_k + 1;
        else % bathroom
            beh.(['sub', num2str(s)]).bathroom(idx_b) = sum(T_block.accuracy, "omitnan")/trials_with_question;
            idx_b = idx_b + 1;
        end

        correct = correct + sum(T_block.accuracy, "omitnan");

    end

    % save data
    beh.mean_acc.kitchen(idx) = mean(beh.(['sub', num2str(s)]).kitchen,2);
    beh.mean_acc.bathroom(idx) = mean(beh.(['sub', num2str(s)]).bathroom,2);
    beh.mean_acc.total(idx) = correct/all_trials_with_question;
    beh.(['sub', num2str(s)]).correct = correct;
    beh.(['sub', num2str(s)]).all_trials_with_question = all_trials_with_question;

    idx=idx+1;
end

%% t-test
% group category mean accuracy > 0.5
for i=1:length(cfg.categories)
    current_cat = cfg.categories{i};
    current_data = beh.mean_acc.(current_cat);
    [~, p, ~, ~] = ttest(current_data, cfg.chance_level, 'tail', 'right');
    p_val.(current_cat) = p;
end
p_vals = [p_val.bathroom, p_val.kitchen];
% fdr
is_sig = fdr_bh(p_vals, 0.05, 'pdep', 'yes'); % [is_sig, ~, ~, adj_p]

%% binom-test
% individual mean accuracy > 0.5
% P(X>=k) = 1 - P(X <= k-1) = 1 - binocdf(k-1, n, p)
p_vals_sub = zeros(1,cfg.n);
for idx = 1:cfg.n
    % data
    s = cfg.subNums(idx);
    current_data = beh.(['sub', num2str(s)]);
    p_vals_sub(idx) = 1 - binocdf(current_data.correct - 1, current_data.all_trials_with_question, 0.5); % 
end
% fdr
is_sig_sub = fdr_bh(p_vals_sub, 0.05, 'pdep', 'yes');

if sum(is_sig_sub)==cfg.n
    disp(['All ', cfg.n, ' subjects are significant better then chance level.']);
else
    idx_not_sig = find(is_sig_sub == 0);
    not_sig_subs = cfg.subNums(idx_not_sig);
    disp(['Only ', num2str(sum(is_sig_sub)), ' of ', num2str(cfg.n), ' subjects are significant better then chance level. No significant result/results: ', num2str(not_sig_subs), '.']);
end


%% plotting
if cfg.plotting

    % graphic parameters
    set(0, 'DefaultTextFontSize', cfg.FontSize);
    set(0, 'DefaultAxesFontSize', cfg.FontSize);
    set(0, 'DefaultTextFontName', cfg.FontName);
    set(0, 'DefaultAxesFontName', cfg.FontName);
    colors = lines(length(cfg.categories));

    % preallocate
    h = gobjects(1, length(cfg.categories));

    fig = figure;
    set(gcf, 'color', [1 1 1]); % white background
    hold on;
    for i=1:length(cfg.categories)
        % data
        current_cat = cfg.categories{i};
        current_data = beh.mean_acc.(current_cat);

        % plot group and individual mean-acc for given cat
        h(i) = bar(i, mean(current_data,2)', 'FaceColor', colors(i,:), 'EdgeColor', 'none');
        xi = i + 0.15 * sign(randn(1,length(cfg.subNums))) .* rand(1,length(cfg.subNums));
        scatter(xi, current_data, 40, 's', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'none', 'LineWidth', 1);

        % mark if significant
        if is_sig(i)
            plot(i, 1, 'k*', 'MarkerSize', 10);
        end

        % mark not significant subjects
        if sum(is_sig_sub)~=cfg.n
            scatter(xi(idx_not_sig), current_data(idx_not_sig), 40, 's', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'none', 'LineWidth', 1);
        end

    end

    xticks(1:2);
    xticklabels(cfg.categories);
    set(gca, 'TickDir', 'out');
    xlabel('Category');
    ylabel('Accuracy');
    ylim([0.5,1]);
    title('Accuracy behavioral task', 'Units', 'normalized', 'Position', [0.5, 1.05, 0]);
    legend(h, cfg.categories);

    % save plot
    if cfg.saving
        save_plot(fig, 'beh_acc', cfg.figPath);
    end

end
end