function PairwiseDecodingPlots(cfg, res)

if ~isfield(cfg, 'makeBetweenComparison'); cfg.makeBetweenComparison = false;end
if ~isfield(cfg, 'labels'); cfg.labels = cfg.subNums;end
if ~isfield(cfg, 'plotMatrix'); cfg.plotMatrix = false;end
if ~isfield(cfg, 'dissimilarity'); cfg.dissimilarity = true;end
if ~isfield(cfg, 'bar'); cfg.bar = false;end

% fileName = fullfile(cfg.outputPath, 'group_level', 'RDM',...
%     'results_RDM_of_pairwise_decoding_reref_pca.mat');
% load(fileName);

% Get subjects and timepoints
fns = fieldnames(res);
subjects = fns(3:end);
timepoints = res.included_time;
timeseries = res.all_time;
x  = timeseries(timepoints);

% Initialize data storage
num_timepoints = numel(timepoints);
all_data = nan(cfg.n, num_timepoints);
all_rdm_data = nan(cfg.nTrials, cfg.nTrials, cfg.n, num_timepoints);
bathroom_accuracy = nan(cfg.n, num_timepoints);
kitchen_accuracy = nan(cfg.n, num_timepoints);

% Collect data
for i_sub = 1:cfg.n
    subID = subjects{i_sub};
    for tp = 1:num_timepoints
        %         tp_idx = timepoints(tp);
        all_data(i_sub, tp) = res.(subID).mean_accuracy(tp);
        all_rdm_data(:, :, i_sub, tp) = res.(subID).rdm(:,:,tp);
        bathroom_accuracy(i_sub, tp)=mean(squareform(res.(subID).rdm(1:50,1:50,tp)));
        kitchen_accuracy(i_sub, tp)=mean(squareform(res.(subID).rdm(51:100,51:100,tp)));
    end
end

% graphic parameters
set(0, 'DefaultTextFontSize', cfg.FontSize);
set(0, 'DefaultAxesFontSize', cfg.FontSize);
set(0, 'DefaultTextFontName', cfg.FontName);
set(0, 'DefaultAxesFontName', cfg.FontName);

%% bar plot
if cfg.bar
    % Compute mean and standard derivation for each tp
    % set parameters
    figure;
    hold on;
    set(gcf, 'color', [1 1 1]); % white background

    % compute standard error of the mean
    sem = std(all_data, 0, 1) ./ sqrt(cfg.n);
    mean_acc = mean(all_data, 1, 'omitnan');

    % Bar plot with error bars
    bar(x,mean_acc, 'FaceColor', 'flat');
    h = errorbar(x, mean_acc, sem, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
    set(h, 'lineWidth', 2);
    set(h, 'linestyle', 'none');
    %set(h, 'barWidth', 0.2);

    % Add horizontal line at chance level
    yline(0.5, '--r', 'LineWidth', 2);
    text(num_timepoints + 0.8, 0.5, 'Chance Level', 'Color', 'r', ...
        'FontSize', 12, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

    % Mark stimulus onset
    xline(0, '--k');

    % Customize plot
    xlim([min(x)-min(x)-0.01, max(x)+0.01])
    ylim([0.45,0.6]);
    set(gca, 'TickDir', 'out');
    xlabel('time [s]');
    ylabel('accuracy');

    % set(gca, 'linewidth', 2);
    title({'Mean pairwise decoding results across timepoints', 'with standard error of the mean'});
    set(gca, 'box', 'off');

    hold off;
end
%% rdm
if cfg.plotMatrix
    % take mean across subjects
    mean_all_rdm_data = squeeze(mean(all_rdm_data, 3));

    fig = figure('position',[1,1,1000,600], 'unit','centimeters');;
    % tiledlayout(4, 4);
    title('Mean pairwise decoding accuracy per timepoint');

    for tp = 1:num_timepoints
        %     if ismember(tp-1, (16:16:num_timepoints))
        %         figure;
        %         tiledlayout(2, 2);
        %     else
        %         nexttile;
        %     end

        nexttile;
        set(gcf, 'color', [1 1 1]); % white background
        % get RDM for tp
        tpRDM = mean_all_rdm_data(:, :, tp);
        alltpRDMs(:, tp) = reshape(tpRDM, [], 1);

        % plot RDM
        imagesc(tpRDM, [0.50, 0.65])
        colorbar;
        title(timeseries(timepoints(tp)));
    end

    % get inter-tp correlation
    nexttile
    corrtps = corr(alltpRDMs, 'type', 'Spearman', 'rows', 'pairwise');

    imagesc(corrtps, [-0.7, 0.7])
    colorbar;
    title('inter-timepoint correlation');
    sgtitle('Mean decoding accuracy', 'FontSize', 18);

    if cfg.saving
        save_plot(fig, 'mean-pairwise-decoding-acc-bar', cfg.figPath);
    end
end

%% Create bar plot with comparison of within and between category
% % funktioniert noch nicht
% diff_per_sub = nan(num_timepoints, cfg.n);
% mean_diff = nan(1, num_timepoints);
% std_diff = mean_diff;
% for i_tp = 1:num_timepoints
%     for i_sub = 1:cfg.n
%         % get within and between category correlation
%         tpRDM = all_rdm_data(:, :, i_sub, i_tp);
%         withinCate = [squareform(tpRDM(1:cfg.nTrials/2, 1:cfg.nTrials/2)), ...
%             squareform(tpRDM(cfg.nTrials/2 + 1:end, cfg.nTrials/2 + 1:end))];
%         betweenCate = reshape(tpRDM(cfg.nTrials/2 + 1:end, 1:cfg.nTrials/2), 1, []);
%         diff_per_sub(i_tp, i_sub) = mean(withinCate) - mean(betweenCate);
%     end
%     % take the mean
%     mean_diff(i_tp) = mean(diff_per_sub(i_tp, :), 'omitnan');
%     std_diff(i_tp) = std(diff_per_sub(i_tp, :));
% end
% 
% if cfg.makeBetweenComparison
%     % plot within - between difference
%     figure;
%     hold on;
% 
%     % Bar plot with error bars
%     bar_handle = bar(mean_diff, 'FaceColor', 'flat');
%     errorbar(1:num_timepoints, mean_diff, std_diff, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
% 
%     % Add horizontal line at chance level
%     yline(0, '--r', 'No Category difference', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'right');
% 
%     % Add jittered individual points
%     rng(0); % For reproducible jitter
%     jitter_amount = 0.1; % Adjust jitter spread
%     for i_tp = 1:num_timepoints
%         x_jitter = i_tp + (rand(cfg.n, 1) - 0.5) * jitter_amount;
%         scatter(x_jitter, diff_per_sub(i_tp, :), 30, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'b');
%         for i_sub = 1:cfg.n
%             text(x_jitter(i_sub) + 0.04, diff_per_sub(i_tp, i_sub), subjects{i_sub}, 'FontSize', 8);
%         end
%     end
% 
%     % Customize plot
%     xticks(1:num_timepoints);
%     xticklabels(timeseries(timepoints));
%     xlabel('time');
%     ylabel('Pariwise decoding diff');
%     ylim([min(min(diff_per_sub))-0.005, max(max(diff_per_sub))+0.005])
%     title('Within - Between category pairwise correaltion');
% 
%     hold off;
% end

% %% is-rdm
% if cfg.plotMatrix
%     figure;
%     tiledlayout(3, 3);
%     sgtitle('IS-RDMs across the whole RDM');
% 
%     for i_tp = 1:num_timepoints
%         if ismember((i_tp-1), (9:9:num_timepoints))
%             figure;
%             tiledlayout(3, 3);
%             sgtitle('IS-RDMs across the whole RDM');
%         end
%         nexttile;
%         set(gcf, 'color', [1 1 1]); % white background
% 
%         % make a matrix with vectorized RDMs
%         for i_sub = 1:cfg.n
%             RDMmat(:, i_sub) = squareform(all_rdm_data(:, :, i_sub, i_tp));
%         end
% 
%         % make and plot IS-RDM
%         cfg.correlation_type = 'spearman';
%         cfg.cell_label_style = 'coef';
%         cfg.new_figure = false;
%         %         if ~cfg.dissimilarity
%         %             cfg.MinColorValue = -0.2;
%         %             cfg.MaxColorValue = 0.2;
%         %         else
%         %             cfg.MinColorValue = 0.8;
%         %             cfg.MaxColorValue = 1.2;
%         %         end
%         
%         [~, mat_out, ~] = make_RDM(RDMmat, cfg);
%         %         mat_out(eye(size(mat_out)) == 1) = 0;
%         medianISC = median(squareform(mat_out), 'omitnan');
%         title(sprintf('%0.3fs (median: %.4f)', x(i_tp), medianISC));
%         %         title([timeseries(timepoints(i_tp)), ' (median: ', num2str(round(medianISC, 4)), ')']);
%     end
% end

%% Lineplot
% without error
fig=figure;
hold on;
set(gcf, 'color', [1 1 1]); % white background
accuracies = {bathroom_accuracy, kitchen_accuracy};
colors = lines(length(cfg.categories));

for i=1:length(cfg.categories)
    plot(x, mean(accuracies{i}, 1, 'omitnan'), 'color', colors(i, :), 'LineWidth', 2);
end

xlim([min(x)-min(x)-0.01, max(x)+0.01])
xlabel('time');
ylabel('mean-accuracy');

% Mark stimulus onset and chance
yline(0.5, '--');
xline(0, '--k');  

legend('bathroom', 'kitchen');
title('pairwise decoding mean accuracy');

if cfg.saving
    save_plot(fig, 'pairwise-decoding-mean-acc', cfg.figPath)
end

% with error
fig=figure; 
hold on;
set(gcf, 'color', [1 1 1]); % white background

for i=1:length(cfg.categories)

    sem = std(accuracies{i}, 0, 1) ./ sqrt(cfg.n);
    mean_acc = mean(accuracies{i}, 1, 'omitnan');
    upper = mean_acc + sem;
    lower = mean_acc - sem;

    fill([x, fliplr(x)], ... % timeseries(timepoints)
        [upper, fliplr(lower)], colors(i, :), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
    %hold on;
    h(i) = plot(x, mean_acc, 'Color', colors(i, :), 'LineWidth', 2);
end

xlim([min(x)-min(x)-0.01, max(x)+0.01])
xlabel('time');
ylabel('mean-accuracy');
% Mark stimulus onset and chance
yline(0.5, '--');
xline(0, '--k');  

legend(h, cfg.categories);
title('pairwise decoding mean accuracy');

if cfg.saving
    save_plot(fig, 'pairwise-decoding-mean-acc-sem', cfg.figPath)
end

end 