function stimuli = compare_pairwiseDecodingRDM_with_stimuliRDM(cfg,d)
% evaluate input
if ~isfield(cfg, 'correlation_type'); cfg.correlation_type = 'pearson';end
if ~isfield(cfg, 'add_legend'); cfg.add_legend = true;end
if ~isfield(cfg, 'show_single_cate'); cfg.show_single_cate = false;end
if ~isfield(cfg, 'partial_cor'); cfg.partial_cor = false;end
if ~isfield(cfg, 'save_name'); cfg.save_name = 'compare_tp_RDMs_to_stimuli_RDMs';end
if ~isfield(cfg, 'xaxis_labels'); cfg.xaxis_labels = true;end
if ~isfield(cfg, 'dnns'); cfg.dnns = {cfg.dnn};end
if ~isfield(cfg, 'smoothing_window'); cfg.smoothing_window = 6;end
cfg.plot_rdm = false;
cfg.ISC_type = 'pairRep';

% time range
timepoints = d.pairRep.included_time;
stimuli = struct;
stimuli.all =struct;
for pr = cfg.predictor_RDMs
    pred = char(pr);
    stimuli.all.(pred) = struct;
    for c=cfg.categories
        category = char(c);
        stimuli.all.(pred).(category) = nan(cfg.n, numel(d.(cfg.ISC_type).included_time));
    end
end


% loop through subjects
for iSub = 1:length(cfg.subNums)
    % progress report
    disp(['subject ',  num2str(cfg.subNums(iSub))]);
    disp('')
    subID = ['sub', num2str(cfg.subNums(iSub))];
    stimuli.(subID)= struct;

    for c = cfg.categories
        category = char(c);

        stimuli.(subID).(category) = nan(numel(cfg.predictor_RDMs), numel(d.(cfg.ISC_type).included_time));

        % get canditate/predictor RDMs
        RDMs = struct;
        RDMs(1).name = ['sub', num2str(cfg.subNums(iSub))];
        RDMs(1).color = [0 0 0];
        RDMs(1).RDM = d.pairRep.(['sub', num2str(cfg.subNums(iSub))]).rdm; % (:, :, 1);
        labels = {RDMs.name};
        [RDMs, cfg.labels] = evaluate_predictor_RDMs(d, RDMs, labels, cfg, category);
        %stimuli_RDM = d.DNN.vgg16_imagenet.stimuli.(category).all_images.RDM;

        for iTp=1:length(timepoints)

            % get time point RDM
            if strcmp(category, 'bathroom')
                RDMs(1).RDM = d.pairRep.(['sub', num2str(cfg.subNums(iSub))]).rdm(1:50, 1:50, iTp);
            else
                RDMs(1).RDM = d.pairRep.(['sub', num2str(cfg.subNums(iSub))]).rdm(51:100, 51:100, iTp);
            end

            %             if cfg.partial_cor
            %                 [~, rMat, ~, cfg] = partial_cor_RDM(cfg, RDMs);
            %             else
            [~, rMat, ~] = cor_RDM(RDMs,cfg);
            %             end

            % save results
            stimuli.(subID).(category)(:, iTp) = rMat(2:end, 1);
        end % tp

        for j = 1:length(cfg.predictor_RDMs)
            pred = char(cfg.predictor_RDMs{j});
            stimuli.all.(pred).(category)(iSub,:) = stimuli.(subID).(category)(j,:);
        end

    end% category
end % sub

end