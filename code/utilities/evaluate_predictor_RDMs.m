function [RDMs,labels] = evaluate_predictor_RDMs(d, RDMs, labels, cfg, category)

% get lenght of RDMs and lables 
n_rdms = length(RDMs);
labels_rdms = length(labels);
if n_rdms == labels_rdms
    disp(' Start evaluating predictor RDMs')
else
    warning('RDMs and labels have not the same length')
end

% get dnn models that are available in matlab
matlab_dnns = cfg.dnns(~ismember(cfg.dnns, {'clip', 'dino'}));

%% human rated simiarlairties
for rdm_num = 1:numel(cfg.predictor_RDMs)

    %% DNN similarities for control drawings based on generated images

    for dnn = matlab_dnns
        dnn = char(dnn);

        % early
        if strcmp('control_early', cfg.predictor_RDMs{rdm_num})
            rdm_index = 1;
            RDMs(end + 1) = d.DNN.(dnn).control.(category).subject_mean(rdm_index);
            labels{end + 1} = cfg.predictor_RDMs{rdm_num};

        % medium
        elseif strcmp('control_medium', cfg.predictor_RDMs{rdm_num})
            rdm_index = ceil(length(d.DNN.(dnn).control.(category).subject_mean)/2);
            RDMs(end + 1) = d.DNN.(dnn).control.(category).subject_mean(rdm_index);
            labels{end + 1} = cfg.predictor_RDMs{rdm_num};

        % late
        elseif strcmp('control_late', cfg.predictor_RDMs{rdm_num})
            rdm_index = length(d.DNN.(dnn).control.(category).subject_mean);
            RDMs(end + 1) = d.DNN.(dnn).control.(category).subject_mean(rdm_index);
            labels{end + 1} = cfg.predictor_RDMs{rdm_num};
        end
    end


    %% DNN similarities for own drawings based on generated images

    for dnn = matlab_dnns
        dnn = char(dnn);

        % early
        if strcmp('typical_early', cfg.predictor_RDMs{rdm_num})
            rdm_index = 1;
            RDMs(end + 1) = d.DNN.(dnn).typical.(category).subject_mean(rdm_index);
            labels{end + 1} = cfg.predictor_RDMs{rdm_num};

        % medium
        elseif strcmp('typical_medium', cfg.predictor_RDMs{rdm_num})
            rdm_index = ceil(length(d.DNN.(dnn).typical.(category).subject_mean)/2);
            RDMs(end + 1) = d.DNN.(dnn).typical.(category).subject_mean(rdm_index);
            labels{end + 1} = cfg.predictor_RDMs{rdm_num};

        % late
        elseif strcmp('typical_late', cfg.predictor_RDMs{rdm_num})
            rdm_index = length(d.DNN.(dnn).typical.(category).subject_mean);
            RDMs(end + 1) = d.DNN.(dnn).typical.(category).subject_mean(rdm_index);
            labels{end + 1} = cfg.predictor_RDMs{rdm_num};
        end
    end

    %% DNN similarities for Photos

    for dnn = matlab_dnns
        dnn = char(dnn);

        % early
        if strcmp('photos_early', cfg.predictor_RDMs{rdm_num})
            rdm_index = 1;
            RDMs(end + 1) = d.DNN.(dnn).photos.(category).subject_mean(rdm_index);
            labels{end + 1} = cfg.predictor_RDMs{rdm_num};

            % medium
        elseif strcmp('photos_medium', cfg.predictor_RDMs{rdm_num})
            rdm_index = ceil(length(d.DNN.(dnn).photos.(category).subject_mean)/2);
            RDMs(end + 1) = d.DNN.(dnn).photos.(category).subject_mean(rdm_index);
            labels{end + 1} = cfg.predictor_RDMs{rdm_num};

            % late
        elseif strcmp('photos_late', cfg.predictor_RDMs{rdm_num})
            rdm_index = length(d.DNN.(dnn).photos.(category).subject_mean);
            RDMs(end + 1) = d.DNN.(dnn).photos.(category).subject_mean(rdm_index);
            labels{end + 1} = cfg.predictor_RDMs{rdm_num};
        end
    end

    %% Other predictorsd

    for dnn = matlab_dnns
        dnn = char(dnn);

        % behavioral responses
        if strcmp('behResponses', cfg.predictor_RDMs{rdm_num})
            RDMs(end + 1) = d.behResponseRDM.(category).IS_RDM;
            labels{end + 1} = cfg.predictor_RDMs{rdm_num};
        elseif strcmp('betweenRDM', cfg.predictor_RDMs{rdm_num})
            RDMs(end + 1).RDM = d.betweenRDM.(cfg.currentVoi).RDM;
            RDMs(end).color = RDMs(end-1).color;
            RDMs(end).name = 'betweenRDM';
            labels{end + 1} = cfg.predictor_RDMs{rdm_num};
        elseif strcmp('withinRDM', cfg.predictor_RDMs{rdm_num})
            RDMs(end + 1).RDM = d.withinRDM.(cfg.currentVoi).RDM;
            RDMs(end).color = RDMs(end-1).color;
            RDMs(end).name = 'withinRDM';
            labels{end + 1} = cfg.predictor_RDMs{rdm_num};

        end
    end


%    
%     %% pyhton models
% 
%     % clip model
%     if strcmp('own_clip_model', cfg.predictor_RDMs{rdm_num})
%         RDMs(end + 1) = d.DNN.clip.Own_drawing_DNN_RDM.(category).subject_mean(end);
%         labels{end + 1} = 'Clip typical image';
%     elseif strcmp('control_clip_model', cfg.predictor_RDMs{rdm_num})
%         RDMs(end + 1) = d.DNN.clip.Control_drawing_DNN_RDM.(category).subject_mean(end);
%         labels{end + 1} = 'Clip control image';
%     elseif strcmp('original_own_clip_model', cfg.predictor_RDMs{rdm_num})
%         RDMs(end + 1) = d.DNN.clip.Original_own_drawing_DNN_RDM.(category).subject_mean(end);
%         labels{end + 1} = 'Original clip typical image';
%     elseif strcmp('original_control_clip_model', cfg.predictor_RDMs{rdm_num})
%         RDMs(end + 1) = d.DNN.clip.Original_control_drawing_DNN_RDM.(category).subject_mean(end);
%         labels{end + 1} = 'Original clip control image';
% 
%     % dino
%     elseif strcmp('own_dino_model', cfg.predictor_RDMs{rdm_num})
%         RDMs(end + 1) = d.DNN.dino.Own_drawing_DNN_RDM.(category).subject_mean(end);
%         labels{end + 1} = 'Dino typical image';
%     elseif strcmp('control_dino_model', cfg.predictor_RDMs{rdm_num})
%         RDMs(end + 1) = d.DNN.dino.Control_drawing_DNN_RDM.(category).subject_mean(end);
%         labels{end + 1} = 'Dino control image';
%     elseif strcmp('original_own_dino_model', cfg.predictor_RDMs{rdm_num})
%         RDMs(end + 1) = d.DNN.dino.Original_own_drawing_DNN_RDM.(category).subject_mean(end);
%         labels{end + 1} = 'Original dino typical image';
%     elseif strcmp('original_control_dino_model', cfg.predictor_RDMs{rdm_num})
%         RDMs(end + 1) = d.DNN.dino.Original_control_drawing_DNN_RDM.(category).subject_mean(end);
%         labels{end + 1} = 'Original dino control image';
%     end

end


% check is RDMs and Lebsl have same length
if length(RDMs) == length(labels)
    disp('Predictor RDMs evaluated')
else
    warning('RDMs and labels have not the same length')
end
end