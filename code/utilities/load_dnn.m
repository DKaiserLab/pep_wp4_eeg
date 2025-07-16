function [cfg, net] = load_dnn(cfg)
if contains(cfg.dnn,'vgg')

    if strcmp(cfg.dnn,'vgg16_imagenet')
        load(fullfile('..', 'vgg16_imagenet', 'vgg16.mat'), 'vgg16');
        net=vgg16;
        clear vgg16
    else
        error('network does not exist.');
    end

    % get layers
    if strcmp(cfg.layer_type,'early_mid_late')
        cfg.loi = [4,21,36];
    elseif strcmp(cfg.layer_type,'late')
        cfg.loi = 36;
    elseif strcmp(cfg.layer_type,'mid_late')
        cfg.loi = [21,36];
    elseif strcmp(cfg.layer_type,'all_output')
        %Get the Conv and FC Layers
        lx=0;
        for layer=1:length(net.Layers)
            if strfind(net.Layers(layer).Name,'conv')==1
                lx=lx+1;
                cfg.loi(lx)=layer;
            elseif strfind(net.Layers(layer).Name,'fc')==1
                lx=lx+1;
                cfg.loi(lx)=layer;
            end
        end
    elseif strcmp(cfg.layer_type,'all')
        cfg.loi = 1:41;
    end

else
    error('network does not exist.');
end
end