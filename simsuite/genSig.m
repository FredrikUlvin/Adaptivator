
function [z, zReal] = genSig(model, steps, x0)
%function [z, zReal] = genSig(phi, H, x0, Qpre, Qpost, Rpre, Rpost, steps)
processNoise = [wgn(size(model.phi,1), steps/2, model.Qpre, 'linear') ...
                    wgn(size(model.phi,1), steps/2, model.Qpost, 'linear')];
                
    measurementNoise = [wgn(size(model.H,1), steps/2, model.Rpre, 'linear') ...
                    wgn(size(model.H,1), steps/2, model.Rpost, 'linear')];
    
    %processNoise = [random('norm', 0, sqrt(Qpre), [size(phi,1), steps]) ...
    %                random('norm', 0, sqrt(Qpost), [size(phi,1), steps])];
    %measurementNoise = [random('norm', 0, sqrt(Rpre), [size(H,1),steps]) ... 
    %                    random('norm', 0, sqrt(Rpost), [size(H,1),steps])];
    meas = [];
    measReal = [];
    x = x0;
    xt = x0;
    
    for i = 1:steps
        xt =  model.phi * xt;
        measReal = [measReal model.H * xt];
        
        x = model.phi * x + [0 0; 0 1]*processNoise(:,i);
        meas = [meas model.H * x + measurementNoise(:,i)];
    end
    
    z = meas;
    zReal = measReal;
end