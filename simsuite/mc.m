function [xEst, S, alpha, inno] = mc(runs, simLen, signal, model, ada, adaSwitch, x0, p0)
    xEst.val = zeros(runs*2, simLen);
    P.val = zeros(runs, simLen);
    
    S.S.val = zeros(runs, simLen);
    S.hapah.val = zeros(runs, simLen);
    S.hqh.val = zeros(runs, simLen);
    S.alphaR.val = zeros(runs, simLen);
    S.wincovar = zeros(runs, simLen);
    
    alpha.val = zeros(runs, simLen);
    
    inno.val = zeros(runs, simLen);
    
    for i = 1:runs
        %[i size(xEst.val( (2*i-1)+1 : 2*i ,:))]
        [z, ~] = genSig(signal, simLen, x0);
        [xEst.val( (2*i-1)+1 : (2*i)+1 ,:), pTemp, sTemp, alpha.val(i,:), S.wincovar(i,:), inno.val(i,:)] = kf(z, x0, p0, model, ada, adaSwitch);
        
        S.S.val(i,:) = sTemp(1,:);
        S.hapah.val(i,:) = sTemp(2,:);
        S.hqh.val(i,:) = sTemp(3,:);
        S.alphaR.val(i,:) = sTemp(4,:);
    end
    
    xEst.mean = mean(xEst.val(1:2:end,:));
    xEst.std = std(xEst.val(1:2:end,:));
    xEst.var = var(xEst.val(1:2:end,:));
    
    S.mean = mean(S.S.val);
    S.std = std(S.S.val);
    S.var = var(S.S.val);
    
    inno.mean = mean(inno.val);
    inno.std = std(inno.val);
    
    alpha.mean = mean(reshape((alpha.val - [ones(runs,simLen/2)*log2(signal.Rpre) ones(runs,simLen/2)*log2(signal.Rpost)]).^2, 1, []));
    alpha.std = std(reshape((alpha.val - [ones(runs,simLen/2)*log2(signal.Rpre) ones(runs,simLen/2)*log2(signal.Rpost)]).^2, 1, []));
    
    
end