function [xest, p, S, alpha, wincovar, inno] = kf(z, x0, p0, model, ada, adaSwitch)
    %% Init
    Q = eye(size(model.phi, 1)) * model.Qpre;
    R = eye(size(model.H, 1)) * model.Rpre;
    H = model.H;
    phi = model.phi;
    gamma = [1 0; 0 1];
    x = zeros(size(phi,1), length(z));
    
    S = zeros(4, length(z));
    
    K = zeros(size(phi,1), length(z));
    inno = zeros(1, length(z));
    alpha = zeros(1, length(z));
    wincovar = zeros(1,length(z));
    
    x(:,1) = x0;
    P = p0;
    
    %% Kalman Filter
    for i = 1:length(z)
        % Update Step
        % Because of the initial estimated x and P, we begin with the
        
        % Update step
        inno(i) = (z(i)) - H*x(:,i);
        
        % Estimator enter from left
        if adaSwitch == 1
            S(1,i) = (H*P*H' + ada.alpha*R);
            ada.update(inno(i), S(1,i));
            alpha(i) = ada.alpha;
            wincovar(i) = ada.covarWindow;
            S(1,i) = (H*P*H' + ada.alpha*R);
            S(2, i) = H*phi*P*phi'*H';
            S(3, i) = H*Q*H';
            S(4, i) = ada.alpha*R;
            
        else
            S(1,i) = (H*P*H' + R);
        end
        
        K(:,i) = P * H' * inv(S(1,i));
        x(:,i) = x(:,i) + K(:,i)*inno(i);
        P = (eye(size(phi,1)) - K(:,i)*H)*P;

        % Prediction step
        if i < length(z)
            x(:,i+1) = phi * x(:,i);
            P = phi * P * phi' + gamma'*Q*gamma;
        end
    end
    
    %% Returning values
    xest = x;
    p = P;
    s = S;
    
end
