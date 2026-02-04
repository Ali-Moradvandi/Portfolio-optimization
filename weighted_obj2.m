function f = weighted_obj2(c, CAP, ROI, S, Phi, alpha, beta)
    % Normalize the weights
    % beta = 1 means linear capital expenditure scaling
    denom = sum(c.^beta .* CAP);            
    if denom == 0
        f = Inf;    % Avoid division by zero                      
        return;
    end
    % Compute weight vector
    w = (c.^beta .* CAP) / denom;
    % Expected return
    meanROI = ROI' * w;
    % Portfolio variance
    variance = w' * S * Phi * S * w;
    % Risk
    stdDev = sqrt(variance);
    % Weighted objective
    f = -alpha * meanROI + (1 - alpha) * stdDev;  
end