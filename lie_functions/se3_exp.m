function [result] = se3_exp(X)
    
    result = [];
    if numel(X) ~= 6
        disp("Bad entry: se3_exp: " + mat2str(X))
        return
    end
    
    X = reshape(X, 6, 1);
    phi = norm(X(4:6));
    zeta_hat = se3_hat(X);
    result = eye(4) + zeta_hat + (1-cos(phi))/(phi) * (zeta_hat^2) ...
             + (phi - sin(phi))/(phi^3) * (zeta_hat^3) ;
    
end

