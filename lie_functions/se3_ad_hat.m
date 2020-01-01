function [result] = se3_ad_hat(X)

    result = [];
    if numel(X) ~= 6
        disp("Bad entry: ad_hat: " + mat2str(X))
        return
    end
    
    rho = X(1:3)';
    phi = X(4:6)';
    
    result = [ so3_hat(phi), so3_skew(rho);
               zeros(3,3)  ,  so3_skew(phi) ];
end
