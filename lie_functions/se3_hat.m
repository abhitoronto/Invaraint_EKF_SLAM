function [result] = se3_hat(X)

    result = [];
    if numel(X) ~= 6
        disp("Bad entry: se3_hat: " + mat2str(X))
        return
    end
    
    rho = reshape(X(1:3), [3,1]);
    phi = reshape(X(4:6), [3,1]);
    
    result = [ so3_skew(phi), rho;
               zeros(1,4)         ];
end

