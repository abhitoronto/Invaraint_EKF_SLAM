function [result] = se3_cd(x)
    %%% This fucntion does the circle dot function
    result = [];
    if numel(x) ~= 4
        disp("Bad entry: se3_cd: " + mat2str(x));
        return
    end
    
    rho = x(1:3)';
    n = x(4);
    result = [ n*eye(3,3), -so3_skew(rho);
               zeros(1,6)];
end

