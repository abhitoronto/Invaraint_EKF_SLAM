function [result] = se3_jaco_inv(phi)
    result = [];
    if numel(phi) ~= 3
        disp("Bad entry: se3_jaco: " + mat2str(phi));
        return
    end

    phi_norm = norm(phi);
    theta = phi_norm/2;
    if(phi_norm < 0.00000001)
        result = eye(3,3);
    else
        result = theta*cot(theta)*eye(3) ...
                 + (1 - theta*cot(theta))*(phi*phi')/(phi_norm^2)...
                 - theta*so3_skew(phi/phi_norm);
    end

end


