function [result] = se3_jaco(phi)
    result = [];
    if numel(phi) ~= 3
        disp("Bad entry: se3_jaco: " + mat2str(phi));
        return
    end

    phi_norm = norm(phi);
    if(phi_norm < 0.00000001)
        result = eye(3,3);
    else
        result = sin(phi_norm)/phi_norm * eye(3,3) ...
                 - (1-cos(phi_norm))/(phi_norm^2) * so3_skew(phi) ...
                 + (phi_norm-sin(phi_norm))/(phi_norm^3) * (phi*phi');
    end

end

