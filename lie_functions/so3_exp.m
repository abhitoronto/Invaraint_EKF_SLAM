function [ result ] = so3_exp( x )
    % compute the exponential mapping of R^3 to SO(3)
	result = [];
    if numel(x) ~= 3
        disp("Bad entry: so3_exp: " + mat2str(x))
        return
    end
    
    theta=norm(x);
    if theta<0.0001
        result=eye(3);
    else
        omega = x/theta;
        result = cos(theta) * eye(3,3) ...
                 - sin(theta) * so3_skew(omega) ...
                 + (1 - cos(theta)) * omega * omega'; 
    end
end

