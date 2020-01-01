function f = se3_log(T)
    
    f = [];
    if size(T, 1)~= 4 || size(T, 2) ~= 4
        disp("Bad entry: se3_log: " + mat2str(T));
        return
    end
    
    phi = so3_log( T(1:3, 1:3) );
    rho = se3_jaco_inv(phi) * T(1:3, 4);
    
    f = [rho; ...
         phi ];
end