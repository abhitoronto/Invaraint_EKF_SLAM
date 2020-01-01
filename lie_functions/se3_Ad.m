function [result] = se3_Ad(mat)
    
    result = [];
    if size(mat, 1) ~= 4 || size(mat, 2) ~= 4
        disp("Bad entry: se3_Ad: " + mat2str(mat));
        return
    end
    
    C = mat(1:3, 1:3);
    r = mat(1:3, 4);
    
    result = [ C         , so3_skew(r) * C;
               zeros(3,3), C               ];
    
end

