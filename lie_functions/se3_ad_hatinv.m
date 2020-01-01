function [result] = se3_ad_hatinv(mat)
    
    result = [];
    if size(mat, 1) ~= 4 || size(mat, 2) ~= 4
        disp("Bad entry: se3_hatinv: " + mat2str(mat));
        return
    end
    
    result = [ so3_hatinv(mat(1:3, 1:3));
               so3_hatinv(mat(1:3, 4:6)) ]; 
    
end

