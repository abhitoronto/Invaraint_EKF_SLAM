function f = se3_inv(T)
    
    f = [];
    if size(T, 1)~= 4 || size(T, 2) ~= 4
        disp("Bad entry: se3_log: " + mat2str(T));
        return
    end
    
    C =  T(1:3, 1:3);
    r = T(1:3, 4);
    
    f = [C'        , -C' * r;
         zeros(1,3),  1      ];
end
