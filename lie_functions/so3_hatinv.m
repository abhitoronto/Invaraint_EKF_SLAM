function v = so3_hatinv(xi)
    
    v = [];
    if size(xi, 1)~= 3 || size(xi, 2) ~= 3
        disp("Bad entry: so3_hatinv: " + mat2str(xi));
        return
    end

    v = [xi(3,2);
         xi(1,3);
         xi(2,1)];
end