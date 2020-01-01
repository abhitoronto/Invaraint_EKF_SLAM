function f = so3_log(R)
    
    f = [];
    if size(R, 1)~= 3 || size(R, 2) ~= 3
        disp("Bad entry: so3_log: " + mat2str(R));
        return
    end
    
    phi = real(acos(1/2*(trace(R)-1)));
%     [V, D] = eig(R);
%     indices = find( abs(real(diag(D)) - 1) < 0.01 );
%     if isempty(indices)
%         disp("oops");
%     end
%     a = real(V(:, indices(1) ));
%     f = phi*a;
%     
%     if so3_exp(f) ~= R
%         f = -phi*a;
%     end
    if phi>0.000001
        f = so3_hatinv(phi/(2*sin(phi))*(R-R'));
    else
        f= so3_hatinv(R - eye(3));
    end
    
end