function y = so3_skew( x )

    y = [];
    if numel(x) ~= 3
        disp("Bad entry: so3_skew: " + mat2str(x));
        return;
    end

    y=[ 0    -x(3)  x(2);...
        x(3)  0    -x(1);...
       -x(2)  x(1)  0];
end

