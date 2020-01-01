function [result] = se3_ad_exp(X)

    result = [];
    if numel(X) ~= 6
        disp("Bad entry: ad_exp: " + mat2str(X))
        return
    end

    result = se3_Ad( se3_exp(X) );
end
