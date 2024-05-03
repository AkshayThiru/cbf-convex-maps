function [diff] = fast_setdiff(X, Y, idx_max)
    % Computes the set difference X \ Y, when X and Y contain positive
    % integers between 1 and idx_max. X and Y are column vectors.
    J_ = true(idx_max, 1);
    J_(Y) = false;
    diff = X(J_(X));
end
