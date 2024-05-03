function [union] = fast_union(X, Y, idx_max)
    % Computes the set union X \cup Y, when X and Y contain positive
    % integers between 1 and idx_max. X and Y are column vectors.
    J = (1: idx_max)';
    J_ = false(idx_max, 1);
    J_(X) = true;
    J_(Y) = true;
    union = J(J_);
end
