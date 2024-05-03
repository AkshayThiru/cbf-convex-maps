function [compl] = fast_2complement(X, Y, idx_max)
    % Computes the set complement (1:idx_max)' \ (X \cup Y), when X and Y
    % contain positive integers between 1 and idx_max. X and Y are column
    % vectors.
    J = (1: idx_max)';
    J_ = true(idx_max, 1);
    J_(X) = false;
    J_(Y) = false;
    compl = J(J_);
end
