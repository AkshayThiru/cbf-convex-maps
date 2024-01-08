classdef RigidEllipsoid < AbstractConvexSet
    % nz-dimensional ellipsoid with parameters: (position, orientation).
    % 
    % (p, R): (position, orientation).
    % Q: Positive-definite weight matrix.
    % C = {z: (z-p)^T R Q R^T (z-p) <= 1}.
    
    properties
        Q
    end
    
    methods
        function obj = RigidEllipsoid(nz, Q)
            obj = obj@AbstractConvexSet(nz + nz * nz, nz, 1);
            
            assert(isequal(size(Q), [nz, nz]));
            try chol(Q);
            catch
                disp('Q matrix is not symmetric positive definite');
            end
            obj.Q = Q;
        end
                
        function cons = A(obj, x, z)
            p = x(1:obj.nz);
            % R is vectorized column-first.
            R = reshape(x(obj.nz+1:end), obj.nz, obj.nz);
            
            cons = (z - p)' * R * obj.Q * R' * (z-p) - 1;
        end
        
        function [A, dAdx, dAdz, d2Adxz_y, d2Adzz_y] = derivatives(obj, x, z, y)
            p = x(1:obj.nz);
            R = reshape(x(obj.nz+1:end), obj.nz, obj.nz);
            
            % Note: A x = \sum_i A^i x_i = \sum_i (x_i I) A^i,
            %           = kron(x', I) vec(A).
            %       A'y = [y' A^1; y' A^2; ...] = kron(I, y') vec(A).
            A = (z - p)' * R * obj.Q * R' * (z-p) - 1;
            dAdz = 2 * (z - p)' * R * obj.Q * R';
            dAdx = [-dAdz, ...
                2 * (z - p)' * R * obj.Q * kron(eye(obj.nz), (z - p)')];
            d2Adzz_y = y(1) * 2 * R * obj.Q * R';
            d2Adxz_y = [-d2Adzz_y, ...
                y(1) * kron(2 * (z - p)' * R * obj.Q, eye(obj.nz)) + ...
                y(1) * 2 * R * obj.Q * kron(eye(obj.nz), (z - p)')];
        end
        
        function [] = plot_surf(obj, x)
            
        end
    end
end
