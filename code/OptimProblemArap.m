%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code implementing the paper "Accelerated Quadratic Proxy for Geometric Optimization", SIGGRAPH 2016.
% Disclaimer: The code is provided as-is for academic use only and without any guarantees. 
%             Please contact the author to report any bugs.
% Written by Shahar Kovalsky (http://www.wisdom.weizmann.ac.il/~shaharko/)
%            Meirav Galun (http://www.wisdom.weizmann.ac.il/~/meirav/)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef OptimProblemArap < OptimProblem
    properties
        % mesh (keep here?)
        V;
        F;
        dim;
        n_vert;
        n_tri;
        areas;
        w;
        
        % internal variables
        Tx;
        R;
        wTxMinusR;
    end
    
    methods
        function obj = OptimProblemArap(V, F, eq_lhs, eq_rhs, V0)
            % copy
            obj.V = V;
            obj.F = F;
            obj.eq_lhs = eq_lhs;
            obj.eq_rhs = eq_rhs;
            obj.dim = size(F,2)-1;
            obj.n_vert = size(V,1);
            obj.n_tri = size(F,1);
            % report
            obj.report(1,'Constructing %s (dim: %d   #vert: %d   #elem: %d)\n', class(obj), obj.dim, obj.n_vert, obj.n_tri);
            % compute transformations
            [obj.T,obj.areas] = computeMeshTranformationCoeffsMex(F, V);
            % set weights
            obj.w = kron(sqrt(obj.areas), ones(obj.dim^2,1));
            % set initial configuration
            obj.x0 = obj.initVertices(V0);
            % set quadratic proxy
            obj.H = obj.setQuadraticProxy();
            % finish construction
            obj.initProblem();
        end
        
        function x0 = initVertices(obj, V0)
            % feasible initialization -- provide input or a single local-global step
            if ~isempty(V0)
                x0 = V0;
            else 
                obj.Tx = obj.T*obj.V(:);
                obj.R = projectRotationMexFast(obj.Tx, obj.dim);
                x0 = solveConstrainedLS(obj.T, obj.R, obj.eq_lhs, obj.eq_rhs);
                x0 = reshape(x0, obj.n_vert, obj.dim);
            end
        end
        
        function [varargout] = evaluateFunctional(obj, x, doVal, doGrad, doHess)
            assert(~doHess, 'Hessian computation is not implementated for this functional');
            % evaluate
            obj.Tx = obj.T*x;
            obj.R = projectRotationMexFast(obj.Tx, obj.dim);
            obj.wTxMinusR = obj.w.*(obj.Tx-obj.R);
            n_arg = 0;
            if doVal
                n_arg = n_arg + 1;
                varargout{n_arg} = sum(obj.wTxMinusR.^2);
            end
            if doGrad
                n_arg = n_arg + 1;
                varargout{n_arg} = 2*((obj.w.*obj.wTxMinusR)'*obj.T)';
            end
        end
        
        function H = setQuadraticProxy(obj)
            wT = spdiag(obj.w)*obj.T;
            H = 2*(wT'*wT);
        end
        
        function t_max = getMaxStep(obj, x, p)
            t_max = inf; % do not enforce a t_max constraint for this functional
        end
        
        function e = evaluatePerElementEnergy(obj, x)
            % evaluate
            obj.Tx = obj.T*x;
            obj.R = projectRotationMexFast(obj.Tx, obj.dim);
            TxMinusR = reshape(obj.Tx-obj.R, obj.dim^2, []);
            e = sqrt(sum(TxMinusR.^2,1))';
        end
    end
    
end


