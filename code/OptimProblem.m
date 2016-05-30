%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code implementing the paper "Accelerated Quadratic Proxy for Geometric Optimization", SIGGRAPH 2016.
% Disclaimer: The code is provided as-is for academic use only and without any guarantees. 
%             Please contact the author to report any bugs.
% Written by Shahar Kovalsky (http://www.wisdom.weizmann.ac.il/~shaharko/)
%            Meirav Galun (http://www.wisdom.weizmann.ac.il/~/meirav/)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef (Abstract) OptimProblem < handle
    properties
        % problem
        T;
        eq_lhs;
        eq_rhs;
        x0;
        n_vars;
        n_eq;
        H;
        
        % parameters
        verbose = 2;
    end
    
    methods (Abstract)
        evaluateFunctional(obj, x, doVal, doGrad, doHess)
        setQuadraticProxy(obj)
        getMaxStep(obj, x, p)
    end
    
    methods
        function obj = OptimProblem
        end
        
        function initProblem(obj)
            obj.n_vars = size(obj.T,2);
            obj.n_eq = size(obj.eq_lhs,1);
        end
        
        function f = evaluateValue(obj, x)
            f = evaluateFunctional(obj, x, true, false, false);
        end
        function f_grad = evaluateGrad(obj, x)
            f_grad = evaluateFunctional(obj, x, false, true, false);
        end
        function [f, f_grad] = evaluateValueGrad(obj, x)
            [f, f_grad] = evaluateFunctional(obj, x, true, true, false);
        end
        
        function [estH,err] = estimateHessian(obj, x)
            [estH,err] = hessian(@(x) obj.evaluateFunctional(x, true, false, false) ,x);
        end
        
        function report(obj,verbosity,varargin)
            if verbosity<=obj.verbose
                fprintf(varargin{:});
            end
        end
    end
    
end