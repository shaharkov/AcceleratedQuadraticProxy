%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code implementing the paper "Accelerated Quadratic Proxy for Geometric Optimization", SIGGRAPH 2016.
% Disclaimer: The code is provided as-is for academic use only and without any guarantees. 
%             Please contact the author to report any bugs.
% Written by Shahar Kovalsky (http://www.wisdom.weizmann.ac.il/~shaharko/)
%            Meirav Galun (http://www.wisdom.weizmann.ac.il/~/meirav/)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef (Abstract) OptimSolver < handle
    properties
        % solver vars
        optimProblem;
        x;
        
        % solver parameters
        
        % log
        tag = '';
        verbose = 2;
        
        % default solver parameters
        default_max_iter = 500;
        default_TolX = 1e-10;
        default_TolFun = 1e-6;       
    end
    
    properties (Dependent)
        X;
    end
    
    methods (Abstract)
        solveTol(obj, TolX, TolFun, max_iter)
    end
    
    methods
        function obj = OptimSolver
            
        end
        
        function value = get.X(obj)
            value = reshape(obj.x,size(obj.optimProblem.x0));
        end
        
        function initSolver(obj, tag, optimProblem)
            % copy
            obj.tag = tag;
            obj.optimProblem = optimProblem;
            obj.x = obj.optimProblem.x0(:);
            % report
            obj.report(1,'Constructing %s::%s (#var: %d)\n', obj.tag, class(obj), obj.optimProblem.n_vars);
        end
        
        function log = solveN(obj,num_iter)
            log = obj.solveTol(0, 0, num_iter);
        end
        
        function log = solveDefault(obj)
            log = obj.solveTol(obj.default_TolX, obj.default_TolFun, obj.default_max_iter);
        end
        
        function report(obj,verbosity,varargin)
            if verbosity<=obj.verbose
                fprintf(varargin{:});
            end
        end
    end
end