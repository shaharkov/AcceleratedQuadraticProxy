%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code implementing the paper "Accelerated Quadratic Proxy for Geometric Optimization", SIGGRAPH 2016.
% Disclaimer: The code is provided as-is for academic use only and without any guarantees. 
%             Please contact the author to report any bugs.
% Written by Shahar Kovalsky (http://www.wisdom.weizmann.ac.il/~shaharko/)
%            Meirav Galun (http://www.wisdom.weizmann.ac.il/~/meirav/)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef OptimSolverAcclQuadProx < OptimSolverIterative
    properties
        % internal variables
        x_prev;
        y;
        p;
        KKT;
        KKT_rhs;
        p_lambda;
        y_f;
        y_fgrad;
        t_start;
        t_init; % store intialization time
        f_count = 0; % store function evaluation count
        
        % solver parameters
        useAccelaration;
        useQuadProxy;
        useLineSearch;
        theta = [];
        
        % step size limiting
        useAccelarationStepSizeLimit = true;
        accelarationStepSizeLimitFactor = 0.5;
        accelarationStepSize;
        useLineSearchStepSizeLimit = true;
        lineSearchStepSizeLimitFactor = 0.5; %0.99;
        useLineSearchStepSizeMemory = true;
        lineSearchStepSizeMemoryFactor = 2^5.5;
        
        % line search parameters
        ls_alpha = 0.2;
        ls_beta = 0.5;
    end
    
    methods
        function obj = OptimSolverAcclQuadProx(tag, optimProblem, useAccelaration, useQuadProxy, useLineSearch)
            t_init_start = tic;
            % copy
            obj.useAccelaration = useAccelaration;
            obj.useQuadProxy = useQuadProxy;
            obj.useLineSearch = useLineSearch;
            % init solver
            obj.initSolver(tag, optimProblem);
            % precomputaions
            if obj.useQuadProxy
                KKT_mat = [obj.optimProblem.H obj.optimProblem.eq_lhs'; obj.optimProblem.eq_lhs sparse(obj.optimProblem.n_eq,obj.optimProblem.n_eq)];
            else
                KKT_mat = [speye(size(optimProblem.H)) obj.optimProblem.eq_lhs'; obj.optimProblem.eq_lhs sparse(obj.optimProblem.n_eq,obj.optimProblem.n_eq)];
            end
            obj.KKT = SparseLU(KKT_mat);
            % init internal variables
            obj.KKT_rhs = zeros(obj.optimProblem.n_vars+obj.optimProblem.n_eq,1);
            obj.x_prev = obj.x;
            obj.t = 1/obj.lineSearchStepSizeMemoryFactor;
            % store init time
            obj.t_init = toc(t_init_start);
        end
        
        function iterate(obj)
            % acceleration
            if obj.useAccelaration
                if obj.useAccelarationStepSizeLimit
                    obj.accelarationStepSize = min(obj.theta, obj.accelarationStepSizeLimitFactor*obj.optimProblem.getMaxStep(obj.x, (obj.x-obj.x_prev)));
                else 
                    obj.accelarationStepSize = obj.theta;
                end
                %obj.y = (1+obj.theta)*obj.x + obj.theta*obj.x_prev;
                obj.y = obj.x + obj.accelarationStepSize*(obj.x-obj.x_prev);
            else
                obj.y = obj.x;
            end
            
            % quadratic proxy minimization
            if obj.useLineSearch
                [obj.y_f, obj.y_fgrad] = obj.optimProblem.evaluateValueGrad(obj.y);
                obj.f_count = obj.f_count + 1;
            else
                obj.y_fgrad = obj.optimProblem.evaluateGrad(obj.y);
                obj.f_count = obj.f_count + 1;
            end
            obj.KKT_rhs(1:obj.optimProblem.n_vars) = -obj.y_fgrad;
            obj.p_lambda = obj.KKT.solve(obj.KKT_rhs);
            obj.p = obj.p_lambda(1:obj.optimProblem.n_vars);
            
            % initialize step size
            if obj.useLineSearchStepSizeMemory
                obj.t_start = min(obj.t * obj.lineSearchStepSizeMemoryFactor, 1);
            else
                obj.t_start = 1;
            end
            if obj.useLineSearchStepSizeLimit
                obj.t = min(obj.t_start, obj.lineSearchStepSizeLimitFactor*obj.optimProblem.getMaxStep(obj.y, obj.p));
            else
                obj.t = obj.t_start;
            end
            
            % line search
            if obj.useLineSearch
                [linesearch_cond_lhs, linesearch_cond_rhs] = computeLineSearchCond;
                while linesearch_cond_lhs>linesearch_cond_rhs
                    obj.t = obj.ls_beta*obj.t;
                    [linesearch_cond_lhs, linesearch_cond_rhs] = computeLineSearchCond;
                end
            end
            % fprintf('%e    %e    %e\n',obj.t_start,obj.lineSearchStepSizeLimitFactor*obj.optimProblem.getMaxStep(obj.y, obj.p),obj.t)
            
            % update
            obj.x_prev = obj.x;
            obj.x = obj.y + obj.t*obj.p;
            
            function [linesearch_cond_lhs, linesearch_cond_rhs] = computeLineSearchCond
                linesearch_cond_lhs = obj.optimProblem.evaluateValue(obj.y + obj.t*obj.p);
                obj.f_count = obj.f_count + 1;
                linesearch_cond_rhs = obj.y_f + obj.ls_alpha*obj.t*obj.y_fgrad'*obj.p;
            end
        end
        
        function setKappa(obj, kappa)
            obj.theta = (1-sqrt(1/kappa))/(1+sqrt(1/kappa));
        end
        
    end
    
end