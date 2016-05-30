%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code implementing the paper "Accelerated Quadratic Proxy for Geometric Optimization", SIGGRAPH 2016.
% Disclaimer: The code is provided as-is for academic use only and without any guarantees. 
%             Please contact the author to report any bugs.
% Written by Shahar Kovalsky (http://www.wisdom.weizmann.ac.il/~shaharko/)
%            Meirav Galun (http://www.wisdom.weizmann.ac.il/~/meirav/)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef (Abstract) OptimSolverIterative < OptimSolver
    properties
        % solver vars
        t = nan;
        
        %
        stopCntAccept = 5;
        tolXCnt = 0;
        tolFunCnt = 0;
        
    end
    
    properties (Dependent)
    end
    
    methods (Abstract)
        iterate(obj)
    end
    
    methods
        function obj = OptimSolverIterative
        end
        
        function log = solveTol(obj, TolX, TolFun, max_iter)
            obj.report(1,'Solving %s::%s  (TolX: %f   TolFun: %f   #max iter: %d)  ', obj.tag, class(obj), TolX, TolFun, max_iter);
            logInit;
            logState;
            % run num_iter iteration
            for iter = 1:max_iter
                t_iter_start = tic;
                obj.iterate();
                t_iter = toc(t_iter_start);
                % log
                logState;
                progBar;
                stop = stopCriteria;
                if stop
                    break;
                end
            end
            % truncate log
            log.iter = log.iter(1:iter+1);
            log.t_iter = log.t_iter(1:iter+1);
            log.f_count = log.f_count(1:iter+1);
            log.f = log.f(1:iter+1);
            log.t = log.t(1:iter+1);
            log.X = log.X(:,:,1:iter+1);

            function logInit
                iter = 0;
                t_iter = 0;
                log.tag = obj.tag;
                log.iter = nan(1, max_iter+1);
                log.t_iter = nan(1, max_iter+1);
                log.f_count = nan(1, max_iter+1);
                log.f = nan(1, max_iter+1);
                log.t = nan(1, max_iter+1);
                log.X = nan([size(obj.X), max_iter+1]);
                log.exitflag = 'MaxIter';
            end
            function logState
                log.iter(iter+1) = iter;
                log.t_iter(iter+1) = t_iter;
                log.f_count(iter+1) = obj.f_count;
                log.f(iter+1) = obj.optimProblem.evaluateValue(obj.x); % do not count, for logging only
                log.t(iter+1) = obj.t;
                log.X(:,:,iter+1) = obj.X;
            end
            function progBar
                if obj.verbose>=1
                    progressbar(iter, max_iter, 40);
                end
            end
            function stop = stopCriteria
                stop = false;
                if iter>1
                    % TolX
                    if norm(log.X(:,:,iter)-log.X(:,:,iter+1),'fro')<TolX*(1+norm(log.X(:,:,iter),'fro'))
                        obj.tolXCnt = obj.tolXCnt + 1;
                    else
                        obj.tolXCnt = 0;
                    end
                    if obj.tolXCnt>=obj.stopCntAccept
                        obj.report(2,' - stopped after %d iterations (TolX)\n', iter);
                        log.exitflag = 'TolX';
                        stop = true;
                    end
                    % TolFun
                    if abs(log.f(iter)-log.f(iter+1))<TolFun*(1+abs(log.f(iter)))
                        obj.tolFunCnt = obj.tolFunCnt + 1;
                    else
                        obj.tolFunCnt = 0;
                    end
                    if obj.tolFunCnt>=obj.stopCntAccept
                        obj.report(2,' - stopped after %d iterations (TolFun)\n', iter);
                        log.exitflag = 'TolFun';
                        stop = true;
                    end
                end
            end
        end
        
    end
    
end