classdef SolverProjectorBD < SolverProjector

    properties
        K;
        lb;
        ub;
        dim;
        distortions;
        flips;
        minsv;
        maxsv;
    end
    
    methods
        function obj = SolverProjectorBD(F, V, eqLHS, eqRHS, K, lb, ub, x0, mode)
            t_start = tic;
            obj.eqLHS = eqLHS;
            obj.eqRHS = eqRHS;
            obj.x0 = x0;
            obj.mode = mode;
            obj.K = K;
            obj.lb = lb;
            obj.ub = ub;
            obj.dim = size(F,2)-1;
            [obj.T,areas] = computeMeshTranformationCoeffsMex(F, V);
            weights = kron(sqrt(areas),ones(obj.dim^2,1));
            obj.W = sparse(1:length(weights),1:length(weights),weights);
            %obj.W = 1;
            obj.initSolver();
            obj.report(1,'SolverProjectorBD is ready (%.3g secs)\n', toc(t_start));
        end
        
        function projectD_(obj)
            [obj.pTx, obj.distortions, obj.flips, obj.minsv, obj.maxsv] = projectBDMex(obj.Tx, obj.dim, obj.K, obj.lb, obj.ub);
        end
        
        function solve(obj, iter_max, tol_err, verbose)
            
            fprintf('-----------------------------------------------------------------\n');
            fprintf('BD PROJECTION (K=%g):\n', obj.K);
            fprintf('-----------------------------------------------------------------\n');
            fprintf('initial max dist %g,  flips %d,  infeasible %d\n', max(obj.distortions), nnz(obj.flips), nnz((obj.distortions>obj.K)|obj.flips));
            fprintf('(min sv %g,  max sv %d)\n', min(abs(obj.minsv)), max(obj.maxsv));
            fprintf('-----------------------------------------------------------------\n');
            for iter = 1:iter_max
                obj.iterate();
                err = mse(obj.tanNormal);
                if verbose
                fprintf('iter %d -- err: %g    dist: %g    flips: %g   time: %g sec (pD:%d%%, pL:%d%%)\n',iter,err,max(obj.distortions),nnz(obj.flips),obj.t_iter,round(100*obj.t_projectD/obj.t_iter), round(100*obj.t_projectLinear/obj.t_iter));
                end
                if (err<tol_err)&&(~nnz(obj.flips))
                    fprintf('err<tol_err --> stopping...\n');
                    break;
                end
            end
            fprintf('-----------------------------------------------------------------\n');
            fprintf('final max dist %g,  flips %d,  infeasible %d\n', max(obj.distortions), nnz(obj.flips), nnz((obj.distortions>obj.K)|obj.flips));
            fprintf('-----------------------------------------------------------------\n');
            
            
        end
    end
    
end

