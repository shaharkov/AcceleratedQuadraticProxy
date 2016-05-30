classdef (Abstract) SolverProjector < handle
    
    properties
        % problem
        T; % lift operator
        W; % norm weights
        eqLHS;
        eqRHS;
        x0;
        % solver
        mode;
        usePreFactorization = true;
        nVars;
        nEq;
        x;
        Tx;
        pTx;
        tanNormal;
        tanLHS;
        tanRHS;
        preFactorization;
        TWW; % T'*W'*W
        TWWT; % T'*W'*W*T
        % temporaries
        % ?
        % log
        % ?
        % display / log
        verbose = 4;
        t_iter;
        t_projectD;
        t_projectLinear;
        t_factorization;
        % aux
        lambdaMultiTangent = 10;
    end
    
    properties (Dependent)
        y;
    end
    
    methods (Abstract)
        projectD_(obj)
    end
    
    methods
        function obj = SolverProjector
            % empty constructor
        end
        
        function value = get.y(obj)
            value = reshape(obj.x,size(obj.x0));
        end
        
        function initSolver(obj)
            obj.nVars = numel(obj.x0);
            obj.nEq = size(obj.eqLHS,1);
            obj.x = colStack(obj.x0);
            obj.Tx  = obj.T*obj.x;
            obj.projectD();
            obj.tanNormal = zeros(obj.nVars,1);
            obj.tanLHS = zeros(obj.nVars,1);
            obj.tanRHS = 0;
            obj.updateProblem();
        end
        
        function updateProblem(obj)
            t_start = tic;
            obj.TWW = obj.T'*obj.W'*obj.W;
            obj.TWWT = obj.TWW*obj.T;
            if (obj.usePreFactorization)
                obj.factorize();
            end
            obj.t_factorization = toc(t_start);
            obj.report(2,'Prectorization took (%.3g secs)\n', obj.t_factorization);
        end
        
        function factorize(obj)
            % construct KKT matrix
            LHS = [obj.TWWT, obj.eqLHS'; obj.eqLHS, sparse(obj.nEq,obj.nEq)];
            % factorize
            obj.preFactorization = SparseLU(LHS);
        end
        
        function projectD(obj)
            t_start = tic;
            obj.projectD_();
            obj.t_projectD = toc(t_start);
        end
        
        function projectLinear(obj)
            t_start = tic;
            switch obj.mode
                case SolverProjectorModeEnum.AltProj
                    temp_RHS = [obj.TWW*obj.pTx; obj.eqRHS];
                    if (obj.usePreFactorization)
                        temp_x_lambda = obj.preFactorization.solve(temp_RHS);
                    else
                        temp_LHS = [obj.TWWT, obj.eqLHS'; obj.eqLHS, sparse(obj.nEq,obj.nEq)];
                        temp_x_lambda = temp_LHS\temp_RHS;
                    end
                    obj.x = temp_x_lambda(1:obj.nVars);
                
                case SolverProjectorModeEnum.Tangent
                    if any(obj.tanNormal) % avoid solving linear system if projecting on D didn't do anything
                        obj.tanLHS = obj.tanNormal'*obj.T;
                        obj.tanRHS = obj.tanNormal'*obj.pTx;
                        if (obj.usePreFactorization)
                            temp_rhs = [obj.TWW*obj.pTx; obj.eqRHS];
                            temp_Au = [obj.tanLHS'; zeros(obj.nEq,1)];
                            %                             % compute lambda update
                            %                             temp_inv_LHS_RHS = obj.preFactorization.solve(temp_rhs);
                            %                             temp_inv_LHS_AuT = obj.preFactorization.solve(temp_Au);
                            %                             temp_lambda_u = (temp_Au'*temp_inv_LHS_RHS - obj.tanRHS) / (temp_Au'*temp_inv_LHS_AuT);
                            %                             % compute x_lambda
                            %                             temp_x_lambda_RHS = temp_rhs - temp_lambda_u*temp_Au;
                            %                             temp_x_lambda = obj.preFactorization.solve(temp_x_lambda_RHS);
                            Fm_c = obj.preFactorization.solve(temp_rhs);
                            Fm_n = obj.preFactorization.solve(temp_Au);
                            temp_x_lambda = Fm_c - (temp_Au'*Fm_c - obj.tanRHS)/(temp_Au'*Fm_n)*Fm_n;
                        else
                            error('Not implemented yet')
                        end
                    else % if tanNormal=0
                        temp_RHS = [obj.TWW*obj.pTx; obj.eqRHS];
                        if (obj.usePreFactorization)
                            temp_x_lambda = obj.preFactorization.solve(temp_RHS);
                        else
                            temp_LHS = [obj.TWWT, obj.eqLHS'; obj.eqLHS, sparse(obj.nEq,obj.nEq)];
                            temp_x_lambda = temp_LHS\temp_RHS;
                        end
                    end
                    obj.x = temp_x_lambda(1:obj.nVars);
                
                case SolverProjectorModeEnum.MultiTangent
                    if isa(obj,'SolverProjectorBD')
                        % split normals
                        manyTanNormals = sparse(1:length(obj.tanNormal), kron(1:(length(obj.tanNormal)/(obj.dim^2)),ones(1,obj.dim^2)), obj.tanNormal);
                        % remove inactive tangents 
                        manyTanNormals = manyTanNormals(:,any(manyTanNormals));
                        % normalize
                        norms = sqrt(sum(manyTanNormals.^2));
                        invNorms = sparse(1:length(norms),1:length(norms),1./norms);
                        manyTanNormals = manyTanNormals*invNorms;
                        % tangents
                        obj.tanLHS = manyTanNormals'*obj.T;
                        obj.tanRHS = manyTanNormals'*obj.pTx;
                        % solve
                        obj.x = solveConstrainedLS([obj.W*obj.T; sqrt(obj.lambdaMultiTangent)*obj.tanLHS],...
                            [obj.W*obj.pTx; sqrt(obj.lambdaMultiTangent)*obj.tanRHS], obj.eqLHS, obj.eqRHS);  
                        %fprintf('# active constraints for multitangents: %d\n',size(obj.tanLHS,1));
                    else
                        error('Problem type not supported');
                    end

                otherwise
                    error('invalid mode');
            end
            obj.t_projectLinear = toc(t_start);
        end
        
        function iterate(obj)
            t_start = tic;
            % lift
            obj.Tx = obj.T*obj.x;
            % project onto D
            obj.projectD();
            % compute normal (=error)
            obj.tanNormal = obj.Tx - obj.pTx;
            % project onto linear constraints
            obj.projectLinear();
            obj.t_iter = toc(t_start);
        end
        
        function iterateN(obj,n)
            for ii = 1:n
                obj.iterate();
            end
        end
        
        function solve(obj)
            error('Not implemented yet')
        end
        
        function report(obj,verbosity,varargin)
            if verbosity<=obj.verbose
                fprintf(varargin{:});
            end
        end
    end    
end

