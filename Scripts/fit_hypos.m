function [par, res] = fit_hypos(data,supp,M,model_type)
    %---------------------------
    % Model specific parameters.
    %---------------------------
    [ fitfn resfn degenfn psize numpar ] = getModelParam(model_type);

    %-----------------
    % Prepare storage.
    %-----------------
    n = size(data,2);
    par = zeros(numpar,M);
    res = zeros(n,M);

    for m=1:M
        % Pick a p-subset.
        psub = data(:,supp(:,m));
        
        if ~any(all(psub==0)) % controlla che non ci siano punti di supporto con le 3 coordinate = 0
            if ~contains(model_type,'fundamental')
                % Fit the model, and check for degeneracy.
                isdegen = feval(degenfn,psub);
            else
                % Check for degeneracy.
                [isdegen, F] = feval(degenfn,psub);
            end
    
            if (isdegen==1)
                error('Degenerate hypotesis!');
            end    
            
            if ~contains(model_type,'fundamental')
                % Fit the model on the p-subset.
                st = feval(fitfn,psub);
                % Compute residuals.
                ds = feval(resfn,st,data);
            else
                % Compute residuals.
                [ds, st] = feval(resfn,F,data);
            end
        else
            st = 0;
            ds = inf;
        end
        
        % Store.
        par(:,m) = st;
        res(:,m) = ds;
    end
end
