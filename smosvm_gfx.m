function [alpha, b] = smosvm_gfx(K, y, C, epsilon, plotfunc)
    if nargin < 5
        epsilon = 1e-3;
    end
    alpha = zeros(size(y));
    err = zeros(size(y));
    b = 0;
    smallval = 1e-8;
    conv = 0;
    plotfunc('done', {alpha, err});
    while ~conv
        conv = 1;
        partialConv = 0;
        fprintf('Going over unbound set...\n');
        while ~partialConv
            partialConv = 1;
            % Can't write it this way, because if empty it has size [0 1]
            % and Matlab executes the loop for one iteration.
            %for im=find(notatbound)
            notatbound = (alpha > smallval) & (alpha < C - smallval);
            nab = find(notatbound);
            nab = nab(randperm(length(nab)));
            %disp(size(findnab));
            if ~isempty(nab)
                for im=1:length(nab)
                    tryone(nab(im), 1);
                end
            end
        end
        fprintf('Partially converged!\n');
        all = randperm(length(alpha));
        for im=1:length(all)
            tryone(all(im), 0);
        end
    end
    
    checkkkt(alpha, y, K, C);
    
    function done = tryone(i, isnab)
        %fprintf('Checking conds for: %d\n', i);
        if ~kktok(i)
            %fprintf('First: %d\n', i);
            %if isnab
                %fprintf('1: from unbound set\n');
            %else
                %fprintf('1: from all\n');
            %end
            done = loop2(i);
            if done
                conv = 0;
                partialConv = 0;
                return;
            end
        end
    end
    
    function done = loop2(i)
        done = 0;
        %[~, allorder] = sort(err(i) - err, 'descend');
        %naborder = allorder(ismember(allorder, find(notatbound)));
        % Had a bug here - didn't have abs()
        %[~, nab2] = sort(double(err ~= 0) .* abs(err(i) - err), 'descend');
        [~, best] = max(double(err ~= 0) .* abs(err(i) - err));
        if err(best) == 0
            best = [];
        end
        %nab2 = nab2(1:length(find(err ~= 0)));
        %disp(nab2);
        nab2 = find(notatbound);
        ab2 = find(~notatbound);
        %nab2 = nab(randperm(length(nab2)));
        nab2 = nab2(randperm(length(nab2)));
        ab2 = ab2(randperm(length(ab2)));
        for i2=[best, nab2', ab2']
            %if ~isempty(best) && i2 == best
            %    fprintf('2: best\n');
            %elseif ismember(i2, nab2)
            %    fprintf('2: from unbound set\n');
            %else
            %    fprintf('2: from bound set\n');
            %end
            if i2 ~= i
                %fprintf('Trying %d\n', i2);
                done = do(i, i2);
                if done
                    break;
                end
            end
        end
    end

    function done = do(i, i2)
        plotfunc('do', [i, i2]);
        done = 0;
        if y(i) ~= y(i2)
            % satisfy "wlog" assumption
            if alpha(i) - alpha(i2) < 0
                temp = i;
                i = i2;
                i2 = temp;
            end
            lb = max(0, alpha(i2) - alpha(i));
            ub = min(C, C + alpha(i2) - alpha(i));
        else
            lb = max(0, alpha(i) + alpha(i2) - C);
            ub = min(C, alpha(i) + alpha(i2));
        end
        if ub == lb
            % Can't make progress, and can't calculate b.
            return;
        end
        %fprintf('Bounds: %g, %g\n', lb, ub);
        eta = 2*K(i,i2) - K(i,i) - K(i2,i2); % 2nd derivative
        err1 = geterr(i);
        err2 = geterr(i2);
        if eta < 0
            % function has a maximum
            newalpha = [0, alpha(i2) - y(i2) * (err1 - err2) / eta];
            newalpha(2) = min(ub, max(lb, newalpha(2)));            % clip
            %fprintf('New alpha is %g\n', newalpha(2));
        else
            % function is either increasing or decreasing (or flat)
            lbobj = obj(lb);
            ubobj = obj(ub);
            if lbobj == ubobj
                % pathological case - try another
                return;
            end
            if lbobj < ubobj
                newalpha = [0, lb];
            else
                newalpha = [0, ub];
            end
        end
        % make this all a bit more sensible ala pseudocode in SMO paper
        if newalpha(2) <= smallval
            newalpha(2) = 0;
        end
        if newalpha(2) >= C - smallval
            newalpha(2) = C;
        end
        if abs(newalpha(2) - alpha(i2)) < ...
                eps*(newalpha(2) + alpha(i2) + eps(newalpha(2) + alpha(i2)))
            % no progress
            return;
        end
        % This was where my mistake was! I didn't have alpha(i) here.
        newalpha(1) = alpha(i) + y(i) * y(i2) * (alpha(i2) - newalpha(2)); % find alpha1
        % find new b, update error cache, then save values
        btemp = zeros(2, 1);
        I=[i, i2];
        for ii=1:2
            btemp(ii) = err(I(ii)) + b;
            for jj=1:2
                btemp(ii) = btemp(ii) + ...
                    y(I(jj)) * (newalpha(jj) - alpha(I(jj))) * K(I(ii),I(jj));
            end
        end
        % always average them - it doesn't make any difference what type
        % they are.
        newb = mean(btemp);
        % update error cache
        for ii=1:find((alpha > smallval) & (alpha < C - smallval))
            err(ii) = err(ii) + b - newb;
            for jj=1:2
                err(ii) = err(ii) + y(I(jj)) * (newalpha(jj) - alpha(I(jj))) * K(ii,I(jj));
            end
        end
        err(i) = 0;
        err(i2) = 0;
        % save values
        alpha(i) = newalpha(1);
        alpha(i2) = newalpha(2);
        b = newb;
        done = 1;
        plotfunc('done', {alpha, err});
        %fprintf('Saved new values (%g, %g) - %g\n', newalpha(1), newalpha(2), newb);
    end

    function e = geterr(i)
        if alpha(i) > smallval && alpha(i) < C - smallval
            e = err(i);
        else
            e = f(i) - y(i);
        end
    end

    function ok = kktok(i)
        %fprintf('i has length %g\n', length(i));
        %fprintf('i: %g\n', i);
        %fprintf('alpha has length %g\n', length(alpha));
        %fprintf('alpha is %g\n', alpha(i));
        %if alpha(i) < -smallval || alpha(i) > C + smallval
        %    error('alpha %d is outside 0 - C', i);
        %end
        ok = 1;
        if alpha(i) > smallval && alpha(i) < C - smallval
            if abs(err(i)) > epsilon
                ok = 0;
            end
        elseif alpha(i) <= smallval
            if y(i) * f(i) < 1 - epsilon
                ok = 0;
            end
        else
            if y(i) * f(i) > 1 + epsilon
                ok = 0;
            end
        end
    end

    function val = f(i)
        val = (alpha .* y)' * K(:,i) - b;
    end
end