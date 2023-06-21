clc; clear all; close all; 
SETUP_PATH = 'setup.m';
run(SETUP_PATH)

 function [domeigval,domeigvec] = get_domeigcomp(J)
            [V,D] = eig(J);
            [sorted_eigval,idx] = sort(diag(D),'descend','ComparisonMethod','real');
            domeigval = sorted_eigval(1);
            domeigidx = idx(1);
            domeigvec = V(:,domeigidx);
        end

        function [Cd, Cr, Cl] = get_consumption_intensity(F, b0, n)
            Cd = F ./ repmat(b0.', n, 1);
            Cr = F ./ repmat(b0, 1, n);
            Cl = F ./ (b0 * b0.');
        end
        
        function J = compute_jacobian_triple(netInfo, b, b_balancing)
            if ~(netInfo.sd>=0 && netInfo.sd<=1 && netInfo.sr>=0 && netInfo.sr<=1 && netInfo.sl>=0 && netInfo.sl<=1)
                disp([netInfo.sd, netInfo.sr, netInfo.sl])
                error('Incorrect value of parameter sd, sr or, sl')
            end
            
            [Cd, Cr, Cl] = BaseClass.get_consumption_intensity(netInfo.F, b_balancing, netInfo.n);
            J_offdiag = netInfo.sd * Cd - netInfo.sr * Cr' ...
                + netInfo.sl * (Cl .* repmat(b, 1, netInfo.n)) ...
                - netInfo.sl * (Cl' .* repmat(b, 1, netInfo.n));
            J_diag = netInfo.sd * (diag(Cd) - diag(netInfo.sd * sum(Cd,1)))...
                + netInfo.sr * (diag(sum(Cr,2)) - diag(Cr))...
                + netInfo.sl * (diag(sum(repmat(b',netInfo.n,1) .* Cl, 2)) - netInfo.sl * diag(sum(repmat(b,1,netInfo.n) .* Cl, 1)))...
                - diag(netInfo.q./b_balancing + netInfo.r./b_balancing);
            J = J_offdiag - diag(J_offdiag) + J_diag;
        end
        
        function db = compute_db(netInfo, b)

            [Cd, Cr, Cl] = BaseClass.get_consumption_intensity(netInfo.F, netInfo.b0, netInfo.n);
            db = netInfo.sd * (sum(Cd .* repmat(b', netInfo.n, 1), 2) - sum(Cd .* repmat(b', netInfo.n, 1), 1)') ...
                + netInfo.sr * (sum(Cr .* repmat(b,1,netInfo.n),2) - sum(Cr .* repmat(b,1,netInfo.n),1)') ...
                + netInfo.sl * (sum(Cl .* repmat(b,1,netInfo.n) .* repmat(b',netInfo.n,1),2) - sum(Cl .* repmat(b',netInfo.n,1) .* repmat(b,1,netInfo.n),1)') ...
                - (netInfo.q./netInfo.b0 + netInfo.r./netInfo.b0).*b + netInfo.p;

            verifyVectorization = false;
            if verifyVectorization
                db2 = zeros(netInfo.n, 1);
                for i = 1:netInfo.n
                    donor = 0;
                    receipment = 0;
                    lv = 0;
                    rest = - (netInfo.q(i)/netInfo.b0(i) + netInfo.r(i)/netInfo.b0(i))*b(i) + netInfo.p(i);
                    for k = 1:netInfo.n
                        donor = donor + netInfo.sd * Cd(i,k) * b(k) - netInfo.sd * Cd(k,i) * b(i);
                        receipment = receipment + netInfo.sr * Cr(i,k) * b(i) - netInfo.sr * Cr(k,i) * b(k);
                        lv = lv + netInfo.sl * Cl(i,k) * b(i) *b(k) - netInfo.sl * Cl(k,i) * b(i) * b(k);
                    end
                    db2(i) = donor + receipment + lv + rest;
                end
                nd = norm(db-db2);
                if nd > 1e-3
                    warning(['Norm of difference ' num2str(nd)]);
                end
            end
            
            %[Mb, Ib] = max(b);
            %if Mb > (100 * max(netInfo.b0))
            %    disp(['Suspicious growth for ', netInfo.name, ' with ID: ', num2str(netInfo.netID), ', Comp: ', num2str(Ib), ' at sd = ', num2str(netInfo.sd), ' and sr = ', num2str(netInfo.sr)])
            %end
        end

  % Evaluate Jacobian at the empirically-observed
                    % biomass
   J = BaseClass.compute_jacobian_triple(netInfo, netInfo.b0, netInfo.b0);
   [domeigval,domeigvec] = BaseClass.get_domeigcomp(J);
