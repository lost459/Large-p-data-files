function [] = G1_General( D, p, lam, q, tau_end, n_stps)
%%% Calculate G2(tau) for a cavity of dimension D using iMPS methods. lam 
%%% is the parameter that detemines the form of the gain and loss
%%% operators. p is the parameter that determines the form of the ss
%%% density matrix. Note that this function requires use of the expv.m
%%% function.

    %%% Constructing the identity matrix on the D-dim and D^2-dim Hilbert
    %%% spaces
    I_D = sparse(1:D,1:D,ones(1,D),D,D);
    I_D2 = sparse(1:D^2,1:D^2,ones(1,D^2),D^2,D^2);
    
    %%% Constructing the identity matrix on the D-dim Hilbert space (this
    %%% is the left eigenmatrix of the T-matrix in the MPS language; i.e.,
    %%% a row vector, which is also the left eigenvector in the flattened
    %%% space notation):
    L1 = I_D;
    L1_flat = reshape(L1,[1,D^2]);

    %%% Create ss density matrix for cavity state. Note that in matlab
    %%% indexing starts with 1, while rho_n starts with n = 0.
    R1 = zeros(D);
    sum = 0;
    for ii = 1:D
        R1(ii,ii) = sin(pi*ii/(D+1))^p; % family of models defined in NP paper, but now to a generic power, p.
        sum = sum + R1(ii,ii);
    end
    rho_ss = R1/sum; % normalise
    rho_ss = sparse(rho_ss);

    %%% Constructing B0 matrix, which corresponds to the lain operator of
    %%% the laser. Note here that care must be taken with indexing.
    B0 = zeros(D);        
    for ii = 2:D
        B0(ii,ii-1) = sqrt(R1(ii,ii)/R1(ii-1,ii-1))^lam;   
    end
    B0 = sparse(B0);

    %%% Constructing B3 matrix, which corresponds to the loss operator of
    %%% the laser. Note here that care must be taken with indexing.
    B3 = zeros(D);        
    for ii = 2:D
        if lam~=0.0
            B3(ii-1,ii) = sqrt(R1(ii-1,ii-1)/R1(ii,ii))^(1-lam);
        else
            B3(ii-1,ii) = sqrt(R1(ii-1,ii-1)/R1(ii,ii))^(1+q/2);
        end
    end
    B3 = sparse(B3);
    
    %%% New ss density matrix for regular pump
    if q ~= 0
        rho_new = zeros(D);
        for ii = 1:D
            if ii == 1
                rho_new(ii,ii) = 1.0;
            end
            if ii == 2
                rho_new(ii,ii) = (B0(ii,ii-1)^2 - (q/2)*B0(ii,ii-1)^4)*rho_new(ii-1,ii-1)/B3(ii-1,ii)^2;
            end
            if ii == 3
                rho_new(ii,ii) = (B0(ii,ii-1)^2 + B3(ii-2,ii-1)^2 - (q/2)*B0(ii,ii-1)^4)*rho_new(ii-1,ii-1)/B3(ii-1,ii)^2 ...
                    + ((q/2)*(B0(ii-1,ii-2)^4 + B0(ii,ii-1)^2*B0(ii-1,ii-2)^2) - B0(ii-1,ii-2)^2)*rho_new(ii-2,ii-2)/B3(ii-1,ii)^2;
            end
            if ii > 3
                rho_new(ii,ii) = (B0(ii,ii-1)^2 + B3(ii-2,ii-1)^2 - (q/2)*B0(ii,ii-1)^4)*rho_new(ii-1,ii-1)/B3(ii-1,ii)^2 ...
                    + ((q/2)*(B0(ii-1,ii-2)^4 + B0(ii,ii-1)^2*B0(ii-1,ii-2)^2) - B0(ii-1,ii-2)^2)*rho_new(ii-2,ii-2)/B3(ii-1,ii)^2 ...
                    - (q/2)*(B0(ii-1,ii-2)^2*B0(ii-2,ii-3)^2)*rho_new(ii-3,ii-3)/B3(ii-1,ii)^2;
            end
        end
        % Get rid of negatives and normalise
        sm = 0;
        for ii = 1:D
            if rho_new(ii,ii) < 0
                rho_new(ii,ii) = 0;
            else
                sm = sm + rho_new(ii,ii);
            end
        end
        rho_ss = rho_new/sm; % normalise
        rho_ss = sparse(rho_ss);
    end
    % Require photon flux <Ldag L> = 1. Rescale B0 and B3 accordingly.
    R_flat = reshape(ctranspose(B3)*B3*rho_ss,[D^2,1]);
    N1 = L1_flat*R_flat;

    % MPS form of Liouvillain superoperator.
    Lg = ctranspose(B0)*B0;
    Ll = ctranspose(B3)*B3;
    LG = kron(conj(B0),B0) - 0.5*(kron(conj(Lg),I_D)+kron(I_D,Lg));
    LL = kron(conj(B3),B3) - 0.5*(kron(conj(Ll),I_D)+kron(I_D,Ll));
    % Normalise master equation
    L = (LG + (q/2)*LG^2 + LL) / N1;
    L = sparse(L)
    
    % MPS calculation of g1(tau).
    R_flat = kron(I_D,B3) * reshape(rho_ss,[D^2,1]) / sqrt(N1);
    
    % tol = 1e-10;
    tol = 1e-18
    SmallC = 1e-18;

    % whos L

    R_NEW = expv(0.0, L+SmallC*I_D2, R_flat, tol, D );
    G1 = L1_flat * kron(conj(B3),I_D) * R_NEW / sqrt(N1);

    M = [full(G1) D p lam q 0];
    dlmwrite('G1.txt', M, 'Delimiter', '\t', '-append', 'precision', 16);

    for ii = 1:n_stps

        ii/n_stps

        R_NEW = expv(tau_end/n_stps, L+SmallC*I_D2, R_NEW, tol, D);
        G1 = L1_flat * kron(conj(B3),I_D) * R_NEW / sqrt(N1);

        M = [full(G1) D p lam q ii*tau_end/n_stps];
        dlmwrite('G1.txt', M, 'Delimiter', '\t', '-append', 'precision', 16);

    end

end


