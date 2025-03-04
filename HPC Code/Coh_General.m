function [ ] = Coh_General( D, p, lam, q )
%%% This function  calculates the coherence in the output field generated
%%% by a cavity with D levels, and Gain and Loss operators determined by
%%% the parameter lam (lam = 0.0 corresponds to a laser with constant
%%% gain, while lam = 1.0 corresponds to a laser with constant
%%% loss). SS density matrix is determined by p (heisenberg scaling is
%%% achieved for p >= 4, however prefactor 'a' in coh = a*mu^4 is maximised
%%% for p = 4).This function outputs the coherence ('coh'), which is 
%%% calculated using iMPS methods; B0 and B3, which are the MPS matricies 
%%% that correspond to the gain and loss operators, respectively; rho_ss, 
%%% is the density matrix of the cavity in the fock basis; L and L0 are 
%%% matricies that are relevant for calculating g1 and g2 using iMPS 
%%% methods.

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
    R1_flat = reshape(rho_ss,[D^2,1]); % flattened density matrix - this is the
    % right eigenmatrix of the T-matrix in MPS language.
    L1_flat0 = L1_flat;
    R1_flat0 = R1_flat;

    % Require photon flux <Ldag L> = 1. Rescale B0 and B3 accordingly.
    R_flat = reshape(ctranspose(B3)*B3*rho_ss,[D^2,1]);
    N1 = L1_flat*R_flat;
    
    % MPS form of Liouvillain superoperator.
    Lg = ctranspose(B0)*B0;
    Ll = ctranspose(B3)*B3;
    LG = kron(conj(B0),B0) - 0.5*(kron(conj(Lg),I_D)+kron(I_D,Lg));
    LL = kron(conj(B3),B3) - 0.5*(kron(conj(Ll),I_D)+kron(I_D,Ll));
    % Normalised master equation
    L = (LG + (q/2)*LG^2 + LL) / N1;
    Q_proj = I_D2 - R1_flat0*L1_flat0;
    L_proj = Q_proj*L*Q_proj;
    R_SigmaMinus = kron(I_D,B3) * R1_flat0  / sqrt(N1);      % this is an \rho_ss-operator-based vector.
    L_SigmaPos = L1_flat0 * kron(conj(B3),I_D) / sqrt(N1);  % this is an I-operator-based vector.

    SmallC = 1e-18; % add small number to L_proj so that the inverse exists
    X = -(L_proj+SmallC*I_D2)\I_D2;
    coh = 2*abs( L_SigmaPos * X * R_SigmaMinus );

    M = [full(coh) D p lam q];
    dlmwrite('coh.txt', M, 'Delimiter', '\t', '-append', 'precision', 16);

end

