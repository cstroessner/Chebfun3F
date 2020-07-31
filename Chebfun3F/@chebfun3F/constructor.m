function cf3F = constructor(cf3F, f, varargin)

% Parse the input
pref             = chebfunpref();
enableLog        = 0;
for k = 1:length(varargin)
    if strcmpi(varargin{k}, 'eps')
        pref.cheb3Prefs.chebfun3eps = varargin{k+1};
    elseif strcmpi(varargin{k}, 'enableLog')
        enableLog = varargin{k+1};
    end
end
if enableLog == 1
    fprintf('\nChebfun3F constructor: \n')
    evalsold = 0;
end

% Set Chebfun parameters:
tech             = pref.tech();
prefStruct       = pref.cheb3Prefs;
tpref            = tech.techPref;
grid             = tpref.minSamples;
maxSample        = tpref.maxLength;     % max polynomialDeg
maxSamplePhase1  = 363;                 % max coarseResolution (not implemented)
maxRank          = prefStruct.maxRank;  % max rank (not implemented)
pseudoLevel      = prefStruct.chebfun3eps;
passSampleTest   = prefStruct.sampleTest;
maxRestarts      = 10;

% Initialize
n                = [grid, grid, grid];  %coarseResolution
m                = n;                   %fineResolution
r                = [6,6,6];             %rank
tol              = pseudoLevel;
cf3F.numEvals    = 0;
cf3F.numRestarts = 0;
cheb             = @(i,n) -cos((i-1).*pi/(n-1));
reffun           = @(n) floor(sqrt(2)^(floor(2*log2(n)) + 1)) + 1;

%% Main Loop
happy = 0;
while ~happy
    
    %% Phase 1
    happyPhase1 = 0;
    while ~happyPhase1
        J = initializeIndexRandomly(r(2), n(2));
        K = initializeIndexRandomly(r(3), n(3));
        
        % Handle to evaluate tensor entries of T_c
        ff = @(i,j,k) f(cheb(i,n(1)),cheb(j,n(2)),cheb(k,n(3)));
        
        for iterations = 1:2
            happyPhase1 = 1;
            
            if enableLog == 1
                fprintf('ACA: tol = %.2d ', tol)
                fprintf('grid %2i %2i %2i ', n)
            end
            
            % ACA 1
            T1 = evalTensor(1:n(1),J,K,ff);
            cf3F.numEvals = cf3F.numEvals + numel(T1);
            T1 = matrizisation(T1,1);
            [~, tol] = getTol(T1, pseudoLevel, tol);
            [Uc, ~, ~, I,I2] = ACA(T1, tol, n(1));
            r(1) = size(I,2);
            JT1 = J;
            KT1 = K;
            
            % ACA 2
            T2 = evalTensor(I,1:n(2),K,ff);
            cf3F.numEvals = cf3F.numEvals + numel(T2);
            T2 = matrizisation(T2,2);
            [~, tol] = getTol(T2, pseudoLevel, tol);
            [Vc, ~, ~, J, J2] = ACA(T2, tol, n(2));
            r(2) = size(J,2);
            KT2 = K;
            IT2 = I;
            
            % ACA 3
            T3 = evalTensor(I,J,1:n(3), ff);
            cf3F.numEvals = cf3F.numEvals + numel(T3);
            T3 = matrizisation(T3,3);
            [reltol, tol] = getTol(T3, pseudoLevel, tol);
            [Wc, ~, ~, K, K2] = ACA(T3, tol, n(3));
            r(3) = size(K,2);
            IT3 = I;
            JT3 = J;
            
            if enableLog == 1
                fprintf('rank %2i %2i %2i evals %8i \n', r(1), r(2), r(3), cf3F.numEvals-evalsold)
                evalsold = cf3F.numEvals;
            end
            
            % Refine n
            breakFlag = 0;
            while r(1)*2*sqrt(2) > n(1)
                n(1) = reffun(n(1));
                breakFlag = 1;
            end
            while r(2)*2*sqrt(2) > n(2)
                n(2) = reffun(n(2));        for i = 1:3
            while r(i)*2*sqrt(2) > n(i)
                n(i) = reffun(n(i));
            end
        end
                breakFlag = 1;
            end
            while r(3)*2*sqrt(2) > n(3)
                n(3) = reffun(n(3));
                breakFlag = 1;
            end
            
            % Reinitialize after refinement
            if breakFlag == 1
                happyPhase1 = 0;
                if enableLog == 1
                    fprintf(' >> refine coarse grid \n')
                end
                break
            % Proceed to refinement if r < 2
            elseif min(r) < 2 
                break
            end
        end
    end
    
    % Catch the rank zero function:
    if size(I,1) == 0 || size(J,1) == 0 || size(K,1) == 0
        cf3F.U = chebfun(zeros([n(1),1]), pref);
        cf3F.V = chebfun(zeros([n(2),1]), pref);
        cf3F.W = chebfun(zeros([n(3),1]), pref);
        cf3F.C = 0;
    else
        %% Refine
        pref.chebfuneps = reltol;
        m = n;
        
        % Check if further refinement is necessary
        %U
        Uf = Uc;
        ct2 = createCT2(Uf);
        resolvedU = happinessCheck(ct2, [], ct2.coeffs, [], pref);
        if ~resolvedU
            m(1) = 2*m(1)-1;
        end
        
        %V
        Vf = Vc;
        ct2 = createCT2(Vf);
        resolvedV = happinessCheck(ct2, [], ct2.coeffs, [], pref);
        if ~resolvedV
            m(2) = 2*m(2)-1;
        end
        
        % W
        Wf = Wc;
        ct2 = createCT2(Wf);
        resolvedW = happinessCheck(ct2, [], ct2.coeffs, [], pref);
        if ~resolvedW
            m(3) = 2*m(3)-1;
        end
        
        % Add function evaluations and check again
        while ~resolvedU || ~resolvedV || ~resolvedW
            % function handle to evaluate T_f
            ff = @(i,j,k) f(cheb(i,m(1)),cheb(j,m(2)),cheb(k,m(3)));
            
            % Map the indices from T_c to T_f
            refFactor = [0, 0, 0];
            iter = 0;
            while min(refFactor) == 0
                iter = iter +1;
                if n(1)*iter-(iter-1) == m(1)
                    refFactor(1) = iter;
                end
                if n(2)*iter-(iter-1) == m(2)
                    refFactor(2) = iter;
                end
                if n(3)*iter-(iter-1) == m(3)
                    refFactor(3) = iter;
                end
            end
            ref = @(i, r) r*i-(r-1);
            
            % U
            Jr = ref(JT1, refFactor(2));
            Kr = ref(KT1, refFactor(3));
            if ~resolvedU
                Uold = Uf;
                [Ij,Ik] = ind2sub([size(Jr,2),size(Kr,2)],I2);
                Uf = zeros([m(1),size(I2,2)]);
                for i = 1:m(1)
                    for j = 1:size(I2,2)
                        if mod(i,2) == 0
                            Uf(i,j) = ff(i, Jr(Ij(j)), Kr(Ik(j)));
                            cf3F.numEvals = cf3F.numEvals+1;
                        else
                            Uf(i,j) = Uold((i+1)/2,j);
                        end
                    end
                end
                ct2 = createCT2(Uf);
                resolvedU = happinessCheck(ct2, [], ct2.coeffs, [], pref);
                if ~resolvedU
                    m(1) = 2*m(1)-1;
                end
            end
            
            % V
            Ir = ref(IT2, refFactor(1));
            Kr = ref(KT2, refFactor(3));
            if ~resolvedV
                Vold = Vf;
                [Ji,Jk] = ind2sub([size(Ir,2),size(Kr,2)],J2);
                Vf = zeros([m(2),size(J2,2)]);
                for i = 1:m(2)
                    for j = 1:size(J2,2)
                        if mod(i,2) == 0
                            Vf(i,j) = ff(Ir(Ji(j)), i, Kr(Jk(j)));
                            cf3F.numEvals = cf3F.numEvals+1;
                        else
                            Vf(i,j) = Vold((i+1)/2,j);
                        end
                    end
                end
                ct2 = createCT2(Vf);
                resolvedV = happinessCheck(ct2, [], ct2.coeffs, [], pref);
                if ~resolvedV
                    m(2) = 2*m(2)-1;
                end
            end
            
            % W
            Ir = ref(IT3, refFactor(1));
            Jr = ref(JT3, refFactor(2));
            if ~resolvedW
                Wold = Wf;
                [Ki,Kj] = ind2sub([size(Ir,2),size(Jr,2)],K2);
                Wf = zeros([m(3),size(K2,2)]);
                for i = 1:m(3)
                    for j = 1:size(K2,2)
                        if mod(i,2) == 0
                            Wf(i,j) = ff(Ir(Ki(j)), Jr(Kj(j)), i);
                            cf3F.numEvals = cf3F.numEvals+1;
                        else
                            Wf(i,j) = Wold((i+1)/2,j);
                        end
                    end
                end
                ct2 = createCT2(Wf);
                resolvedW = happinessCheck(ct2, [], ct2.coeffs, [], pref);
                if ~resolvedW
                    m(3) = 2*m(3)-1;
                end
            end
        end
        if enableLog == 1
            fprintf('Refinement: %3.i %3.i %3.i                         evals %8.i \n', m, cf3F.numEvals-evalsold)
            evalsold = cf3F.numEvals;
        end
        [~, tol] = getTol(Uf, pseudoLevel, tol);
        [~, tol] = getTol(Vf, pseudoLevel, tol);
        [~, tol] = getTol(Wf, pseudoLevel, tol);
        
        %% Phase 3
        
        % Compute factor matrices
        [Q1,~] = qr(Uf,0);
        I = DEIM(Q1);
        U = Q1/Q1(I,:);
        
        [Q2,~] = qr(Vf,0);
        J = DEIM(Q2);
        V = Q2/Q2(J,:);
        
        [Q3,~] = qr(Wf,0);
        K = DEIM(Q3);
        W = Q3/Q3(K,:);
        
        % Store U(x), V(y), W(z)
        cf3F.U = chebfun(U, pref);
        cf3F.V = chebfun(V, pref);
        cf3F.W = chebfun(W, pref);
        
        % Compute the core
        cf3F.C = evalTensor(I,J,K, ff);
        cf3F.numEvals = cf3F.numEvals + numel(cf3F.C);
        
        if enableLog == 1
            fprintf('Core & Assembling: max cond %.2d            evals %8.i \n', max([cond(Uf), cond(Vf), cond(Wf)]), cf3F.numEvals-evalsold)
        end
        evalsold = cf3F.numEvals;
    end
    
    % Sample Test
    if ( passSampleTest && cf3F.numRestarts < maxRestarts )
        cf3Fhandle = @(x,y,z) cf3F.feval(x,y,z);
        happy = sampleTest(f, cf3Fhandle, tol, enableLog);
        cf3F.numEvals = cf3F.numEvals + 30;
        if enableLog == 1
            fprintf('       evals %8.i\n', cf3F.numEvals - evalsold)
            evalsold = cf3F.numEvals;
            if happy == 1
                fprintf('                                                total %8.i\n', cf3F.numEvals)
            end
        end
    else
        happy = 1;
        if enableLog == 1
            fprintf('                                                total %8.i\n', cf3F.numEvals)
        end
    end
    
    % Restart
    if ~happy
        % Increase n
        n(1) = floor(sqrt(2)^(floor(2*log2(n(1))) + 1)) + 1;
        n(2) = floor(sqrt(2)^(floor(2*log2(n(2))) + 1)) + 1;
        n(3) = floor(sqrt(2)^(floor(2*log2(n(3))) + 1)) + 1;
        cf3F.numRestarts = cf3F.numRestarts + 1;
        
        if enableLog == 1
            fprintf('>> restart \n')
        end
        
        % Ensure r is large enougth for
        % (1,r,r) functions
        if r(1) > 1 || r(2) > 1 || r(2) > 1
            if r(1) < 3
                r(3) = max(6,2*r(3));
                r(2) = max(6,2*r(2));
            elseif r(2) < 2
                r(1) = max(6,2*r(1));
                r(3) = max(6,2*r(3));
            elseif r(3) < 2
                r(1) = max(6,2*r(1));
                r(2) = max(6,2*r(2));
            end
        end
        
        % very low-rank functions
        r(1) = max(r(1),3);
        r(2) = max(r(2),3);
        r(3) = max(r(3),3);
        
    end
end