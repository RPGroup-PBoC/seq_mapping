clear all
close all

%Do a single energy matrix first

% Create a random emat
L = 20;
e_bj = rand(4,L);

% Plot betas vs es
betas = linspace(-20, 20, 1000);
log_rhos = [];
for i=1:numel(betas)
    beta = betas(i);
    Evals(i,1) = sum(sum(e_bj.*exp(-beta.*e_bj))./sum(exp(-beta.*e_bj)));
    logZ = sum(log(sum(exp(-beta*e_bj))));
    expect_e = sum(sum(e_bj.*exp(-beta.*e_bj))./sum(exp(-beta.*e_bj)));
    expect_esq = sum(sum(e_bj.^2.*exp(-beta.*e_bj))./sum(exp(-beta.*e_bj)));
    expect_desq = expect_esq - sum(sum(e_bj.*exp(-beta.*e_bj)).^2./sum(exp(-beta.*e_bj)).^2);
    log_rhos(i,1) = beta*expect_e + logZ - .5*log(2*pi*expect_desq);
end
minE = min(Evals);
maxE = max(Evals);

figure('windowstyle', 'docked')
plot(Evals, log_rhos);

% Initialize sequence at minimum energy
s0 = zeros(4,L);
s0 = e_bj == repmat(min(e_bj), 4, 1);
%for i=1:L
%    s0(randsample(4,1),i) = 1;
%end
%s0 = (1/4)*ones(4,L);

% Montecarlo in space of sequences
N = 50000;
s_old = s0;
logp_old = -Inf;
beta_old = 0;
es = zeros(N,1);
E_old = 0;
figure('windowstyle', 'docked')
for n=1:N
    
    % Perturb sequence
    s = s_old;
    for k=1:1
        i = randsample(L,1);
        b = randsample(4,1);
        s(:,i) = 0;
        s(b,i) = 1;
    end
    
    % Compute sequence energy
    E = sum(s(:).*e_bj(:));
    
    % Determine beta
    %fun = @(beta) beta_conditions(beta, E, e_bj);
    %beta = fzero(fun, beta_old);
    k = dsearchn(Evals,E);
    beta = betas(k);
    
    % Compute density of states
    Z = prod(sum(exp(-beta*e_bj)));
    expect_e = sum(sum(e_bj.*exp(-beta.*e_bj))./sum(exp(-beta.*e_bj)));
    expect_esq = sum(sum(e_bj.^2.*exp(-beta.*e_bj))./sum(exp(-beta.*e_bj)));
    expect_desq = expect_esq - sum(sum(e_bj.*exp(-beta.*e_bj)).^2./sum(exp(-beta.*e_bj)).^2);
    log_rho = beta*E + log(Z) - .5*log(2*pi*expect_desq);
    
    logp_new = -log_rho;
    
    if rand() < exp(logp_new - logp_old) && E > minE && E < maxE;
        E_old = E;
        s_old = s;
        beta_old = beta;
        logp_old = logp_new;
    end 
    es(n,1) = E;
    
    if mod(n,1000)==0
        disp(['MCMC, step ' num2str(n)])
        hist(es(1:n),100)
        drawnow;
    end
end
    