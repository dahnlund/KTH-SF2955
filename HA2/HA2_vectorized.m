clear variables; clc

%% Loading given data

Germany_inf = load("germany_infected.csv");Germany_rem= load("germany_removed.csv"); % German data

Iran_inf= load("iran_infected.csv");Iran_rem= load("iran_removed.csv"); % Iran data

Iran_Population= 84*10^6;
Germany_Population= 83*10^6;


Iran_Pop_vec= ones(length(Iran_inf),1)*Iran_Population;
Germany_Pop_vec= ones(length(Germany_inf),1)*Germany_Population; 

Iran_Sus= Iran_Pop_vec-Iran_inf - Iran_rem;
Germany_Sus= Germany_Pop_vec-Germany_inf-Germany_rem;

%% Implementation of Hybrid Sampler - Germany

% Parameters
phi = 0.995;
alfa = 2;
beta = 2;
sigma = 0.02;
M = 3;
a = 1;
b = 1;
P = Germany_Population;
p_si = @(lambdat, i) 1- exp(-lambdat.*i/P );

kappa = @(p_si, s) max(0,(1/phi -1)*s.*p_si);  % Max to avoid negative values

loglikl = @(kappa, x) sum(gammaln(x+kappa)-gammaln(x+1)-gammaln(kappa)+kappa*log(1-phi)+x*log(phi)); %Log-likelihood func without DeltaR since it is not used in the functions below;
% For the expression above the trick was used that log(factorial(x)) =
% log(gamma(x+1))

log_fullcond_t = @(kappa, x) loglikl(kappa, x); % The indicator condition is instead used in the loop

log_fullcond_lambda = @(lambdas,kappa, x) sum(alfa*log(beta)-gammaln(alfa)+(alfa-1)*log(lambdas)-beta*lambdas+loglikl(kappa, x));

acc_prob_l = @(lambda_s, lambda, kappa_cand, kappa, x) exp(log_fullcond_lambda(lambda_s, kappa_cand, x)-log_fullcond_lambda(lambda, kappa, x));
acc_prob_t = @(kappa_cand, kappa, x) exp(log_fullcond_t(kappa_cand, x)-log_fullcond_t(kappa, x));

N = 10000;
B = 2000; % Burn in
NB = N+B;
T = length(Germany_Sus);
num_t = 1;  % Number of breakpoints

% Define variables to save every iteration
t = zeros(NB,num_t);
lambdas = zeros(NB, num_t+1);
p_ir = zeros(NB,1);

% Initialization
saved_t(1,:) = sort(randsample(T,num_t));
lambdas(1,:) = 1;

delta_s = abs(diff(Germany_Sus));
delta_s= [0;delta_s];
delta_i = Germany_inf(2:end)-Germany_inf(1:end-1);
delta_i= [0;delta_i];

for m = 1:NB-1

    % Lambda update
    
    epslambda= randn(1,num_t+1);
    lambdas_cand = lambdas(m,:) + sigma *epslambda;
    if lambdas_cand < 0
        lambdas_cand = lambdas(m,:);  % If generating negative lambda, return to previous value
    end
    lambda_time_cand = lambda_t(lambdas_cand, saved_t(m,:), T);
    
    lambda_time = lambda_t(lambdas(m,:), saved_t(m,:),T);

    p_si_cand = p_si(lambda_time_cand, Germany_inf);
    p_si_1 = p_si(lambda_time, Germany_inf);

    kappa_cand = kappa(p_si_cand, Germany_Sus);
    kappa1 = kappa(p_si_1, Germany_Sus);
   

    prob_l = acc_prob_l(lambdas_cand, lambdas(m,:), kappa_cand, kappa1, delta_s);

    if rand < min(1,prob_l)
        lambdas(m+1,:) = lambdas_cand;
    else
        lambdas(m+1,:) = lambdas(m,:);
    end
    
    % Update t

    epst = randi([-M M], 1, num_t);
    t_cand = saved_t(m,:) + epst; % Create breakpoint candidate (i.e. t*)

    if t_condition(t_cand, T) == 0
        t_cand = saved_t(m,:);
    end

    lambda_time = lambda_t(lambdas(m+1,:), saved_t(m,:),T);
    lambda_time_cand = lambda_t(lambdas(m+1,:), t_cand, T);

    p_si_cand = p_si(lambda_time_cand, Germany_inf);
    p_si_1 = p_si(lambda_time, Germany_inf);

    kappa_cand = kappa(p_si_cand, Germany_Sus);
    kappa1 = kappa(p_si_1, Germany_Sus);
    
    prob_t = acc_prob_t(kappa_cand, kappa1, delta_s);
    if rand < min(1,prob_t)
        saved_t(m+1,:) = t_cand;
    else
        saved_t(m+1,:) = saved_t(m,:);
    end


    % Update p_ir
    
    p_ir(m) = betarnd(sum(delta_i+delta_s) - a,sum(Germany_inf-delta_s-delta_i)+b);
end

% Final estimation

Germany_est_p_ir = mean(p_ir(B:end));
Germany_est_lambdas = mean(lambdas(B:end,:), 1);
Germany_est_t = mean(saved_t(B:end,:), 1);
plot(lambdas)
title("Germany, Estimation of \lambda (one breakpoint)")
xlabel("Iterations")
legend(["\lambda_i before breakpoint" "\lambda_i after breakpoint"])
figure
plot(saved_t)
title("Germany, Estimation of breakpoint occuring")
xlabel("Iterations")
ylabel("Days")
legend("Breakpoint estimation")

% Estimation of lambda(t)
figure
est_lambda_t = lambda_t(Germany_est_lambdas, Germany_est_t, T);
plot(est_lambda_t(1:end-1), linewidth= 3)

title("Estimated \lambda as a function of time for Germany")
legend("\lambda(t)")
ylabel("\lambda")
xlabel("Days")

%% Implementation of Hybrid Sampler - Iran


% Parameters
phi = 0.995;
alfa = 2;
beta = 2;
sigma = 0.02;
M = 5;
a = 1;
b = 1;
P = Iran_Population;
p_si = @(lambdat, i) 1- exp(-lambdat.*i/P );

kappa = @(p_si, s) max(0,(1/phi -1)*s.*p_si);  % Max to avoid negative values

loglikl = @(kappa, x) sum(gammaln(x+kappa)-gammaln(x+1)-gammaln(kappa)+kappa*log(1-phi)+x*log(phi)); %Log-likelihood func;

log_fullcond_t = @(kappa, x) loglikl(kappa, x); % The indicator condition is instead used in the loop

log_fullcond_lambda = @(lambdas,kappa, x) sum(alfa*log(beta)-gammaln(alfa)+(alfa-1)*log(lambdas)-beta*lambdas+loglikl(kappa, x));

acc_prob_l = @(lambda_s, lambda, kappa_cand, kappa, x) exp(log_fullcond_lambda(lambda_s, kappa_cand, x)-log_fullcond_lambda(lambda, kappa, x));
acc_prob_t = @(kappa_cand, kappa, x) exp(log_fullcond_t(kappa_cand, x)-log_fullcond_t(kappa, x));

N = 10000;
B = 2000; % Burn in
NB = N+B;
T = length(Iran_Sus);
num_t = 1;

% Define variables to save every iteration
saved_t = zeros(NB,num_t);
lambdas = zeros(NB, num_t+1);  % Number of lambdas will be one more than the number of breakpoints
p_ir = zeros(NB,1);

% Initialization
saved_t(1,:) = sort(randsample(T,num_t));
lambdas(1,:) = 1;

delta_s = abs(diff(Iran_Sus));
delta_s= [0;delta_s];
delta_i = Iran_inf(2:end)-Iran_inf(1:end-1);
delta_i= [0;delta_i];

for m = 1:NB-1

    % Lambda update
    
    epslambda= randn(1,num_t+1);
    lambdas_cand = lambdas(m,:) + sigma *epslambda;
    if lambdas_cand < 0
        lambdas_cand = lambdas(m,:);  % If generating negative lambda, return to previous value
    end
    lambda_time_cand = lambda_t(lambdas_cand, saved_t(m,:), T);
    
    lambda_time = lambda_t(lambdas(m,:), saved_t(m,:),T);

    p_si_cand = p_si(lambda_time_cand, Iran_inf);
    p_si_1 = p_si(lambda_time, Iran_inf);

    kappa_cand = kappa(p_si_cand, Iran_Sus);
    kappa1 = kappa(p_si_1, Iran_Sus);
   

    prob_l = acc_prob_l(lambdas_cand, lambdas(m,:), kappa_cand, kappa1, delta_s);

    if rand < min(1,prob_l)
        lambdas(m+1,:) = lambdas_cand;
    else
        lambdas(m+1,:) = lambdas(m,:);
    end
    
    % Update t

    epst = randi([-M M], 1, num_t);
    t_cand = saved_t(m,:) + epst; % Create breakpoint candidate (i.e. t*)

    if t_condition(t_cand, T) == 0
        t_cand = saved_t(m,:);
    end

    lambda_time = lambda_t(lambdas(m+1,:), saved_t(m,:),T);
    lambda_time_cand = lambda_t(lambdas(m+1,:), t_cand, T);

    p_si_cand = p_si(lambda_time_cand, Iran_inf);
    p_si_1 = p_si(lambda_time, Iran_inf);

    kappa_cand = kappa(p_si_cand, Iran_Sus);
    kappa1 = kappa(p_si_1, Iran_Sus);
    
    prob_t = acc_prob_t(kappa_cand, kappa1, delta_s);
    if rand < min(1,prob_t)
        %disp("yes")
        saved_t(m+1,:) = t_cand;
    else
        %disp("no")
        saved_t(m+1,:) = saved_t(m,:);
    end


    % Update p_ir
    
    p_ir(m) = betarnd(sum(delta_i+delta_s) - a,sum(Iran_inf-delta_s-delta_i)+b);

end

% Final estimation

Iran_est_p_ir = mean(p_ir(B:end));
Iran_est_lambdas = mean(lambdas(B:end,:), 1);
Iran_est_t = mean(saved_t(B:end,:), 1);
figure
plot(lambdas)
title("Iran, Estimation of \lambda (one breakpoints)")
xlabel("Iterations")
legend(["\lambda_i before breakpoint" "\lambda_i after first breakpoint"])
figure
plot(saved_t)
title("Iran, Estimation of breakpoint occuring")
xlabel("Iterations")
ylabel("Days")
legend("Breakpoint 1 estimation")

% Estimation of lambda(t)
figure
est_lambda_t = lambda_t(Iran_est_lambdas, Iran_est_t, T);
plot(est_lambda_t(1:end-1), linewidth= 3)

title("Estimated \lambda as a function of time for Iran")
legend("\lambda(t)")
ylabel("\lambda")
xlabel("Days")