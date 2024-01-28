clear variables; clc

load("stations.mat"); load("RSSI-measurements.mat")

%% DEFINE PROBLEM

resampling = false;  % If true, perform resampling

%% Define parameters

dt = 0.5;
a = 0.6;
sigma = 0.5;
Z_space = [0 3.5 0 0 -3.5; 0 0 3.5 -3.5 0];
P = 1/20 * (ones(5,5) + 15*eye(5,5));

ksit_z = [dt^2/2; dt; 0];
ksit_w = [dt^2/2; dt; 1];

ksi_z = zeros(6,2);
ksi_z(1:3,1) = ksit_z;
ksi_z(4:6,2) = ksit_z;

ksi_w = zeros(6,2);
ksi_w(1:3,1) = ksit_w;
ksi_w(4:6,2) = ksit_w;

phit = [1 dt dt^2/2; 0 1 dt; 0 0 a];
phi = zeros(6,6);
phi(1:3,1:3) = phit;
phi(4:6,4:6) = phit;

%% TEST
m = 501;
X = zeros(6,m);

Z_indices = randsample(5,1,true);
Z = Z_space(:,Z_indices);
for n = 2:m
    X(:,n) = phi * X(:,n-1) + ksi_z * Z + ksi_w * randn(2,1)*sigma;
    Z_indices = P_sample(Z_indices);
    Z = Z_space(:, Z_indices);
end

N = 20000;
n = length(Y(1,:));

u = 90;
eta = 3;
zeta = 1.5;

Y_test = [u - 10*eta*log10(eucl([X(1,:); X(4,:)], pos_vec(:,1))); ...
    u - 10*eta*log10(eucl([X(1,:); X(4,:)], pos_vec(:,2)));...
    u - 10*eta*log10(eucl([X(1,:); X(4,:)], pos_vec(:,3)));...
    u - 10*eta*log10(eucl([X(1,:); X(4,:)], pos_vec(:,4)));...
    u - 10*eta*log10(eucl([X(1,:); X(4,:)], pos_vec(:,5)));...
    u - 10*eta*log10(eucl([X(1,:); X(4,:)], pos_vec(:,6)))] + zeta *randn(6,1);

part = diag([500 5 5 200 5 5]).^(1/2) * randn(6, N); % Initialization

tau = zeros(6,n);

p =  @(x,y) mvnpdf(y', [u - 10*eta*log10(eucl([x(1,:); x(4,:)], pos_vec(:,1))); ...
    u - 10*eta*log10(eucl([x(1,:); x(4,:)], pos_vec(:,2)));...
    u - 10*eta*log10(eucl([x(1,:); x(4,:)], pos_vec(:,3)));...
    u - 10*eta*log10(eucl([x(1,:); x(4,:)], pos_vec(:,4)));...
    u - 10*eta*log10(eucl([x(1,:); x(4,:)], pos_vec(:,5)));...
    u - 10*eta*log10(eucl([x(1,:); x(4,:)], pos_vec(:,6)))]', zeta^2*eye(6,6));

saved_w = zeros(N,n);
    
w = p(part, Y_test(:,1));
saved_w(:,1) = w;

tau(:,1) = [sum(part(1,:) .* w')./sum(w);...
            sum(part(2,:) .* w')./sum(w);...
            sum(part(3,:) .* w')./sum(w);...
            sum(part(4,:) .* w')./sum(w);...
            sum(part(5,:) .* w')./sum(w);...
            sum(part(6,:) .* w')./sum(w)];

for t = 2:n
    
    Z = Z_space(:,randsample(5, N,true));

    if resampling == true
        ind = randsample(N,N, true, w);
        part = phi * part(:,ind) + ksi_z * Z + ksi_w * randn(2,N)*sigma;
        l = log(p(part, Y_test(:,t)));
    else
        part = phi * part+ ksi_z * Z + ksi_w * randn(2,N)*sigma;
        l = log(p(part, Y_test(:,t))) + log(w);
    end

    w = exp(l - max(l));
    saved_w(:,t) = w;
    tau(:,t) = [sum(part(1,:) .* w')./sum(w);...
            sum(part(2,:) .* w')./sum(w);...
            sum(part(3,:) .* w')./sum(w);...
            sum(part(4,:) .* w')./sum(w);...
            sum(part(5,:) .* w')./sum(w);...
            sum(part(6,:) .* w')./sum(w)];
end

plot(pos_vec(1,:), pos_vec(2,:), 'o', LineWidth = 4)
hold on
plot(X(1,:), X(4,:), LineWidth = 2)
plot(tau(1,:), tau(4,:))
title("Prediction from generated X and corresponding Y")
legend(["Positions" "Generated motion" "Predicted motion"])

%% SIS with given Y data
N = 10000;
n = length(Y(1,:));

u = 90;
eta = 3;
zeta = 1.5;
part = diag([500 5 5 200 5 5]).^(1/2) * randn(6, N); % Initialization

tau = zeros(6,n);false
    
p =  @(x,y) mvnpdf(y', [u - 10*eta*log10(eucl([x(1,:); x(4,:)], pos_vec(:,1))); ...
    u - 10*eta*log10(eucl([x(1,:); x(4,:)], pos_vec(:,2)));...
    u - 10*eta*log10(eucl([x(1,:); x(4,:)], pos_vec(:,3)));...
    u - 10*eta*log10(eucl([x(1,:); x(4,:)], pos_vec(:,4)));...
    u - 10*eta*log10(eucl([x(1,:); x(4,:)], pos_vec(:,5)));...
    u - 10*eta*log10(eucl([x(1,:); x(4,:)], pos_vec(:,6)))]', zeta^2*eye(6,6));

saved_w = zeros(N,n);

w = p(part, Y(:,1));
saved_w(:,1) = w;

tau(:,1) = [sum(part(1,:) .* w')./sum(w);...
            sum(part(2,:) .* w')./sum(w);...
            sum(part(3,:) .* w')./sum(w);...
            sum(part(4,:) .* w')./sum(w);...
            sum(part(5,:) .* w')./sum(w);...
            sum(part(6,:) .* w')./sum(w)];

for t = 2:n
    
    Z = Z_space(:,randsample(5, N,true));

    if resampling == true
        ind = randsample(N,N, true, w);
        part = phi * part(:,ind) + ksi_z * Z + ksi_w * randn(2,N)*sigma;
        l = log(p(part, Y(:,t)));
    else
        part = phi * part+ ksi_z * Z + ksi_w * randn(2,N)*sigma;
        l = log(p(part, Y(:,t))) + log(w);
    end

    w = exp(l - max(l));
    saved_w(:,t) = w;
    tau(:,t) = [sum(part(1,:) .* w')./sum(w);...
            sum(part(2,:) .* w')./sum(w);...
            sum(part(3,:) .* w')./sum(w);...
            sum(part(4,:) .* w')./sum(w);...
            sum(part(5,:) .* w')./sum(w);...
            sum(part(6,:) .* w')./sum(w)];
end


%% PLOT
figure
plot(pos_vec(1,:), pos_vec(2,:), 'o', LineWidth = 4)
hold on
plot(tau(1,:), tau(4,:), LineWidth=2)


figure
for i = 1:5
    indexes = [2; 10; 20; 30; 200];
    histogram(saved_w(:,indexes(i)));
    legend(["n=2" "n=10" "n=20" "n=30" "n=200"])
    ylim([0 10000])
    hold on
end
figure
ESS_saved = zeros(n,1);
for i = 1:n
samplesize = 1/sum((saved_w(:,i)/sum(saved_w(:,i))).^2);
ESS_saved(i) = samplesize;
end
plot(ESS_saved)
title("Efficient sample size for different weights")
ylabel("ESS")
xlabel("Time steps (n)")