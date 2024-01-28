clear variables; clc

load("stations.mat"); load("RSSI-measurements-unknown-sigma.mat")

%% DEFINE PROBLEM

resampling = true;  % If true, perform resampling

%% Define parameters

dt = 0.5;
a = 0.6;
sigma = 0.5;
Z_space = [0 3.5 0 0 -3.5; 0 0 3.5 -3.5 0];
P = 1/20 * (ones(5,5) + 15*eye(5,5));  % This is not used in the iterations, since the chain already is stationary in a uniform distribution

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

%% SIS
N = 10000;
n = length(Y(1,:));

u = 90;
eta = 3;

liklist = [];
likl = [];
for zeta =  1.5:0.05:2.7

    part = diag([500 5 5 200 5 5]).^(1/2) * randn(6, N); % Initialization
    
    tau = zeros(6,n);
    
    p =  @(x,y) mvnpdf(y', [u - 10*eta*log10(eucl([x(1,:); x(4,:)], pos_vec(:,1))); ...
        u - 10*eta*log10(eucl([x(1,:); x(4,:)], pos_vec(:,2)));...
        u - 10*eta*log10(eucl([x(1,:); x(4,:)], pos_vec(:,3)));...
        u - 10*eta*log10(eucl([x(1,:); x(4,:)], pos_vec(:,4)));...
        u - 10*eta*log10(eucl([x(1,:); x(4,:)], pos_vec(:,5)));...
        u - 10*eta*log10(eucl([x(1,:); x(4,:)], pos_vec(:,6)))]', zeta^2*eye(6,6));
    
    saved_w = zeros(N,n);
    
    w = p(part, Y(:,1));
    saved_w(:,1) = w;
    
    logOmega = log(sum(saved_w(:,1)));

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
        logOmega = logOmega + (sum(l));
    end

    disp(logOmega)
    likl = [likl; logOmega];
end
likl = 1/n*(-(n+1)*log(N)+likl);

zetas = 1.5:0.05:2.7;
plot(zetas,likl)
hold on
plot(zetas(likl == max(likl)),max(likl),'*')
figure
disp("Zeta_opt =")
disp(zetas(likl == max(likl)))
%% PLOT
plot(pos_vec(1,:), pos_vec(2,:), 'o', LineWidth = 4)
hold on
plot(tau(1,:), tau(4,:), LineWidth=2)