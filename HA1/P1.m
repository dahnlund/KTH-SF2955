%% Problem 1

dt = 0.5;
a = 0.6;
sigma = 0.5;
Z_space = [0 3.5 0 0 -3.5; 0 0 3.5 -3.5 0];

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

m = 500;
X = zeros(6,m);
X(:,1) = diag([500 5 5 200 5 5]).^(1/2) * randn(6, 1);

Z_indices = randsample(5,1,true); % Z_initialization
Z = Z_space(:, Z_indices);
for n = 2:m
    X(:,n) = phi * X(:,n-1) + ksi_z * Z + ksi_w * randn(2,1)*sigma;
    Z_indices = P_sample(Z_indices);
    Z = Z_space(:,Z_indices);
end
plot(X(1,:), X(4,:), LineWidth = 2)
title("Generated trajectory")
xlabel("X_1")
ylabel("X_2")


