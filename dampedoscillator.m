% EKF application to parameter estimation for dampled harmonic oscilltor 

% Here we make use of jacobian of the system matrix as we here linearize
% the system 

% In UKF we don't linearize the system matrix instead we make use of true
% dynamics of the system 


clc ; close all ; clear all ; 

m = 10 ; k = 5 ; b = 3 ; 

omega = sqrt(k/m) ; 

bspring = b/m ; 

A = [ 0 1 ; -k/m -b/m] ; % Linearized system matrix 

rhs = @(t,x) A*x ; 

% Simulate system 

xinit =  [0;1] ; h = 0.01 ; T = 50 ; 

time = 0:h:T ; 

[t, trueTrajectory ] = ode45(rhs , time , xinit ) ; 

obsNoise = 0.1^2 ; 

obs = trueTrajectory(:,1) ; % system output thing 1st state as a output we have 

obs = obs + obsNoise*randn(size(obs)) ; 


omega = k/m ; 
gamma = b/m ; 



%EKF 

xbar = [2; 3 ;6; 8] ;

xest = xbar ; 

P = diag([0.1 0.1 0.1 0.1]) ; 

varest = diag(P) ; 

A = [0 1 0 0 ; -xbar(3)/m -xbar(4)/m -xbar(1)/m -xbar(2)/m ; 0 0 0 0 ; 0 0 0 0];

% F = eye(length(xbar)) + A*h + (A^2*h^2)/factorial(2) 

F = eye(length(xbar)) + A*h + (A^2*h^2)/factorial(2) + (A^3*h^3)/factorial(3) + (A^4*h^4)/factorial(4) + (A^5*h^5)/factorial(5) + (A^6*h^6)/factorial(6)
% F is state transition matrix 

H = [1 0 0 0 ] ; 

Q_variance = 0.01^2 ; 

Q = piecewise_white_noise(4,Q_variance , h ) ; 

R = obsNoise ; 

for i = 2:length(obs) 

    % Predict step 
    f = { @(x1 , x2 , x3 , x4 , t ) (x1+x2*t) ; 
          @(x1 , x2 , x3 , x4 , t ) (x2 + (-x4/m*x2 - x3/m*x1)*t)  ;
          @(x1 , x2 , x3 , x4 , t ) (x3)  ; 
          @(x1 , x2 , x3 , x4 , t ) (x4)} ;

    xbar(1,1) = f{1}(xbar(1) , xbar(2) , xbar(3) , xbar(4) , h) ; 
    xbar(2,1) = f{2}(xbar(1) , xbar(2) , xbar(3) , xbar(4) , h) ; 
    xbar(3,1) = f{3}(xbar(1) , xbar(2) , xbar(3) , xbar(4) , h) ; 
    xbar(4,1) = f{4}(xbar(1) , xbar(2) , xbar(3) , xbar(4) , h) ;

    A = [0 1 0 0 ; -xbar(3)/m -xbar(4)/m -xbar(1)/m -xbar(2)/m ; 0 0 0 0 ; 0 0 0 0];

    F = eye(length(xbar)) + A*h + (A^2*h^2)/factorial(2) + (A^3*h^3)/factorial(3) + (A^4*h^4)/factorial(4) + (A^5*h^5)/factorial(5) + (A^6*h^6)/factorial(6)

    P = F*P*F' + Q ; 

    % Correction step 

    K = P*H'*inv(H*P*H' + R ) ; 
    y = obs(i) - H*xbar ; 

    xbar = xbar + K*(y) ; 
    P = P - K*H*P ; 

    xbarEstimate(:,i) = xbar(:,1) ; 
    varEstimate(:,i) = diag(P) ; 
    KalmanGain(:,i) = K ; 
    Residual(:,i) = y ; 

end


% Position and velocity plots 

% figure('Renderer','painters','Position',[15 15 900 700]) 

figure(1) 

plot(time , trueTrajectory(:,1) , 'kx' , 'LineWidth',1) ; 
hold on ; 

plot(time , obs , 'bo' , 'MarkerFaceColor','b' , 'MarkerSize',1.5) ;

plot(time , xbarEstimate(1,:) , 'r' , 'LineWidth', 1.5) ; 

legend('TrueTrajectory' , 'noisy output of position' , ' EKF estimate ') ; 
title('Estimation of position ') ; 
xlabel('Time') ; 
ylabel('Position') ;
% lower_x1 = xbarEstimate(1,:) + 2*sqrt(varEstimate(1,:)) ; 
% upper_x1 = xbarEstimate(1,:) - 2*sqrt(varEstimate(1,:)) ; 
% 
% ciplot(lower_x1 , upper_x1 ,'b')

figure(2) 
plot(time , trueTrajectory(:,2) , 'k' , 'LineWidth',2) ; 
hold on ; 

% plot(time , obs , 'bo' , 'MarkerFaceColor','b' , 'MarkerSize',2) ;
%Cause we are measuring only one state 

plot(time , xbarEstimate(2,:) , 'r' , 'LineWidth',2); 

legend('TrueTrajectory' , ' EKF estimate '); 
title('Estimation of velocity '); 
xlabel('Time'); 
ylabel('velocity');


figure(3) 
% plot(time , trueTrajectory(:,3) , 'k' , 'LineWidth',1.5) ; 
% hold on ; 

% plot(time , obs , 'bo' , 'MarkerFaceColor','b' , 'MarkerSize',2) ;
%Cause we are measuring only one state 

plot(time , xbarEstimate(3,:) , 'r' , 'LineWidth',2) ; 
title('Spring constant ') ; 
xlabel('Time') ; 
ylabel('value of k ') ;
figure(4) 

% 
% plot(time , trueTrajectory(:,4) , 'k' , 'LineWidth',1.5) ; 
% hold on ; 

% plot(time , obs , 'bo' , 'MarkerFaceColor','b' , 'MarkerSize',2) ;
%Cause we are measuring only one state 

plot(time , xbarEstimate(4,:) , 'r' , 'LineWidth',2) ; 

title('Estimation of damping factor') ; 
xlabel('Time') ; 
ylabel('damping factor') ;




