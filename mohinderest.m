% EKF estimation of states ( position and velocity ) and parameter ( zeta )
% Mohinder grewal book 

% Here we make use of jacobian of the system matrix as we here linearize
% the system 

% In UKF we don't linearize the system matrix instead we make use of true
% dynamics of the system 

clc ; close all ; clear all ; 

% m = 10 ; k = 11 ; b = 5 ; 

% omega = sqrt(k/m) ; 

omega = 10 ; 

zeta = 0.1 ; 

% bspring = b/m ; 

A = [ 0 1 ; -omega^2 -2*zeta*omega] ; % Linearized system matrix 

rhs = @(t,x) A*x ; 

% Simulate system 

xinit =  [0;1] ; h = .02 ; T = 20 ; 

time = 0:h:T ; 

[t, trueTrajectory ] = ode45(rhs , time , xinit ) ; 

obsNoise = 0.001^2 ; 

obs = trueTrajectory(:,1) ; % system output thing 1st state as a output we have 

obs = obs + obsNoise*randn(size(obs)) ; 


% omega = k/m ; 
% gamma = b/m ; 



%EKF 

xbar = [ 0; 1 ; 0.7] ;

xest = xbar ; 

P = diag([2 2 2]) ; 

varest = diag(P) ; 

A = [0 1 0 ; -omega^2 -2*xbar(3)*omega -2*xbar(2)*omega ; 0 0 0 ]  ; 

% F = eye(length(xbar)) + A*h + (A^2*h^2)/factorial(2) 

F = eye(length(xbar)) + A*h + (A^2*h^2)/factorial(2) + (A^3*h^3)/factorial(3) + (A^4*h^4)/factorial(4) + (A^5*h^5)/factorial(5) + (A^6*h^6)/factorial(6)
% F is state transition matrix 

H = [1 0 0 ] ; 

Q_variance = sqrt(4.47) ; 

Q = piecewise_white_noise(3,Q_variance , h ) ; 

R = obsNoise ; 

for i = 2:length(obs) 

    % Predict step 
    f = { @(x1 , x2 , x3 , t ) (x1+x2*t) ; 
          @(x1 , x2 , x3  , t ) (-(omega^2)*x1*t + (1 - 10*t*x3)*x2  )  ;
          @(x1 , x2 , x3  , t ) (x3) } ;

    xbar(1,1) = f{1}(xbar(1) , xbar(2) , xbar(3) , h) ; 
    xbar(2,1) = f{2}(xbar(1) , xbar(2) , xbar(3) , h) ; 
    xbar(3,1) = f{3}(xbar(1) , xbar(2) , xbar(3) , h) ; 
%     xbar(4,1) = f{4}(xbar(1) , xbar(2) , xbar(3) , xbar(4) , h) ;

    A = [0 1 0 ; -omega^2 -2*xbar(3)*omega -2*xbar(2)*omega ; 0 0 0 ]  ; 

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

plot(time , trueTrajectory(:,1) , 'kx' , 'LineWidth',2) ; 
hold on ; 

plot(time , obs , 'bo' , 'MarkerFaceColor','b' , 'MarkerSize',3) ;

plot(time , xbarEstimate(1,:) , 'r' , 'LineWidth',2) ; 
legend('TrueTrajectory' , 'noisy output of position' , ' EKF estimate ') ; 
title('Estimation of position ') ; 
xlabel('Time') ; 
ylabel('Position') ;
% plot(xbarEstimate(1,:)' - trueTrajectory(:,1) , 'y')  % error is of order 10^(-3)

% lower_x1 = xbarEstimate(1,:) + 2*sqrt(varEstimate(1,:)) ; 
% upper_x1 = xbarEstimate(1,:) - 2*sqrt(varEstimate(1,:)) ; 
% 
% ciplot(lower_x1 , upper_x1 ,'b')

figure(2) 
plot(time , trueTrajectory(:,2) , 'k' , 'LineWidth',1.5) ; 
hold on ; 

% plot(time , obs , 'bo' , 'MarkerFaceColor','b' , 'MarkerSize',2) ;
%Cause we are measuring only one state 

plot(time , xbarEstimate(2,:) , 'r' , 'LineWidth',2) ;

% plot(xbarEstimate(2,:)' - trueTrajectory(:,2) , 'y')  %0.07 error 
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
title('Estimation of damping factor') ; 
xlabel('Time') ; 
ylabel('damping factor') ;



% figure(4) 
% 
% % 
% % plot(time , trueTrajectory(:,4) , 'k' , 'LineWidth',1.5) ; 
% % hold on ; 
% 
% % plot(time , obs , 'bo' , 'MarkerFaceColor','b' , 'MarkerSize',2) ;
% %Cause we are measuring only one state 
% 
% plot(time , xbarEstimate(4,:) , 'r' , 'LineWidth',2) ; 






