clc;
clear;
% close all;
%%
k=1; % spring stiffness
m=1; % Mass
c=1; % Damping
t_0 = 0; % Start Time 
t_f = 10; % Final Time
Ts = 0.05; % Time Step
t = t_0:Ts:t_f ; % Time Vector

%% State-Space Model
A=[0 1;-k/m -c/m];  
B=[0;1/m];
C=[1 0];
tran=ss(A,B,C,0);

%% State-Space Model in discrete-time
trand = c2d(tran,Ts);
Am = trand.a;
Bm = trand.b;
Cm = trand.c;
Dm = trand.d;
n = size(Am ,1)  ;  % number of states
q = size(Cm , 1) ;  % number of outputs

%% Construction of MPC matrices
N = 10 ; % Prediction Horizon
G = zeros(q*N , size(Am , 1)) ;
for i = 1:N
    G(q*i-q+1:q*i , :) = Cm * Am^i ;
end
for i = 1:N
    for j = 1:i     
        H(q*i-q+1:q*i , j) = Cm * Am^(i-j) * Bm ;
    end
end
H=[zeros(q,N);H(1:q*(N-1),1:N)];
F = zeros(N*q , size(Bm , 2)) ;
for i = 1:N
    F(q*i-q+1:q*i , :) = Cm * Am^(i-1) * Bm ;
end
DiffM=eye(N);
for i=1:N-1
DiffM(i,i+1)=-1;
end
DiffM(end)=0;
r=0.001;
Q=r*DiffM'*DiffM+H'*H;

%% Refrence trajectory
Nt = numel(t) ;
Rs(1 , :) =[ones(floor(Nt/4) , 1)' , 2*ones(floor(Nt/4) , 1)' , -ones(floor(Nt/4) , 1)' , zeros(floor(Nt/4+1) , 1)'] ;
 
%% Prediction variables
y = zeros(q , Nt) ;
u = zeros(1 , Nt) ;
x = zeros(size(Am , 1) , Nt) ;

%% Initial conditions
x(:,1) = [1 ; 0];
y(:,1) = Cm*x(:,1);

%% Optimization inequality matrices
Aineq=[eye(N);-eye(N)];
Bineq=5*ones(2*N,1);
options = optimoptions('quadprog','Display','off'); % Solver

%% Simulation loop
for  i = 1:Nt-N-1
    x(:,i+1)=Am*x(:,i)+Bm*u(:,i);
    y(:,i+1)=Cm*x(:,i);

    %% Calculating input using analytical solutions
    % (comment this or the next section)
    f=H'*(G*x(:,i)+F*u(:,i)-Rs(:,i:i+N-1)');
    U=-Q\f;

    %% Calculating input using quadprog
    %U=quadprog(Q,f',Aineq,Bineq,[],[],[],[],[],options);

    %% Update input for next iteration
    u(i+1)=U(1,:);
end

%%  Plot Results
figure(1);
subplot(2,1,1)
plot( t(1:Nt-N-1) , Rs(1 , 1:Nt-N-1) ,'--b', 'LineWidth' , 2) ; hold on
plot( t(1:Nt-N-1) , y(1 , 1:Nt-N-1) ,'-g', 'LineWidth' , 2) ; hold on
xlabel('Time  (second)') ;
ylabel('Amp - Y1') ;
title('PFC control - SISO') ;
grid on
legend('Ref','Out1') ;
subplot(2,1,2)
plot( t(1:Nt-N-1) , u(1 , 1:Nt-N-1) ,'-g', 'LineWidth' , 2) ; hold on
xlabel('Time  (second)') ;
ylabel('Amp - u') ;
title('PFC control - SISO') ;
grid on
legend('Control Effort') ;