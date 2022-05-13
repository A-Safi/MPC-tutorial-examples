yalmip('clear')
clear

k=1; % spring stiffness
m=10; % Masses
c=0.1; % Damping
t_0 = 0 ; % Start Time
t_f = 100; % Final Time
Ts = 0.5; % Time Step
t = t_0:Ts:t_f ; % Time Vector
Nt=length(t);

%% State-Space Model
A=[0 0 1 0;
0 0 0 1
-k/m k/m -c/m 0
k/m -k/m 0 -c/m];
B=[0 0;0 0;1/m 0;0 1/m];
C=[1 0 0 0;0 1 0 0];
tran=ss(A,B,C,0);

%% State-Space Model in discrete-time
trand = c2d(tran,Ts);
Am = trand.a;
Bm = trand.b;
Cm = trand.c;
Dm = trand.d;

nx = size(Am,1); % Number of states
nu = size(Bm,2); % Number of inputs
ny = size(Cm,1);

%% Construction of MPC matrices
R = 1e-3*eye(nu);
Q = eye(ny);

N = 10; % Prediction Horizon

%% Prediction variables
u = sdpvar(repmat(nu,1,N),repmat(1,1,N));
x = sdpvar(repmat(nx,1,N+1),repmat(1,1,N+1));
r = sdpvar(repmat(ny,1,N+1),repmat(1,1,N+1));

%% Optimization Problem
objective = 0;
constraints = [];
for k = 1:N
    objective = objective + (C*x{k}-r{k})'*Q*(C*x{k}-r{k}) + u{k}'*R*u{k};
    constraints = [constraints, x{k+1} == Am*x{k}+Bm*u{k}];
    constraints = [constraints, -[5;5] <= u{k} <= [5;5],...
                            -[6;6;6;6] <= x{k+1} <= [6;6;6;6]];
end
objective = objective + (C*x{N+1}-r{N+1})'*(C*x{N+1}-r{N+1});

parameters_in = {x{1},[r{:}]};
solutions_out = {[u{:}], [x{:}]};

controller = optimizer(constraints, objective,[],parameters_in,solutions_out);

%% Buffer
uu = zeros(nu , Nt) ;
xx = zeros(size(Am , 1) , Nt) ;

%% Initial conditions
xx(:,1) = [0;0;0;0];
y(:,1) = Cm*xx(:,1);

%% Refrence trajectory
Rs(1 , :) =[ones(floor(Nt/4) , 1)' , 2*ones(floor(Nt/4) , 1)' , -ones(floor(Nt/4) , 1)' , zeros(floor(Nt/4+1) , 1)'] ;
Rs(2 , :) =-[ones(floor(Nt/4) , 1)' , 2*ones(floor(Nt/4) , 1)' , -ones(floor(Nt/4) , 1)' , zeros(floor(Nt/4+1) , 1)'] ;

%% Simulation Loop
for i = 1:Nt-N-1
    future_r = Rs(: , i:i+N);    
    inputs = {xx(:,i),future_r};
    [solutions,~] = controller{inputs};    
    U = solutions{1};
    uu(:,i)=U(:,1);
    xx(:,i+1) = Am*xx(:,i) + Bm*uu(:,i);
    y(:,i+1) = Cm*xx(:,i+1);
end

%%  Plot Results
figure(1);
subplot(2,1,1)
plot(t(1:Nt-N-1) , Rs(1 , 1:Nt-N-1),'--b', 'LineWidth' , 1.5); hold on
plot(t(1:Nt-N-1) , y(1 , 1:Nt-N-1) ,'-r', 'LineWidth' , 2)
plot(t(1:Nt-N-1) , Rs(2 , 1:Nt-N-1),'--b', 'LineWidth' , 1.5); hold on
plot(t(1:Nt-N-1) , y(2 , 1:Nt-N-1) ,'-r', 'LineWidth' , 2) 
xlabel('Time  (second)') ;
ylabel('Output') ;
title('PFC control - MIMO') ;
grid on

subplot(2,1,2)
plot( t(1:Nt-N-1) , uu(1 , 1:Nt-N-1) ,'-r', 'LineWidth' , 2) ; hold on
plot( t(1:Nt-N-1) , uu(2 , 1:Nt-N-1) ,'-b', 'LineWidth' , 2) ; hold on

xlabel('Time  (second)') ;
ylabel('Input') ;
title('PFC control - MIMO') ;
grid on