%{ 
   Valuing European Put options using the Least Squares Monte Carlo
   Determining the expected exposure and credit value adjustment (CVA)
   (longstaff & schwartz algorithm, 2001)
   @ Aniruddha Dutta 
   ---------------------------------------------------------------------
   Europeans Put options valued using the Least Squares Monte Carlo and
   Black Scholes formula. The results show that LSM method is not a 
   good way to price european options at each time steps which leads to 
   incorrect expected exposure and CVA values.
  %}

%{ -------- Parameters -------------%}


T = 2.0;            % Time to maturity
r = 0.06;           % Risk free interest rate
sigma = 0.2;        % Volatility
K = 40;             % Strike price
S0 = 36;            % Underlying asset price in-the-money
%S0 = 44;           % Underlying asset price out-of-money
N = 50;             % number of time steps
M = 100000;         % 100k price paths
k = 3;              % Number of Basis Functions used for Laguerre Polyomials

dt = T/N;           % value of time step
t = 0:dt:T;         % time vector

% Generate stock price matrix (100k price paths) 

S = ones(M, N+1);    
R = exp((r - sigma^2/2)*dt+sigma*sqrt(dt)*randn(N,M));
SS = cumprod([S0*ones(1,M); R]);
S = SS';

%xlswrite('myFile_S044.xlsx',S)   % Write stock paths in a excel file
%filename = 'myFile_S036.csv';   % filename for S0 equals 36
%filename = 'myFile_S044.csv';   % filename for S0 equals 44
%S = csvread(filename);          % Read the file
%size(S)                         % size of stock paths

P = max((K - S(:,N+1)),0);      %payoff at time T

X1 = []; 
C1 = []; 
P1 = [];
S1 = [];

for i = N:-1:2
   
    itmP = find(S(:,i));        % All paths
    
    X =  S(itmP,i);             % All price paths
    
    Y = P(itmP)*exp(-r*dt);     % discounted payoffs from continuation
    
    A = BasisFunctions(X,k);
    beta = A\Y;                 % regression coeff
    
    C = A * beta;               % Estimated continuation value
    E = max(K - X,0);           %  Immediate exercise value
    
    exP = itmP(K - X>0);        % Paths where strike greater than spot
    
    %rest = setdiff(1:M,exP);   % Rest of the paths
    
    %P(exP) = E(C<E);   % Better to exrcise? Insert value in payoff vector
    
    P(exP) = P(exP)*exp(-r*dt); %Insert payoffs and discount back one step
    
    u = mean(P * exp(-r*dt));   % Value of the option
    opt_val(i) = u;
    
    X1 = [X1 X]; 
    C1 = [C1 C]; 
    P1 = [P1 P];
    
end

% Value of the option from LSM
fprintf('European Put Option value from LSM Regression is $%.3f\n\n', u);   

% Calculate the initial option value from Black Scholes Method

[Call,Put] = blsprice(S0,K,r,T,sigma);    % initial option value
fprintf('Black Scholes European Put Option value is $%.3f\n\n', Put);

init_opt_val = Put;   

% X1_val, C1_val, P1_val, z1_val would collectively store specific.... 
....distributed values of X1,C1,P1 and z1 for all i - 1 to 49
X1_val = []; 
C1_val = [];
P1_val = [];
z1_val = [];


for i = 1:1:49
    
    % X_val,C_val,P_val stores distributed values for individual X vectors.. 
    % ... for a specific set of spot prices
    
    X_val = [];
    C_val = [];
    P_val = [];
    
    %z1 and z2 are the range of z values for +, -  1.96 standard...
    ....deviations of the underlying S0
    
    sp = 14;                         % picking 15 values of spot price
    z1 = S0 - 1.96 * std(S(:,i)) * sqrt(T - i * dt);     
    z2 = S0 + 1.96 * std(S(:,i)) * sqrt(T - i * dt);
    z = z1:(z2 - z1)/sp:z2;                          
    
    for j = 1:length(z)
    
        X_val = [X_val ; X1(j,i)];
        C_val = [C_val ; C1(j,i)];
        P_val = [P_val ; P1(j,i)];
        
        
    end
    
    X1_val = [X1_val X_val];
    C1_val = [C1_val C_val];
    P1_val = [P1_val P_val];
    z1_val = [z1_val z'];
end
  

 %{Generating regression values vs spot prices plot for each specific...
...time step after sorting.Time step ranges from 1-49.

X1_val_sort=sort(X1_val(:,45),'ascend'); % X1_val spot price at time step 25
C1_val_sort=sort(C1_val(:,45),'descend');% C1_val reg values at time step 25
title('European Put Option vs Spot Price','Fontsize',15,....
    'Fontweight','bold','Color','k');
%scatter(X1_val(:,45),C1_val(:,45),35,'b');    % scatter plot if needed

plot(X1_val_sort,C1_val_sort,'-o');       
xt = get(gca, 'XTick');                   
yt = get(gca, 'YTick');
set(gca, 'FontSize', 14)
hold on
%}
%{ 
The following code will generate plots for black 
scholes values for different spot 
spot prices at individual specific time steps. So we can have plots for 
each mentioned time step.The following is an example for time step 25.
%}


Put1 = [];        % vector stores all the BS option val at diff time steps
Put_Col = [];
[Call,Put] = blsprice(X1_val(:,45),40,.06,(T - 45 * dt),.2);   % BS value at time step 25
Put_sort = sort(Put,'descend');            % Sorting the values
plot(X1_val_sort, Put_sort,'-*');   % BS option val vs spot price for time step 25
xlabel('Spot Price','FontSize',16,'FontWeight','bold','Color','k')
ylabel('Option Value','FontSize',16,'FontWeight','bold','Color','k')
hold off
legend('Regression Function','Black Scholes Value');
diff = abs(Put_sort - C1_val_sort);        % absolute error between black 
                                           % scholes and regression values
err = sum(diff);                           % total error

%}

%{-------------------------------------------------------------------------------------------------------------%}
% Black Scholes price calculated for each time step

Put1 = [];                  % vector to store Put val at all time steps
Put_Col = [];

for iCol = 1:1:48           % all columns of X1 
    Put_Col = [];           % init a temp col to store intermediate results
    for iRow = 1:1:15     % 1000 rows of X1

        %X1(iRow, iCol) takes all values of X1 
        %iCol *dt updates t for each colum of X1. for column1 T-dt, for
        %col2 T-2dt ... and so on 
        [Call,Put] = blsprice(X1_val(iRow,iCol),40,.06,T-iCol*dt,.2);    
        Put_Col = [Put_Col; Put]; 
        
    end 
    Put1 = [Put1 Put_Col]; 
end

%{
The loop below is for caculating the collateral at each step which is 
set at different values.
%}
%  

for i = 2:5:49
    %coll = mean(Put1(:,i-1));    % collateral = option value previous step
    coll = init_opt_val;          % collateral = initial option value
    %coll = 2 * init_opt_val;     % collateral = twice option value
    exp_reg = max(C1  - coll,0);  % exposure from LSM regression values
    exp_bls = max(Put1 - coll, 0);% exposure from black scholes values
    
end

% Find expected exposure from exp_bls and exp_reg for all time steps

expt_reg = mean(exp_reg);   % expected exposure from LSM regression values
pfe_reg = prctile(exp_reg, 95);  % PFE from regression
expt_bls = mean(exp_bls);   % expected exposure from black scholes values

tme = [];

for i = 1:1:49
    time = T - i*dt;
    tme = [tme;time];
    pfe_reg = prctile(exp_reg, 95);
    pfe_bls = prctile(exp_bls, 95);

end


% Plotting the expected exposure from regression and black scholes with
% time

plot(tme(1:3:49),expt_bls(1:3:49),'-o');
title('European Put Expected Exposure vs Time','Fontsize',12,....
'Fontweight','bold','Color','k');
xt = get(gca, 'XTick');
yt = get(gca, 'YTick');
set(gca, 'FontSize', 12)
xlabel('Time','FontSize',12,'FontWeight','bold','Color','k')
ylabel('Expected Exposure','FontSize',12,'FontWeight','bold','Color','k')
%legend('Black Scholes','Regression');
hold on
plot(tme(1:3:49),expt_reg(1:3:49),'-o');
plot(tme(1:3:49),pfe_reg(1:3:49), '-o');
plot(tme(1:3:49),pfe_bls(1:3:49), '-o');
legend('Black Scholes','Regression','PFE Regression','PFE Black Scholes');
hold on


    
%{ ---------------- end of program --------------------------------- end of program --------------------------------- %}
