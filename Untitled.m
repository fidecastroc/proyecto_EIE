clear, close all,clc
%

I = 3, 8, 26, 76, 225, 298, 258, 233, 189, 128, 68, 29, 14, 4;
n = 0:length(I)-1;
Ii        = 1; Ri        = 0;N        = 763;
y0        = [N-(Ii+Ri), Ii, Ri]; %Initial conditions (SIR)
%LSE
[cc,fval] = fminsearch(@(c) sir_ssq(c, N, y0, I),[0.5 0.5]);
beta      = cc(1); g =cc(2);
%ML
[cc_ml,fval] = fminsearch(@(c) sir_mle(c, N, y0, I),[0.5 0.5]);
beta_ml      = cc_ml(1); g_ml =cc_ml(2);
%Final Solution LS & ML
tspan     = linspace(0,14,100);
[t,y_ls]  = ode45(@(t,y) sir(t, y, beta/N, g),tspan,y0); 
[t,y_ml]  = ode45(@(t,y) sir(t, y, beta_ml/N, g_ml),tspan,y0); 
%plot
plot(n,I,'b*', 'MarkerSize', 12)
title('Influenza in a boarding school (British Medical Journal, 4 March 1978)')
hold on
plot(t,y_ls(:,2),'r','LineWidth', 4.0,...
     t,y_ml(:,2),'g','LineWidth', 4.0)
xlabel('Time/days since first case')
ylabel('Number of confirmed people per day')
grid
h = legend ('Confirmed', strcat('LS :',num2str(cc)), strcat('ML :',num2str(cc_ml)));
legend (h, 'location', 'northwest');

function dydt=sir(t,y,beta,g)
  S        = y(1);  I = y(2);  R= y(3); %Initial conditions
%Diferences  
  dS       = -beta*I.*S;
  dI       =  beta*I.*S  - g*I;
  dR       =  g*I;
  dydt     = [dS;dI;dR;];
end
%
function f = sir_ssq(coeff, N, y0, I)
    beta   = coeff(1);  g = coeff(2);%Initial Free parameters
    tspan  = 0:length(I)-1;
    [t,y]  = ode45(@(t,y) sir(t, y, beta/N, g),tspan,y0); 
    Y_fun  = y(:,2)';
    f      = sum((Y_fun-I).^2);
end
function f =sir_mle(coeff, N, y0, I)
    beta   = coeff(1);  g = coeff(2); %Initial Free parameters
    tspan  = 0:length(I)-1;
    [t,y]  = ode45(@(t,y) sir(t, y, beta/N, g),tspan,y0); 
    Y_fun  = y(:,2)';
    f      = -sum(log(poisspdf(I,mean(Y_fun))));
end   

