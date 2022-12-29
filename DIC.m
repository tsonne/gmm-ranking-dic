% Ranking Ground Motion Models (GMMs) based on the Deviance Information
% Criterion (DIC) by Milad Kowsari

% Inputs: A matrix that includes ground motion observations (the first
% column) and ground motion predictions from different GMMs (the second column until the last one)

% Outputs: 
% Two cases are assumed for the modeling: 
% (1) when the standard deviation of the model is assumed to be known (denoted as DIC1). Only the original sigma is used. 
% (2) when the standard deviation of the model is unknown (denoted as
% DIC2). The posterior sigam will be estimated by Bayesian statistics.

% The models with smaller DIC scores perform better compared to models with higher scores.


% Reference:
% Kowsari, M., Halldorsson, B., Hrafnkelsson, B., & Jónsson, S. (2019).
% Selection of earthquake ground motion models using the deviance information criterion.
% Soil Dynamics and Earthquake Engineering, 117, 288-299.

%% Inputs
DATA=load('inputs.txt'); % 
obs=DATA(:,1); % Natural logarithm of observed data
N=length(obs); % Number of observations
pred=DATA(:,2:end); % Predictions from the selected GMMs

ss=0.6; % Scaling parameter; suggested values: 0.6 if the inputs are in natural logarithm; 0.3 if the inputs are in base 10 logarithm
ss2=ss.^2; 
nu0=1; % the number of degrees of freedom for the chi-squared distribution
nu=N+nu0;% the number of degrees of freedom for the inverse chi-squared distribution

org_sigma=[0.665 0.642 0.66	0.723	0.667	0.561	0.524	0.792]; % sigma of the selected GMMs in natural logarithm

%% Case 1. if sigma is fixed (original sigma) then D_theta_hat=D_theta_ave;
  for i=1:length(org_sigma)
      
D_theta_hat1(:,i)=N*log(2*pi)+N*log((org_sigma(:,i)).^2)+(1./(org_sigma(:,i)).^2).*sum((obs-pred(:,i)).^2); % The deviance of normal (D_theta_hat)

  end
  
D_theta_ave1=D_theta_hat1;

DIC1=2*D_theta_ave1-D_theta_hat1

%% Case2. if sigma is unknown 
  for i=1:length(org_sigma)
   tau2(:,i)=(nu^-1)*((nu0*ss2)+sum((obs-pred(:,i)).^2)); % scaling parameter for the inverse chi-squared distribution
   y(:,i) =(tau2(:,i)*nu)./sichi2rnd(nu,10000,1); % can be also obtained by MCMC 
   y(:,i) =sqrt(y(:,i));
  
mean_y(:,i)=sqrt((nu*tau2(:,i))/(nu-2));

D_theta_ave(:,i)=mean(N*(log(2*pi))+N*log(y(:,i).^2)+(1./(y(:,i).^2)).*sum((obs-pred(:,i)).^2)); % D_theta_average

D_theta_hat(:,i)=N*log(2*pi)+N*log(mean_y(:,i).^2)+(1./(mean_y(:,i).^2)).*sum((obs-pred(:,i)).^2); % The deviance of normal (D_theta_hat)


  end
  
  
DIC2=2*D_theta_ave-D_theta_hat

