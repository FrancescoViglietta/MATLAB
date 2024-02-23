%% Likelihood for heteroskedastic SVAR
%author: Francesco Viglietta
function [logLik]= Likelihood_Fanelli_Heterosk(teta)
%the "teta" will be the input value of this function, which in the
%minimization command will be the randomly generated initial values
global Sigma_u_1
global Sigma_u_2
global T_2
global T_1
global A
global C
global M
global ParamNumberA
global ParamNumberC

%Declaring the structure of the matrix A




for c_par = 1 : size(ParamNumberA,1)
    A(ParamNumberA(c_par,1))=teta(c_par);   
end

     



    for c_par = 1 : size(ParamNumberC,1)
    C(ParamNumberC(c_par,1))=teta(c_par+ size(ParamNumberA,1));   %this for loop puts correctly the imput values in the matrix A according to the positions declared in ParamNumberB
    end

    SIGMA_1=inv(A)*inv(A');
    SIGMA_2=inv(A+C)*inv(A+C)'; 
    SIGMA_1_inv=inv(SIGMA_1);
    SIGMA_2_inv=inv(SIGMA_2);

    
logLik = +0.5*(log(2*pi)*(T_1*M+T_2*M)) +0.5*T_1*(log(det(SIGMA_1))+trace(Sigma_u_1*SIGMA_1_inv  ) ) +0.5*T_2*(log(det(SIGMA_2))+trace(Sigma_u_2*SIGMA_2_inv  ) )  ;                                                                              

%the loglikelihood formula usually used to estimate the variance-covariance matrix after having estimated the reduced form parameters 
                                                                                                     
%+log(2*pi)*(T_1*M+T_2*M) to chek if it's correct to add this term. However
%it's a constant and  doesn't affect the maximization

% putting the '+' sign in front of each term because the fminun command
% minimizes. The original log likelihood has the '-' sign in front of each
% term
    
 end