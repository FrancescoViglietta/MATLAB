%% Likelihood for heteroskedastic SVAR
%author: Francesco Viglietta
function [logLik]= Likelihood_Fanelli_Heterosk_boot(teta)
%the "teta" will be the input value of this function, which in the
%minimization command will be the randomly generated initial values
global Sigma_u_1_boot
global Sigma_u_2_boot
global T_2_boot
global T_1_boot

global M
global ParamNumberA
global ParamNumberC
global A_boot
global C_boot

%Declaring the structure of the matrix A




for c_par = 1 : size(ParamNumberA,1)
    A_boot(ParamNumberA(c_par,1))=teta(c_par);   
end

     


    for c_par = 1 : size(ParamNumberC,1)
    C_boot(ParamNumberC(c_par,1))=teta(c_par+ size(ParamNumberA,1));   %this for loop puts correctly the imput values in the matrix A according to the positions declared in ParamNumberB
    end

    SIGMA_1=inv(A_boot)*inv(A_boot');
    SIGMA_2=inv(A_boot+C_boot)*inv(A_boot+C_boot)'; 
    SIGMA_1_inv=inv(SIGMA_1);
    SIGMA_2_inv=inv(SIGMA_2);

    
logLik = -(-0.5*(log(2*pi)*(T_1_boot*M+T_2_boot*M)) -0.5*T_1_boot*(log(det(SIGMA_1))+trace(Sigma_u_1_boot*SIGMA_1_inv  ) ) -0.5*T_2_boot*(log(det(SIGMA_2))+trace(Sigma_u_2_boot*SIGMA_2_inv  ) ) ) ;                                                                              

 %the loglikelihhod formula usually used to estimate the variance-covariance matrix after having estimated the reduced form parameters 
                                                                                                     
%+log(2*pi)*(T_1*M+T_2*M) to chek if it's correct to add this term. However
%it's a constant and  doesn't affect the maximization

% putting the '+' sign in front of each term because the fminun command
% minimizes. The original log likelihood has the '-' sign in front of each
% term
    
 end