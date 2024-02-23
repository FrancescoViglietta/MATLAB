%% 
clear 
clc
close all
clear global     % cleans all gloabl variables
%% Global Variables 
global p             % # of VAR lags
global T_overall  
global T_1 
global T_2 % # of observations in dataset   
global Sigma_u_1  
global Sigma_u_2 % error covariance matrix
global ParamNumberA  %this vector states the positions of the free parameters to be estimated in the vectorized matrix A
global A
global C
global ParamNumberC %this vector states the positions of the free parameters to be estimated in the vectorized matrix B
global M 
global T_break %it corresponds with the row number associated with the break in the original dataset (see excel file)


 %% Data set


load DATI_FISCAL_CK.txt;
%cd 'C:\Users\Utente\Desktop\Project_structural Macro'





% Variables per column in the dataset: 
%  1-GDP (real gdp per capita)
%  4-TAX (total tax revenues per capita)
%  5-G (Fiscal spending in real per capita term)
% Real interest rate
%real interest rate
RR=DATI_FISCAL_CK(:,3)-DATI_FISCAL_CK(:,2);




%selecting the variables for the VAR
Variablesselection=[4,5,1]; %the first variable is TAX, the second is G and the last is GDP

Var_dataset=[DATI_FISCAL_CK(:,Variablesselection),RR]; %this line extract the variables of interest from the original dataset

%% VAR Estimation

M=4;
p=4;

VAR_Const = [NaN; NaN; NaN; NaN]; %setting the vector of the parameters associated with th constant

 Pi1 =[NaN NaN NaN NaN; 
       NaN NaN NaN NaN;
       NaN NaN NaN NaN;
       NaN NaN NaN NaN] ;      % Each matrix is a 4x4 and will collect the parameters for each lag
                          
 Pi2 =[NaN NaN NaN NaN; 
       NaN NaN NaN NaN;
       NaN NaN NaN NaN;
       NaN NaN NaN NaN] ; 


Pi3 =[NaN NaN NaN NaN; 
       NaN NaN NaN NaN;
       NaN NaN NaN NaN;
       NaN NaN NaN NaN] ;

Pi4 =[NaN NaN NaN NaN; 
       NaN NaN NaN NaN;
       NaN NaN NaN NaN;
       NaN NaN NaN NaN] ;




 


Linear_trend=[NaN; NaN; NaN; NaN];

VAR_Pi = {Pi1 Pi2 Pi3 Pi4};  % collecting all the matrices of parameters in one


 %Setting VAR (4)
 
VAR = varm('Constant',VAR_Const,'AR',VAR_Pi,'Trend', Linear_trend);%specification of the VAR(4) with a constant and a linear trend

[EstVAR,EstSE,logLikVAR,Residuals] = estimate(VAR,Var_dataset); % Estimation of the above specified var



Const = EstVAR.Constant;
Trend= EstVAR.Trend;
mP1 = EstVAR.AR{1,1};
mP2 = EstVAR.AR{1,2};
mP3 = EstVAR.AR{1,3};
mP4 = EstVAR.AR{1,4};



%% Analyzing the residuals
T_overall=(size(Var_dataset,1)-p);

T_break=136;

%T_break=118 (1979 break)
 
%First we obtain the variance covariance matrices from 1-Overall sample; 2-From the start to 1979-Q2 ("pre Volker era"); 3_From 1972_Q3 to the end of the sample (Great-Moderation era) 

 %Overall residuals covariance matrix
 Sigma_u_overall=(Residuals'*Residuals)/T_overall;

 % From 1950-Q1 to 1979-Q2
 Residuals_preV=(Residuals(1:(T_break-p),1:4));
 Sigma_u_preV=(Residuals_preV'*Residuals_preV)/(T_break-p);


 % From 1979-Q2 to end
 Residuals_GM=(Residuals((T_break-p+1):end, 1:4));
 Sigma_u_GM=(Residuals_GM'*Residuals_GM)/(T_overall-size(Residuals_preV,1));

Sigma_u_1=Sigma_u_preV;
Sigma_u_2= Sigma_u_GM;

%% Second we want to test if the variance covariance parameters in the two
%regimes are statistically different

  %in order to do this, we first compute the log likelihood of the model by
  %keeping the reduced form parameters equal to those of the overall sample
  %but changing the varaince coviariance matrix estimator, one time with the
  %pre volker one and, the other time with the great moderation one

 T_1=size(Residuals_preV, 1);
  
 logLik_preV = (-0.5*T_1*M*(log(2*pi)) + 0.5*T_1*log(det(inv(Sigma_u_preV))) - 0.5*T_1*trace(inv(Sigma_u_preV)*Sigma_u_overall));
 
 
 T_2=size(Residuals_GM,1);
 logLik_GM = (-0.5*T_2*M*(log(2*pi)) + 0.5*T_2*log(det(inv(Sigma_u_GM))) - 0.5*T_2*trace(inv(Sigma_u_GM)*Sigma_u_overall));
 
 
 loglik_sum=(logLik_preV+logLik_GM);



 LR_test=(-2*(logLikVAR-loglik_sum));
 


 



%% Identification scheme
 rng(4) %use this seed
 
 A=[NaN  0 NaN  0; 
    0  NaN   0 0;
    NaN NaN  NaN 0;
    0 NaN NaN  NaN];
 
  C =[NaN  0 NaN 0; 
   0    NaN   0 0;
   NaN NaN  NaN 0;
   0 NaN NaN  NaN];
ParamNumberA=[ 1, 3,  6, 7, 8, 9, 11, 12, 16]';
ParamNumberC=[ 1, 3,  6, 7, 8, 9, 11, 12, 16]';

 
 StructuralParam_A = size(ParamNumberA,1);
 StructuralParam_C = size(ParamNumberC,1);


 Initial_AC=(randn(StructuralParam_A+StructuralParam_C,1)/10);



SE_A=zeros(size(A,1),size(A,1));
HSelectionA = zeros(M*M,StructuralParam_A);

SE_C=zeros(size(C,1),size(C,1));
HSelectionC = zeros(M*M,StructuralParam_C);


Sigma_u_1=Sigma_u_preV;
Sigma_u_2=Sigma_u_GM;

%%Likelihood Maximization:the function is the one specified in Lanne & Lutekpohl but with the Fanelli methodology
options = optimset('MaxFunEvals',20000000,'TolFun',1e-500,'MaxIter',200000000,'TolX',1e-50);   
[StructuralParam_Estimation_MATRIX,Likelihood_SVAR,exitflag, output, grad,Hessian_MATRIX] = fminunc('Likelihood_Fanelli_Heterosk', Initial_AC', options);

SE_Hessian_MATRIX = diag(inv(Hessian_MATRIX)).^0.5;  % computes Hessian-based standard errors



%%Filling the matrices with the estimated parameters
for c_par= 1:  size(ParamNumberA,1)
     A(ParamNumberA(c_par,1))= StructuralParam_Estimation_MATRIX(c_par); %filling the A  matrix with the estimated parameters
     SE_A(ParamNumberA(c_par,1))= SE_Hessian_MATRIX(c_par); %filling the standard errors matrix
     HSelectionA(ParamNumberA(c_par,1),c_par) = 1;
end  


for c_par= 1:  size(ParamNumberC,1)
     C(ParamNumberC(c_par,1))= StructuralParam_Estimation_MATRIX(c_par+size(ParamNumberA,1)); %filling the C  matrix with the estimated parameters
     SE_C(ParamNumberC(c_par,1))= SE_Hessian_MATRIX(c_par+size(ParamNumberA,1)); %filling the standard errors matrix
     HSelectionC(ParamNumberC(c_par,1),c_par) = 1;
end  

A
AC=(A+C)
A_inv=inv(A)
AC_inv=inv(AC)

%% Sign normalization  for A nd A+C
   
    if A(1,1)<0
    A(:,1)=-A(:,1);
    end
    if A(2,2)<0
    A(:,2)=-A(:,2);
    end
    if A(3,3)<0                 %this sign normalization allows to reason always with an increase in all the variables on the matrices used to compute IRF
   A(:,3)=-A(:,3);
    end
    if A(4,4)<0
    A(:,4)=-A(:,4);
    end



    if AC(1,1)<0
    AC(:,1)=-AC(:,1);
    end
    if AC(2,2)<0
    AC(:,2)=-AC(:,2);
    end
    if AC(3,3)<0
    AC(:,3)=-AC(:,3);
    end
    if AC(4,4)<0
    AC(:,4)=-AC(:,4);
    end

A
AC
A_inv=inv(A)
AC_inv=inv(AC)


%% Rank condition without cross-equation restrictions (only for heterosk approach)   
DuplicationMatrix = DuplicationMatrixFunction(M);
mDD=(DuplicationMatrix'*DuplicationMatrix)^(-1)*DuplicationMatrix';  
mNN=DuplicationMatrix*mDD;


Left_matrix=kron(eye(2), mDD);

I_block=kron(A_inv,eye(M));
II_block=zeros(M*M);
III_block=kron(AC_inv,eye(M));
Middle_matrix= [I_block, II_block;
                  III_block, III_block];

block_ofzeros_C=zeros(M*M,size(ParamNumberC,1));
block_ofzeros_A=zeros(M*M,size(ParamNumberA,1));
Right_matrix= [HSelectionA, block_ofzeros_C; block_ofzeros_A, HSelectionC];
                     
Jacobian= (Left_matrix*Middle_matrix*Right_matrix) ;



rank(Jacobian)

(size(ParamNumberA,1)+size(ParamNumberC,1))


 
 
%% Generating the selection vectors and average ratios for the fiscal and the tax multiplier

e_i=[1,0,0,0]'; %this vector select the tax shock
e_j=[0,0,1,0]';%this vector select GDP
b_i=[0,1,0,0]';%this vector selects the G shock

unlogged_gdp= exp(DATI_FISCAL_CK(:,1)); 
unlogged_tax=exp(DATI_FISCAL_CK(:,4));
unlogged_g=exp(DATI_FISCAL_CK(:,5));

Dynamic_tax_multiplier_ratio_1=mean(unlogged_gdp(1:T_1+p,1)./unlogged_tax(1:T_1+p,1));
Dynamic_tax_multiplier_ratio_2=mean(unlogged_gdp(T_2+p:end,1)./unlogged_tax(T_2+p:end,1));%average of the unlogged gdp and unlogged tax ratio

FiscalSp_multiplier_ratio_1=mean(unlogged_gdp(1:T_1+p,1)./unlogged_g(1:T_1+p,1));
FiscalSp_multiplier_ratio_2=mean(unlogged_gdp(T_2+p:end,1)./unlogged_g(T_2+p:end,1));
%% Moving block bootstrap for heteroskedastic FanelliBachiocchi approach



global A_boot
global C_boot

global T_1_boot
global T_2_boot
global Sigma_u_1_boot
global Sigma_u_2_boot


BootstrapIterations=1000;
HorizonIRF = 20; %horizon 
quant = [5,95];


for boot=1:BootstrapIterations

    Res_star_1=MMB_res(Residuals_preV, T_1, 4, M);
    Res_star_2=MMB_res(Residuals_GM, T_2, 4, M);
    Residuals_Boot=[Res_star_1;Res_star_2];


 DataSet_Bootstrap=zeros(T_overall+p,M);

  %DataSet_Bootstrap(1:p,:)=Var_dataset(1:p,:); % set the first p elements
  %equal to the original sample values, for the MMB the first p observation are setted to zero...... 

       
       for t = p+1 : T_overall+p
             DataSet_Bootstrap(t,:)=Const + Trend + mP1 * DataSet_Bootstrap(t-1,:)' +...
                                       mP2 * DataSet_Bootstrap(t-2,:)' + ...
                                       mP3 * DataSet_Bootstrap(t-3,:)' + ...
                                       mP4 * DataSet_Bootstrap(t-4,:)' + ...
                                       Residuals_Boot(t-p,:)';
       end

      DataSet_Bootstrap=DataSet_Bootstrap(end-T_overall+1:end,:);   %......and then deleted


 T_1_boot=(T_break-p-p);
 T_2_boot=((T_overall-p)-(T_break-p-p));


  

[EstVAR_Boot,EstSE_Boot,logLikVAR_Boot,Residuals_Boot] = estimate(VAR,DataSet_Bootstrap); 
    mP1_Boot = EstVAR_Boot.AR{1,1};
    mP2_Boot = EstVAR_Boot.AR{1,2};                                          
    mP3_Boot = EstVAR_Boot.AR{1,3};
    mP4_Boot = EstVAR_Boot.AR{1,4};
    
 Sigma_u_overall_boot=(Residuals_Boot'*Residuals_Boot)/T_overall-p;

 %First period bootstrapped variance-covariance matrix
 
 Residuals_preV_boot=(Residuals_Boot(1:(T_break-p-p),1:4));
 Sigma_u_preV_boot=(Residuals_preV_boot'*Residuals_preV_boot)/(T_break-p-p);
 

 %Second period bootstrapped variance-covariance matrix
 
 Residuals_GM_boot=(Residuals_Boot((T_break-p-p+1):220, 1:4));
 Sigma_u_GM_boot=(Residuals_GM_boot'*Residuals_GM_boot)/((T_overall-p)-(T_break-p-p));

 
 Sigma_u_1_boot=Sigma_u_preV_boot;
 Sigma_u_2_boot= Sigma_u_GM_boot;


 
 StructuralParam_A_boot = size(ParamNumberA,1);
 StructuralParam_C_boot = size(ParamNumberC,1);
 Initial_AC_boot=StructuralParam_Estimation_MATRIX';


  A_boot=[NaN  0 NaN  0; 
   0  NaN   0 0;
   NaN NaN  NaN   0;
   0 NaN NaN  NaN];
 
  C_boot =[NaN  0 NaN  0; 
   0    NaN   0 0;
   NaN NaN  NaN   0;
   0 NaN NaN  NaN];



options = optimset('MaxFunEvals',200000,'TolFun',1e-12,'MaxIter',200000,'TolX',1e-12);   
[StructuralParam_Estimation_MATRIX_boot,Likelihood_SVAR,exitflag, output, grad,Hessian_MATRIX] = fminunc('Likelihood_Fanelli_Heterosk_boot', Initial_AC_boot', options);

 



%%Filling the matrices with the estimated parameters
for c_par= 1:  size(ParamNumberA,1)
     A_boot(ParamNumberA(c_par,1))= StructuralParam_Estimation_MATRIX_boot(c_par); %filling the A  matrix with the estimated parameters
     
end  


for c_par= 1:  size(ParamNumberC,1)
     C_boot(ParamNumberC(c_par,1))= StructuralParam_Estimation_MATRIX_boot(c_par+size(ParamNumberA,1)); %filling the C  matrix with the estimated parameters
     
end  

AC_boot=(A_boot+C_boot);

   if A_boot(1,1)<0
    A_boot(:,1)=-A_boot(:,1);
    end
    if A_boot(2,2)<0
    A_boot(:,2)=-A_boot(:,2);
    end
    if A_boot(3,3)<0
   A_boot(:,3)=-A_boot(:,3);
    end
    if A_boot(4,4)<0
    A_boot(:,4)=-A_boot(:,4);
    end

   

    if AC_boot(1,1)<0
    AC_boot(:,1)=-AC_boot(:,1);
    end
    if AC_boot(2,2)<0
    AC_boot(:,2)=-AC_boot(:,2);
    end
    if AC_boot(3,3)<0
    AC_boot(:,3)=-AC_boot(:,3);
    end
    if AC_boot(4,4)<0
    AC_boot(:,4)=-AC_boot(:,4);
    end

     A_inv_boot=inv(A_boot);
     AC_inv_boot=inv(AC_boot);

J=[eye(M) zeros(M,M*(p-1))]; %the selection matrix for the companion form
    CompanionMatrix_Boot = [mP1_Boot mP2_Boot mP3_Boot mP4_Boot;
                            eye(M*(p-1)) zeros(M*(p-1),M)];  %getting the  bootstrapped companion matrix for the first period

  


    for h = 0 : HorizonIRF
    TETA_Boot_1(:,:,h+1,boot)=J*CompanionMatrix_Boot^h*J'*A_inv_boot;
    TETA_Boot_2(:,:,h+1,boot)=J*CompanionMatrix_Boot^h*J'*AC_inv_boot; %obtaining the impulse response functions(.,.) for each bootstrapped sample
    end 


    for h = 0 : HorizonIRF %all the following computations refers to the bootstrap case in order to get confidence intervals. The point estimates for the two mulpipliers will be in the  "IRF FOR THE MULTIPLIERS" section
         
    
    Impulse_tax_to_GDPh_boot_1= e_j'*J*CompanionMatrix_Boot^h*J'*A_inv_boot*e_i; %getting the impulse response of taxes on GDP for each of the 20 steps ahead
    Impulse_tax_to_GDPh_boot_2= e_j'*J*CompanionMatrix_Boot^h*J'*AC_inv_boot*e_i;
    
     
    Impulse_tax_to_TAX_boot_1=e_i'*A_inv_boot*e_i; %getting the on impact response of a shock in taxes on taxes
    Impulse_tax_to_TAX_boot_2=e_i'*AC_inv_boot*e_i;
    
    Impulse_g_to_GDPh_boot_1=e_j'*J*CompanionMatrix_Boot^h*J'*A_inv_boot*b_i; %getting the impulse response of fiscal spending on GDP for each of the 20 steps ahead 
    Impulse_g_to_GDPh_boot_2=e_j'*J*CompanionMatrix_Boot^h*J'*AC_inv_boot*b_i;
    
    Impulse_g_to_g_boot_1=b_i'*A_inv_boot*b_i; %getting the on impact response of a shock in fiscal spending  on fiscal spending
    Impulse_g_to_g_boot_2=b_i'*AC_inv_boot*b_i;
    
    
    Dynamic_tax_multiplier_boot_1(1,1,h+1,boot)=(Impulse_tax_to_GDPh_boot_1./Impulse_tax_to_TAX_boot_1)*Dynamic_tax_multiplier_ratio_1; %getting the bootrapped  dynamic tax multiplier. The round parenthesis for the dimensions have been placed because the selection vectors declared above only put the non interesting elements equal to zero and this affects the calculations
    Dynamic_tax_multiplier_boot_2(1,1,h+1,boot)=(Impulse_tax_to_GDPh_boot_2./Impulse_tax_to_TAX_boot_2)*Dynamic_tax_multiplier_ratio_2;
    
    
    
    Fiscal_spending_multiplier_boot_1(1,1,h+1,boot)=(Impulse_g_to_GDPh_boot_1./Impulse_g_to_g_boot_1)*FiscalSp_multiplier_ratio_1; %getting the bootrapped  fical spending multiplier.
    Fiscal_spending_multiplier_boot_2(1,1,h+1,boot)=(Impulse_g_to_GDPh_boot_2./Impulse_g_to_g_boot_2)*FiscalSp_multiplier_ratio_2;


     end 


   boot
 end

IRF_Inf_Boot_1 = prctile(TETA_Boot_1,quant(1),4); %this command selects the inferior bound for the selected confidence level on the bootstrapped IRFs
IRF_Sup_Boot_1 = prctile(TETA_Boot_1,quant(2),4);

IRF_Inf_Boot_2 = prctile(TETA_Boot_2,quant(1),4); %this command selects the inferior bound for the selected confidence level on the bootstrapped IRFs
IRF_Sup_Boot_2 = prctile(TETA_Boot_2,quant(2),4);


IRF_Inf_Boot_DynamTax_multi_1 = prctile(Dynamic_tax_multiplier_boot_1,quant(1),4); %the following lines select the inferior and the superior bounds for the selected confidence level on the bootstrapped multipliers obtianed above
IRF_Sup_Boot_DynamTax_multi_1 = prctile(Dynamic_tax_multiplier_boot_1,quant(2),4);

IRF_Inf_Boot_DynamTax_multi_2 = prctile(Dynamic_tax_multiplier_boot_2,quant(1),4); %the following lines select the inferior and the superior bounds for the selected confidence level on the bootstrapped multipliers obtianed above
IRF_Sup_Boot_DynamTax_multi_2 = prctile(Dynamic_tax_multiplier_boot_2,quant(2),4);

IRF_Inf_Boot_Fiscalsp_multi_1 = prctile(Fiscal_spending_multiplier_boot_1,quant(1),4);
IRF_Sup_Boot_Fiscalsp_multi_1 = prctile(Fiscal_spending_multiplier_boot_1,quant(2),4);

IRF_Inf_Boot_Fiscalsp_multi_2 = prctile(Fiscal_spending_multiplier_boot_2,quant(1),4);
IRF_Sup_Boot_Fiscalsp_multi_2 = prctile(Fiscal_spending_multiplier_boot_2,quant(2),4);


%% Impulse response for the preVolker era
HorizonIRF=20;
LineWidth_IRF = 1.5;
FontSizeIRFGraph = 14;

C_IRF_1 = A_inv;     % instantaneous impact at h=0

J=[eye(M) zeros(M,M*(p-1))];                          % selection matrix J used in IRF computation 

CompanionMatrix = [mP1 mP2 mP3 mP4;           % VAR companion matrix
                   eye(M*(p-1)) zeros(M*(p-1),M)];
                       
    for h = 0 : HorizonIRF
    TETA_1(:,:,h+1)=J*CompanionMatrix^h*J'*C_IRF_1; %computing the point estimates on the original sample (notice that these are different fro the bootstrapped ones)
    end
  
IRF_Inf_Boot_1 =(2*TETA_1)- prctile(TETA_Boot_1,quant(1),4); %Hall's percentile intervals
IRF_Sup_Boot_1 =(2*TETA_1)- prctile(TETA_Boot_1,quant(2),4);



Titles=cell(1,4);
Titles{1,1}='$$\varepsilon_{TAX}$$ $$Shock$$';
Titles{1,2}='$$\varepsilon_{G}$$ $$Shock$$';
Titles{1,3}='$$\varepsilon_{GDP}$$ $$Shock$$';          %setting properties of the graph
Titles{1,4}='$$\varepsilon_{RR}$$ $$Shock$$';



YLabel=cell(3,1);
YLabel{1,1}='$$TAX$$';
YLabel{2,1}='$$G$$';
YLabel{3,1}='$$GDP$$';
YLabel{4,1}='$$RR$$';
index = 1;
Shock_1=[1 1 1 1];

M_IRF = 4;

figure(1)
for jr = 1 : M_IRF
    for jc = 1 : M_IRF
    TETA_Iter_Sample = squeeze(TETA_1(jr,jc,:));
    TETA_Iter_Boot_Inf = squeeze(IRF_Inf_Boot_1(jr,jc,:));
    TETA_Iter_Boot_Sup = squeeze(IRF_Sup_Boot_1(jr,jc,:));
    subplot(M_IRF,M_IRF,index)  
    x = 1:1:HorizonIRF+1;
    y= TETA_Iter_Sample'*Shock_1(jr);
    plot(y,'Color',[0 0.4470 0.7410], 'LineWidth',LineWidth_IRF);     
    hold all
    plot(TETA_Iter_Boot_Inf,'--r', 'LineWidth',LineWidth_IRF);
    plot(TETA_Iter_Boot_Sup,'--r', 'LineWidth',LineWidth_IRF);           
    plot(zeros(HorizonIRF+1,1),'k','LineWidth',1);
    ylabel(YLabel{jr,1},'interpreter','latex');
    title(Titles{1,jc},'interpreter','latex');
    set(gca,'FontSize',FontSizeIRFGraph);
    axis tight
    index=index+1;
    
    end


end



%% Impulse response for the great moderation

LineWidth_IRF = 1.5;
FontSizeIRFGraph = 14;

C_IRF_2 = AC_inv;     % instantaneous impact at h=0

J=[eye(M) zeros(M,M*(p-1))];                          % selection matrix J used in IRF computation 

CompanionMatrix = [mP1 mP2 mP3 mP4;           % VAR companion matrix
                   eye(M*(p-1)) zeros(M*(p-1),M)];
                       
    for h = 0 : HorizonIRF
    TETA_2(:,:,h+1)=J*CompanionMatrix^h*J'*C_IRF_2; %computing the point estimates on the original sample (notice that these are different fro the bootstrapped ones)
    end
IRF_Inf_Boot_2 =(2*TETA_2)-prctile(TETA_Boot_2,quant(1),4); %Hall's percentile interval
IRF_Sup_Boot_2 =(2*TETA_2)-prctile(TETA_Boot_2,quant(2),4);    




Titles=cell(1,4);
Titles{1,1}='$$\varepsilon_{TAX}$$ $$Shock$$';
Titles{1,2}='$$\varepsilon_{G}$$ $$Shock$$';
Titles{1,3}='$$\varepsilon_{GDP}$$ $$Shock$$';          %setting properties of the graph
Titles{1,4}='$$\varepsilon_{RR}$$ $$Shock$$';



YLabel=cell(3,1);
YLabel{1,1}='$$TAX$$';
YLabel{2,1}='$$G$$';
YLabel{3,1}='$$GDP$$';
YLabel{4,1}='$$RR$$';
index = 1;
Shock_1=[1 1 1 1];

M_IRF = 4;

figure(2)
for jr = 1 : M_IRF
    for jc = 1 : M_IRF
    TETA_Iter_Sample = squeeze(TETA_2(jr,jc,:));
    TETA_Iter_Boot_Inf = squeeze(IRF_Inf_Boot_2(jr,jc,:));
    TETA_Iter_Boot_Sup = squeeze(IRF_Sup_Boot_2(jr,jc,:));
    subplot(M_IRF,M_IRF,index)  
    x = 1:1:HorizonIRF+1;
    y= TETA_Iter_Sample'*Shock_1(jr);
    plot(y,'Color',[0 0.4470 0.7410], 'LineWidth',LineWidth_IRF);      
    hold all
    plot(TETA_Iter_Boot_Inf,'--r', 'LineWidth',LineWidth_IRF);
    plot(TETA_Iter_Boot_Sup,'--r', 'LineWidth',LineWidth_IRF);           
    plot(zeros(HorizonIRF+1,1),'k','LineWidth',1);
    ylabel(YLabel{jr,1},'interpreter','latex');
    title(Titles{1,jc},'interpreter','latex');
    set(gca,'FontSize',FontSizeIRFGraph);
    axis tight
    index=index+1;
    
    end
end






%% IRFs FOR THE MULTIPLIERS
HorizonIRF=20; %number of time steps ahead the shock

    % instantaneous impact at h=0

J=[eye(M) zeros(M,M*(p-1))];
CompanionMatrix = [mP1 mP2 mP3 mP4;           
                   eye(M*(p-1)) zeros(M*(p-1),M)]; %selection and companion matrix for the VAR VMA representation


                                             
    for h = 0 : HorizonIRF
    
    Impulse_tax_to_GDPh_1=e_j'*J*CompanionMatrix^h*J'*C_IRF_1*e_i;
    Impulse_tax_to_GDPh_2=e_j'*J*CompanionMatrix^h*J'*C_IRF_2*e_i;%properly selecting the shocks as did for the bootstrap case (here the multipliers are computed on the whole sample to obtain the point estimates)
    
    Impulse_tax_to_TAX_1=e_i'*C_IRF_1*e_i;
    Impulse_tax_to_TAX_2=e_i'*C_IRF_2*e_i;
   
    
    Impulse_g_to_GDPh_1=e_j'*J*CompanionMatrix^h*J'*C_IRF_1*b_i;
    Impulse_g_to_GDPh_2=e_j'*J*CompanionMatrix^h*J'*C_IRF_2*b_i;
    
    Impulse_g_to_g_1=b_i'*C_IRF_1*b_i;
    Impulse_g_to_g_2=b_i'*C_IRF_2*b_i;


        
        
    Dynamic_tax_multiplier_1(h+1,1)=(Impulse_tax_to_GDPh_1./Impulse_tax_to_TAX_1)*Dynamic_tax_multiplier_ratio_1;
    Dynamic_tax_multiplier_2(h+1,1)=(Impulse_tax_to_GDPh_2./Impulse_tax_to_TAX_2)*Dynamic_tax_multiplier_ratio_2;

    Fiscal_spending_multiplier_1(h+1,1)=(Impulse_g_to_GDPh_1./Impulse_g_to_g_1)*FiscalSp_multiplier_ratio_1;  
    Fiscal_spending_multiplier_2(h+1,1)=(Impulse_g_to_GDPh_2./Impulse_g_to_g_2)*FiscalSp_multiplier_ratio_2; %obtaining the multipliers
    end

    
    
    











INF_IRF_BOOT_Dynamic_TAX_multiplier_1=(2*Dynamic_tax_multiplier_1)-squeeze(IRF_Inf_Boot_DynamTax_multi_1);  %for a matter of usefulness, the  command "squeeze" reduces an object of three dimensions in a vector (easier to plot)
SUP_IRF_BOOT_Dynamic_TAX_multiplier_1=(2*Dynamic_tax_multiplier_1)-squeeze(IRF_Sup_Boot_DynamTax_multi_1);
INF_IRF_BOOT_Dynamic_TAX_multiplier_2=(2*Dynamic_tax_multiplier_2)-squeeze(IRF_Inf_Boot_DynamTax_multi_2);  
SUP_IRF_BOOT_Dynamic_TAX_multiplier_2=(2*Dynamic_tax_multiplier_2)-squeeze(IRF_Sup_Boot_DynamTax_multi_2);

INF_IRF_BOOT_Fiscal_spending_multiplier_1=(2*Fiscal_spending_multiplier_1)-squeeze(IRF_Inf_Boot_Fiscalsp_multi_1);
SUP_IRF_BOOT_Fiscal_spending_multiplier_1=(2*Fiscal_spending_multiplier_1)-squeeze(IRF_Sup_Boot_Fiscalsp_multi_1);
INF_IRF_BOOT_Fiscal_spending_multiplier_2=(2*Fiscal_spending_multiplier_2)-squeeze(IRF_Inf_Boot_Fiscalsp_multi_2);
SUP_IRF_BOOT_Fiscal_spending_multiplier_2=(2*Fiscal_spending_multiplier_2)-squeeze(IRF_Sup_Boot_Fiscalsp_multi_2);

%% Multipliers 
% Define time steps
time_steps = 0:20;


figure(3);


% Subplot 1
subplot(2,2,1);
hold on;
fill([time_steps, fliplr(time_steps)], [INF_IRF_BOOT_Dynamic_TAX_multiplier_1', fliplr(SUP_IRF_BOOT_Dynamic_TAX_multiplier_1')], 'blue', 'FaceAlpha', 0.3, 'EdgeColor','none');
plot(time_steps, Dynamic_tax_multiplier_1, '-', 'Color', [0.6350, 0.0780, 0.1840], 'LineWidth', 1.8);
line([min(time_steps), max(time_steps)], [0, 0], 'Color', 'k', 'LineStyle', '-', 'LineWidth', 1);
hold off;
xlabel('Time Horizon');
ylabel('');
title('Dynamic Tax multiplier Pre-Volker');
ax1 = gca; 

% Subplot 2
subplot(2,2,2);
hold on;
fill([time_steps, fliplr(time_steps)], [INF_IRF_BOOT_Dynamic_TAX_multiplier_2', fliplr(SUP_IRF_BOOT_Dynamic_TAX_multiplier_2')], 'blue', 'FaceAlpha', 0.3,'EdgeColor','none');
plot(time_steps, Dynamic_tax_multiplier_2, '-', 'Color', [0.6350, 0.0780, 0.1840], 'LineWidth', 1.8);
line([min(time_steps), max(time_steps)], [0, 0], 'Color', 'k', 'LineStyle', '-', 'LineWidth', 1);
hold off;
xlabel('Time Horizon');
ylabel('');
title('Dynamic Tax multiplier Great Moderation');
ax2 = gca; 

% Subplot 3
subplot(2,2,3);
hold on;
fill([time_steps, fliplr(time_steps)], [INF_IRF_BOOT_Fiscal_spending_multiplier_1', fliplr(SUP_IRF_BOOT_Fiscal_spending_multiplier_1')], 'blue', 'FaceAlpha', 0.3,'EdgeColor','none');
plot(time_steps, Fiscal_spending_multiplier_1, '-', 'Color', [0.6350, 0.0780, 0.1840], 'LineWidth', 1.8);
line([min(time_steps), max(time_steps)], [0, 0], 'Color', 'k', 'LineStyle', '-', 'LineWidth', 1);
hold off;
xlabel('Time Horizon');
ylabel('');
title('Fiscal Spending multiplier pre-Volker');
ax3 = gca; 

% Subplot 4
subplot(2,2,4);
hold on;
fill([time_steps, fliplr(time_steps)], [INF_IRF_BOOT_Fiscal_spending_multiplier_2', fliplr(SUP_IRF_BOOT_Fiscal_spending_multiplier_2')], 'blue', 'FaceAlpha', 0.3,'EdgeColor','none');
plot(time_steps, Fiscal_spending_multiplier_2, '-', 'Color', [0.6350, 0.0780, 0.1840], 'LineWidth', 1.8);
line([min(time_steps), max(time_steps)], [0, 0], 'Color', 'k', 'LineStyle', '-', 'LineWidth', 1);
hold off;
xlabel('Time Horizon');
ylabel('');
title('Fiscal Spending multiplier Great Moderation');
ax4 = gca; 

% Ensure y-axes have the same scale
y_limits1 = [min([ax1.YLim, ax2.YLim]), max([ax1.YLim, ax2.YLim])];
y_limits2 = [min([ax3.YLim, ax4.YLim]), max([ax3.YLim, ax4.YLim])];
ax1.YLim = y_limits1;
ax2.YLim = y_limits1;
ax3.YLim = y_limits2;
ax4.YLim = y_limits2;


