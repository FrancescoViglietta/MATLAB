function [MMB_residuals]= MMB_res(Residuals_for_blocks, T, l, M) %% function to be applied on original residuals.Yields as output the rsiduals ready to build new datasets.
%the classical residuals coming from a var model should be putted with the
%Residuals'
%l should be a factor of T
 N=T/l;
 i=0:((T)-l);
 
 Residuals_for_blocks=Residuals_for_blocks';
 
 for  i=0:((T)-l)
             B=Residuals_for_blocks(1:M,i+1:i+l);
             B_blocks{:,i+1}=B;
 end
 i=0:((T)-l);
 r_i=randi([1,T-l+1],size((i+1:N),2),1);
 
 for par=1:size((r_i),1)
             B_rand(:,par)=B_blocks(:,r_i(par));
    
 end

 RAND_blocks=cell2mat(B_rand);

  for r=0:(T-l)
             for s=1:l
              Exp_u=(1/T-l+1)*sum(RAND_blocks(:,s+r), 2);   %centering all the residuals
             end   
  end


  for d=1:size(RAND_blocks,2)
       RAND_cent_blocks(:,d)=RAND_blocks(:,d)-Exp_u;

  end

         
    MMB_residuals=RAND_cent_blocks';

%author: Francesco Viglietta based on Bruggeman et al 2016
