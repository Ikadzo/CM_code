%This is the set of ODES for exploring age dependent FOI,i.e.
%stepwise function,for the MSIR model and maybe MSIRSIR1,MSIRSIR2 would be 
%redundant. The steps are for ages 0-12m,12-144m

function dyda=msirsir15(t,y,x,delta,p,n)

sigma=x(1);       
if t<=1
    lambda=x(2);
else
    lambda=x(3);
end
lambda1=lambda;
        
%M1
dyda(1,1)=-n*sigma*y(1);
%M2...Mn
for k=2:n
    dyda(k,1)=n*sigma*y(k-1)-n*sigma*y(k);
end
%So1
dyda(n+1,1)=n*sigma*y(n)-lambda*y(n+1);
%Fo
dyda(n+2,1)=lambda*(1-p)*y(n+1)+lambda1*y(n+4);
%F1
dyda(n+3,1)=lambda*p*y(n+1)-delta*y(n+3);
%S1
dyda(n+4,1)=delta*y(n+3)-lambda1*y(n+4);

return