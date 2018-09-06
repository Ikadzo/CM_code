%This is the likelihood function for exploring age dependen FOI,i.e.
%stepwise function,for the MSIR model and maybe MSIRSIR1,MSIRSIR2 would be 
%redundant. The steps are for ages 0-6m,6-12m,12-48m, 48-144m


function L=msirsirmax15(x,age,N,R,delta,p,n)

y0=zeros(n+4,1);y0(1)=0.999;y0(n+1)=0.001;
options = odeset('AbsTol',1e-40);


%Solving the ODEs
[a,y]=ode23(@msirsir15,age,y0,[],x,delta,p,n);


%y does not strictly add up to 1 across the ages and this brings problems
%when getting the log, so below, i normalize it to avoid the issues

Sum=sum(y,2);
for i=1:length(age)
    Y(i,:)=y(i,:)./Sum(i);
end

%The proportions seropositive
F=sum(Y(:,1:n),2)+Y(:,n+2)+Y(:,n+3);

%Getting the log likelihood
C=(factorial(N))./((factorial(R)).*(factorial(N-R)));%the constant I usually ignore
logl=log(C)+(R.*log(F)+(N-R).*log(1-F));
L=sum(logl);
L=-L;


end
