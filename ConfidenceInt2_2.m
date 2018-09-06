%This function is for the confidence intervals. it should return the
%confidence intervals for the data and the parameters.
%Unlike 'ConfidenceInt2' this function estimates CIs for a stepwise FOI
%function that is not defined by contact probabilities.

%5th August 2014

function results=ConfidenceInt2_2(Data,Num,cutoff,n,x,p)

%Data=the data set to be resampled. contains ages and titre respectively
%Num=the number of samples to reproduce
%cutoff=the cutoff for seropositivity
%n=number of M classes
%x=initial guess for optimization
%p=proportion that loses seropositivity post primary infection

titre=reshape(Data(:,2),48,20);
%__________________________________________________________________________
%________________________Bootstrapping_____________________________________
iRandom=zeros(48,20,Num);
for i = 1:Num
    iRandom(:,:,i) = ceil(length(titre(:,1))*rand(length(titre(:,1)),length(titre(1,:))));
    %length(data(:,1)) sets the upper bound for the random numbers, and
    %ceil rounds off towards infinity, ie ceil(a) is >=a.eg ceil(0.01)=1 so
    %the random numbers will be between 1 and 48. every third dimension
    %element of iRandom is a different sample. 
end
% Initialize where you store the sampled data,preallocation before looping 
% saves on running time
sampledT = 0*iRandom;
ungrouped = zeros(960,Num);
Pos=zeros(20,Num);

for k=1:Num% Num samples
    for j=1:length(iRandom(1,:,1))%20 age groups
        y=zeros(48,1);
        for i=1:length(iRandom(:,1,1))%48 individuals
            sampledT(i,j,k)=titre(iRandom(i,j,k),j);
            if sampledT(i,j,k)>=cutoff
                    y(i)=1;
            end
        end
        Pos(j,k)=nnz(y); 
    end
    ungrouped(:,k)=reshape(sampledT(:,:,k),960,1);
end
%__________________________________________________________________________
%________________________CI's for the data_________________________________
CId=zeros(20,2); 
Prop=(1/48)*Pos;%Proportions seropositive at every age group for every sample
for m=1:20
    z=Prop(m,:);%the proportions positive for the mth age group from all the samples
    Z=sort(z);
    CId(m,:)=quantile(Z,[0.05,0.95]);%confidence interval for the mth age group
    clear z Z
end
R=[48;44;36;35;22;13;12;19;17;21;19;20;25;27;40;45;48;48;48;48];
RealProp=R/48;
%figure;errorbar(1:1:20,RealProp,(RealProp-CId(:,1)),(CId(:,2)-RealProp),'r*')
%__________________________________________________________________________
%__________________________Restimating_____________________________________
opt1=zeros(Num,3);fval1=zeros(Num,1);exitflag1=zeros(Num,1);%preallocation of variables 
F1=zeros(20,Num);
data=zeros(624,3,Num);
age=[0.0417;0.125;0.2083;0.2917;0.375;0.4583;0.542;0.625;0.708;0.792;0.875;...
    0.958;1.125;1.375;1.75;2.5;4;6;8;10.542]; 

lb=0*x;
options=optimset('MaxFunEvals',1000*length(x));

for k=1:Num
    data(:,:,k)=Regroup(Data(:,1),ungrouped(:,k),cutoff);
    Age=(data(:,1,k))/12;% age in years
    Nr=data(:,2,k);% Numbers at every age
    Rr=data(:,3,k);% Numbers seropositive at every age
    delta=4;%Rate of loss of primary seropositive status
 
    [opt1(k,:),fval1(k,:),exitflag1(k,:)]=fmincon(@msirsirmax15,x,[],[],[],[],lb,[],[],options,Age,Nr,Rr,delta,p,n);
    %solving with the optimised values
    y0=zeros(n+4,1);y0(1)=0.999;y0(n+1)=0.001;
    [a,y]=ode45(@msirsir15,age,y0,options,opt1(k,:),delta,p,n);
    F1(:,k)=sum(y(:,1:n),2)+y(:,n+2)+y(:,n+3);
end

figure;
plot(age,F1,'b');xlabel('Age in years');ylabel('Proportions Seropositive')
hold on;errorbar(age,RealProp,(RealProp-CId(:,1)),(CId(:,2)-RealProp),'r*')
title('the data and the MSIR confidence region')

results=opt1;
return