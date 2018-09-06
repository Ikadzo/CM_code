% This is a function similar to Catalyticmodel13 but now refitting the
% different models using subsets of the data:
%                   1. Minus those that were positive for RSV antigen
%                   2. Minus those that had LRTIs
% Written on 19th May 2016
%__________________________________________________________________________
cd('/Users/ikombe/Desktop/IKombe/CM_revisited/Codes')

% Parameters/ quantities
cutoff=1.5; % for seropositivity
delta=4; % rate of loss of Abs
lower=[(0:1:12),15,18,24,36,60,84,108];% lower bound gor age group
upper=[(1:1:12),15,18,24,36,60,84,108,145];
age=(lower+((upper-lower)/2))/12;% Mid point of age group
Num=1000;% Number of times to resample the data
%--------------------
% The ungrouped data
%--------------------
% 1. All the data(fitted to make sure code is working fine)
    All_data = xlsread('Short_version.xls');
    Age=All_data(:,1);
    LRTI=All_data(:,2);%346
    RSVpos=All_data(:,3);%64
    RSVlog=All_data(:,4);

    % preparing for fitting
    [data]=Regroup(Age,RSVlog,1.5);
    Age_1=(data(:,1))/12;% age in years
    Nr_1=data(:,2);% Numbers at every age
    Rr_1=data(:,3);% Numbers seropositive at every age
% 2. Minus RSV positives
    Age_pos=Age(RSVpos==0);
    RSVlog_pos=RSVlog(RSVpos==0);

    % preparing for fitting
    [data2]=Regroup(Age_pos,RSVlog_pos,1.5);
    Age_2=(data2(:,1))/12;% age in years
    Nr_2=data2(:,2);% Numbers at every age
    Rr_2=data2(:,3);% Numbers seropositive at every age
% 3. Minus LRTIs 
    Age_lrti=Age(LRTI==0);
    RSVlog_lrti=RSVlog(LRTI==0);

    % preparing for fitting
    [data3]=Regroup(Age_lrti,RSVlog_lrti,1.5);
    Age_3=(data3(:,1))/12;% age in years
    Nr_3=data3(:,2);% Numbers at every age
    Rr_3=data3(:,3);% Numbers seropositive at every age


%--------------------
% The grouped data
%--------------------
% 1. All the data
    %sorting the ages in ascending order and matching the respective titres
    [sortedage,IX]=sort(Age);sortedlog=RSVlog(IX);
    % Clustering the seropositives
    for k=1:20
        clear check
        check=(lower(k)<sortedage)&(sortedage<=upper(k));%which elements are in the kth group
        num(k)=nnz(check);%how many elements are in the kth group
        y=zeros(length(check),1);
        for j=1:length(check)
            if check(j)==1 && sortedlog(j)>=cutoff
               y(j)=1;
            end
        end
        pos(k)=nnz(y);        
    end
    Obs_1=pos./num;
    clear [sortedage,IX,sortedlog]
% 2. Minus RSV positives
    %sorting the ages in ascending order and matching the respective titres
    [sortedage,IX]=sort(Age_pos);sortedlog=RSVlog_pos(IX);
    % Clustering the seropositives
    for k=1:20
        clear check
        check=(lower(k)<sortedage)&(sortedage<=upper(k));%which elements are in the kth group
        num(k)=nnz(check);%how many elements are in the kth group
        y=zeros(length(check),1);
        for j=1:length(check)
            if check(j)==1 && sortedlog(j)>=cutoff
               y(j)=1;
            end
        end
        pos(k)=nnz(y);        
    end
    Obs_2=pos./num;
    clear [sortedage,IX,sortedlog]    
% 3. Minus LRTIs
    %sorting the ages in ascending order and matching the respective titres
    [sortedage,IX]=sort(Age_lrti);sortedlog=RSVlog_lrti(IX);
    % Clustering the seropositives
    for k=1:20
        clear check
        check=(lower(k)<sortedage)&(sortedage<=upper(k));%which elements are in the kth group
        num(k)=nnz(check);%how many elements are in the kth group
        y=zeros(length(check),1);
        for j=1:length(check)
            if check(j)==1 && sortedlog(j)>=cutoff
               y(j)=1;
            end
        end
        pos(k)=nnz(y);        
    end
    Obs_3=pos./num;
    clear [sortedage,IX,sortedlog]

% Plotting the grouped data
figure;scatter(age,Obs_1,'ro','filled');hold on;plot(age,Obs_1,'r')
hold on;scatter(age,Obs_2,'b.');hold on;plot(age,Obs_2,'b-')
hold on;scatter(age,Obs_3,'g*');hold on;plot(age,Obs_3,'g-')

%__________________________________________________________________________
%_______________________Estimating STEPWISE FUNCTION_______________________

% 1. All the data(just MSIR to check)
    % MSIR
    %Estimating
    x=[1.789,0.0035,2.1791]; % Initial guess
    p=0;%proportion that lose primary seropositive status
    n=2;% number on M classes
    lb=[0,0,0];
    options=optimset('MaxFunEvals',1000*length(x));
    [opt1,fval1,exitflag1]=fmincon(@msirsirmax15,x,[],[],[],[],lb,[],[],options,Age_1,Nr_1,Rr_1,delta,p,n);

    %Solving the ODEs with the new opt
    y0=zeros(n+4,1);y0(1)=0.999;y0(n+1)=0.001;
    [a1,y1]=ode45(@msirsir15,age,y0,[],opt1,delta,p,n);
    F1=sum(y1(:,1:n),2)+y1(:,n+2)+y1(:,n+3);

% 2. Minus RSV positives
    % MSIR
    %Estimating
    p=0;n=2;
    [opt2,fval2,exitflag2]=fmincon(@msirsirmax15,x,[],[],[],[],lb,[],[],options,Age_2,Nr_2,Rr_2,delta,p,n);
    %Solving the ODEs with the new opt
    y0=zeros(n+4,1);y0(1)=0.999;y0(n+1)=0.001;
    [a2,y2]=ode45(@msirsir15,age,y0,[],opt2,delta,p,n);
    F2_1=sum(y2(:,1:n),2)+y2(:,n+2)+y2(:,n+3);

    % MSIRSIR_1
    %Estimating
    p=1;n=2;
    [opt22,fval22,exitflag22]=fmincon(@msirsirmax15,x,[],[],[],[],lb,[],[],options,Age_2,Nr_2,Rr_2,delta,p,n);
    %Solving the ODEs with the new opt
    %The ODEs
    y0=zeros(n+4,1);y0(1)=0.999;y0(n+1)=0.001;
    [a22,y22]=ode45(@msirsir15,age,y0,[],opt22,delta,p,n);
    F2_2=sum(y22(:,1:n),2)+y22(:,n+2)+y22(:,n+3);

% 3. Minus LRTIs
    % MSIR
    %Estimating
    p=0;n=2;
    [opt3,fval3,exitflag3]=fmincon(@msirsirmax15,x,[],[],[],[],lb,[],[],options,Age_3,Nr_3,Rr_3,delta,p,n);
    %Solving the ODEs with the new opt
    y0=zeros(n+4,1);y0(1)=0.999;y0(n+1)=0.001;
    [a3,y3]=ode45(@msirsir15,age,y0,[],opt3,delta,p,n);
    F3_1=sum(y3(:,1:n),2)+y3(:,n+2)+y3(:,n+3);

    % MSIRSIR_1
    %Estimating
    p=1;n=2;
    [opt32,fval32,exitflag32]=fmincon(@msirsirmax15,x,[],[],[],[],lb,[],[],options,Age_3,Nr_3,Rr_3,delta,p,n);
    %Solving the ODEs with the new opt
    %The ODEs
    y0=zeros(n+4,1);y0(1)=0.999;y0(n+1)=0.001;
    [a32,y32]=ode45(@msirsir15,age,y0,[],opt32,delta,p,n);
    F3_2=sum(y32(:,1:n),2)+y32(:,n+2)+y32(:,n+3);

%__________________________________________________________________________
%________________Confidence region for the MSIR(yet to make it work)____________________________

% 1. All the data
    p=0;n=2;
    Paras_1=ConfidenceInt2_2([Age,RSVlog],Num,cutoff,n,x,p);
    % The Conidence regions of the fit
    Fci=zeros(20,Num);
    for k=1:Num
          %solving with the optimised values
        y0=zeros(n+4,1);y0(1)=0.999;y0(n+1)=0.001;
        [a,y]=ode45(@msirsir15,age,y0,options,Paras_1(k,:),delta,p,n);
        Fci(:,k)=sum(y(:,1:n),2)+y(:,n+2)+y(:,n+3);
    end
    figure;plot(age,Fci,'Color',[0.8,0.8,0.8]);xlabel('Age in years');ylabel('Proportions Seropositive')
    hold on;scatter(age,Obs_1,'ro','filled')
    hold on;plot(age,F1,'b')
% 2. Minus RSV positives
    p=0;n=1;
    Paras_2=ConfidenceInt2_2([Age_pos,RSVlog_pos],Num,cutoff,n,x,p);
    %Rate of loss of Mat Abs
    SIGMA_2=sort(Paras_2(:,1));
    CIsigma_2=quantile(SIGMA_2,[0.05,0.95]);
    %Force of infection for <=1  yr
    FOI1_2=sort(Paras_2(:,2));
    CIfoi1_2=quantile(FOI1_2,[0.05,0.95]);
    %Force of infection for >1  yr
    FOI2_2=sort(Paras_2(:,3));
    CIfoi2_2=quantile(FOI2_2,[0.05,0.95]);

% 3. Minus LRTIs
    p=0;n=1;
    Paras_3=ConfidenceInt2_2([Age_lrti,RSVlog_lrti],Num,cutoff,n,x,p);
    %Rate of loss of Mat Abs
    SIGMA_3=sort(Paras_3(:,1));
    CIsigma_3=quantile(SIGMA_3,[0.05,0.95]);
    %Force of infection for <=1  yr
    FOI1_3=sort(Paras_3(:,2));
    CIfoi1_3=quantile(FOI1_3,[0.05,0.95]);
    %Force of infection for >1  yr
    FOI2_3=sort(Paras_3(:,3));
    CIfoi2_3=quantile(FOI2_3,[0.05,0.95]);

%__________________________________________________________________________
%______________________________Plotting____________________________________

% Plotting the grouped data alone
figure;scatter(age,Obs_1,'ro','filled');hold on;plot(age,Obs_1,'r')
hold on;scatter(age,Obs_2,'b.');hold on;plot(age,Obs_2,'b-')
hold on;scatter(age,Obs_3,'g*');hold on;plot(age,Obs_3,'g-')


% Option 1: 3 separate figures
figure;scatter(age,Obs_1,'ro','filled');title('Using all the data')
hold on; plot(age,F1,'k')

figure;scatter(age,Obs_2,'bo','filled');title('Data minus RSV antigen positive')
hold on; plot(age,F2_1,'k')

figure;scatter(age,Obs_3,'go','filled');title('Data minus LTRI positive')
hold on; plot(age,F3_1,'k')

% Option 2:All 3 together 
figure;scatter(age,Obs_1,'ro','filled')
hold on;plot(age,F1,'r')

hold on;scatter(age,Obs_2,'b.')
hold on;plot(age,F2_1,'b')

hold on;scatter(age,Obs_3,'g*')
hold on;plot(age,F3_1,'g')

%____________________________________
% Plotting model predictions and FOI
M= sum(y1(:,1:n),2);
S= y1(:,n+1)+y1(:,n+4);
F= y1(:,n+2)+y1(:,n+3);

Pred=[M,F];
Foi= [ones(12,1)*opt1(2);ones(8,1)*opt1(3)];

figure;
[ax, h1] = plotyy(age, Pred, age,Foi, @bar, @line);
[ax, h1] = plotyy(age, Foi, age,Pred, @plot, @bar);

figure;
bar(Pred,'stacked')
figure;
plot(age,Foi)

