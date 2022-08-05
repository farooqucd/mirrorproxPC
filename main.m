clc
clear all
close all
warning('off','all')
rng(1)
profile on

%Consider a rectangular area with DxD m^2
%nAPs distributed APs serves nUsers terminals, they all randomly located in the
%area

nAPs=150; %number of APs
nUsers=15; %number of terminals
nTx=1; %number of antennas/AP
B=20; %bandwidth in Mhz

Tc=200; %coherence samples length
Tp=20; %pilot samples lenght
D=1; %area length in kilometer
d0=0.01; %reference distance in km
d1=0.05; %reference distance in km

Hb = 15; %base station height in m
Hm = 1.65; %mobile height in m
f = 1900; %frequency in MHz
aL = (1.1*log10(f)-0.7)*Hm-(1.56*log10(f)-0.8);
L = 46.3+33.9*log10(f)-13.82*log10(Hb)-aL; %path loss
sigma_sh=8; %shadowing in dB

noise_p = 10^((-203.975+10*log10(B*10^6)+9)/10); %noise power
power_p=0.2; %pilot power
zeta_p=power_p; %power of each pilot (not normalized)

myeps=1e-8; %stopping criteria
itt=100;

opts=sdpsettings('solver','mosek','verbose',0,'dualize',0); %options for YALMIP solver
maxBiter=20000; %max iterations for bisection
maxMPiter=20000; %max iterations for bisection
maxAPGiter = 20000; %max iteations for APG
myalpha = 0.2; %step size for APG
mytau = 100; %smoothing parameter

myomega=nAPs;
eta_max=myomega*0.2; %maximum power allocation of each user
theta_max=log(eta_max );


%% Pilot Asignment: (random choice)
[U,S,V11]=svd(randn(Tp,Tp)); %U includes Tp orthogonal sequences
%U=ones(Tp,Tp);
myPsi=zeros(Tp,nUsers); %pilot set of cell-free systems
if Tp<nUsers
    myPsi(:,[1:1:Tp])=U;
    for k=(Tp+1):nUsers
        Point=randi([1,Tp]);
        myPsi(:,k)=U(:,Point);
    end
else
    myPsi=U(:,[1:1:nUsers]);
end

%% Randomly enerate feasible power allocation (without normalization)
% For this case, the constraint is 0<=eta_k<=eta_max
% myEta = rand(nUsers,1)*eta_max;
Uopt=ones(nAPs,nUsers);
%% Large-scale fading matrix
myBetaBar=get_slow_fading(nAPs,nUsers,L,D,d0,d1,sigma_sh); %large-scale fading

%% Channel normalization wrt to the noise power
myBeta  = myBetaBar/(noise_p); %normalized fading channel

%% Channel estimaton
mau  = zeros(nAPs,nUsers);
for m=1:nAPs
    for k=1:nUsers
        mau(m,k)=norm( (myBeta(m,:).^(1/2)).*(myPsi(:,k)'*myPsi))^2;
    end
end
myNu=Tp*zeta_p*(myBeta.^2)./(Tp*zeta_p*mau + 1); %second moment of channel estimate

%% Precalculation
aki=zeros(nUsers,nUsers);
bki=zeros(nUsers,nUsers);
ck=zeros(nUsers,1);
for kUser=1:nUsers
    aki(kUser,kUser)=(Uopt(:,kUser)'*myNu(:,kUser))^2;
    for iUser=1:nUsers
        if(iUser~=kUser)
            myNutilde=abs(myPsi(:,kUser)'*myPsi(:,iUser))*myNu(:,kUser)./myBeta(:,kUser).*myBeta(:,iUser);
            aki(kUser,iUser)=(Uopt(:,kUser)'*myNutilde)^2;
        end
        Dki=myNu(:,kUser).*myBeta(:,iUser);
        bki(kUser,iUser)=(1/nTx)*(Uopt(:,kUser).^2)'*Dki;
    end
    ck(kUser) = (1/nTx)*(Uopt(:,kUser).^2)'*myNu(:,kUser);
    ck(kUser)=myomega*ck(kUser)/aki(kUser,kUser);
    bki(kUser,:)=bki(kUser,:)/aki(kUser,kUser);
    aki(kUser,:)=aki(kUser,:)/aki(kUser,kUser);
end

%% Solve convex power allocation using Bisection in CVX
tmin=0;
tmax=max(eta_max/ck);
etaB=zeros(nUsers,1);
ti = cputime;
for iB=1:maxBiter
    t0=(tmin+tmax)/2;
    etaY=sdpvar(nUsers,1);
    obj=[];
    F=[etaY>=0,etaY(kUser)<=eta_max];
    for kUser=1:nUsers
        jUser=1:nUsers~=kUser;
        F=[F,(aki(kUser,jUser)*etaY(jUser)+bki(kUser,:)*etaY+ck(kUser))*t0<=etaY(kUser)];
    end
    diagnotics = optimize(F,-obj,opts);

    if(diagnotics.problem==0)
        tmin=t0;
        etaB=value(etaY);
    else
        tmax=t0;
    end
    if(tmax-tmin<myeps)
        break;
    end
end
tf = cputime;
objB=uplink_userRate(nTx,nUsers,myNu,myBeta,myPsi,etaB,Uopt,myomega);
tB = tf-ti;

%% Accelerated Projected Gradient
x_prev = rand(nUsers,1);
t_now = 1;
t_prev = 1;
x_now = x_prev;
prevcobj = 0;
bestcobj=zeros(maxAPGiter,1);
ti = cputime;
for itercAPG=1:maxAPGiter
    y = x_now+((t_prev-1)/t_now)*(x_now-x_prev);
    Dy = gradConvexPC(nUsers,y,aki,bki,ck,mytau);
    y_next= y-myalpha*Dy;

    x_now = proj_hyperplane(y_next,-Inf,theta_max);

    tempobj = uplink_userRate(nTx,nUsers,myNu,myBeta,myPsi,exp(x_now),Uopt,myomega);
    if itercAPG==1
        bestcobj(itercAPG)=tempobj;
    else
        bestcobj(itercAPG)=max(bestcobj(itercAPG-1),tempobj);
    end
    x_prev = x_now;
    t_prev = t_now;
    t_now = (sqrt(4*t_now^2+1)+1)/2;
    prevcobj=bestcobj(itercAPG);
    if(itercAPG>itt && abs((bestcobj(itercAPG)-bestcobj(itercAPG-itt))/bestcobj(itercAPG))<myeps)
        break;
    end
end
tf = cputime;
bestcobj(itercAPG+1:maxAPGiter)=bestcobj(itercAPG);
objAPG(:,1)=bestcobj;
tAPG = tf-ti;

%% Mirror Prox Algorithm
theta_prev=rand(nUsers,1);
lambda_prev=rand(nUsers,1);
stepsize=2;
rr=0.5;
objMPiter=zeros(maxMPiter,1);
tt = 0;
for iMPiter=1:maxMPiter
    tt0 = cputime;
    [Dtheta_temp,Dlambda_temp]=gradminmaxboth(nUsers,aki,bki,ck,theta_prev,lambda_prev);
    for ii=1:5
        theta_temp_next=proj_theta(theta_prev-(stepsize)*Dtheta_temp,theta_max);
        lambda_temp_next=projsplx(lambda_prev+(stepsize)*Dlambda_temp);
        [Dtheta,Dlambda,~]=gradminmaxboth(nUsers,aki,bki,ck,theta_temp_next,lambda_temp_next);
        theta_next=proj_theta(theta_prev-(stepsize)*Dtheta,theta_max);
        lambda_next=projsplx(lambda_prev+(stepsize)*Dlambda);
        delta=stepsize*((Dtheta)'*(theta_temp_next-theta_next)) ...
            -((Dlambda)'*(lambda_temp_next-lambda_next)) ...
            -0.5*norm(theta_next-theta_prev)^2 ...
            -0.5*norm(lambda_next-lambda_prev)^2;
        if delta<=0
            break;
        else
            stepsize=stepsize*rr;
        end
    end
    if(mod(iMPiter,101)==0)
        stepsize = 1;
    end

    tt1 = cputime;
    tempobj=uplink_userRate(nTx,nUsers,myNu,myBeta,myPsi,exp(theta_next),Uopt,myomega);
    if iMPiter==1
        objMPiter(iMPiter)=tempobj;
    else
        objMPiter(iMPiter)=max([objMPiter(iMPiter-1);tempobj]);

    end
    tt2 = cputime;

    theta_prev=theta_next;
    lambda_prev=lambda_next;
    if(iMPiter>itt && abs((objMPiter(iMPiter)-objMPiter(iMPiter-itt))/objMPiter(iMPiter))<myeps)
        break;
    end
    tt3 = cputime;
    tt = tt+tt1-tt0 + tt3-tt2;
end
objMPiter(iMPiter+1:maxMPiter)=objMPiter(iMPiter);
objMP(:,1)=objMPiter;
tMP = tt;

profile viewer
%% Plot Objectives
figure
plotIter=3000;
iterationsB=objB*ones(plotIter,1); %take a few iterations
iterationsMP=objMP(1:plotIter); %take a few iterations
iterationsAPG=objAPG(1:plotIter,:); %take a few iterations

iterationsB(end) %bisection final objective
iterationsMP(end) %MP final objective
iterationsAPG(end) %APG final objective

plot(1:plotIter,iterationsB,'g'); %plot bisection iterations
hold on
plot(1:plotIter,iterationsMP,'b'); %plot MP iterations
plot(1:plotIter,iterationsAPG,'k'); %plot APG iterations
hold off
