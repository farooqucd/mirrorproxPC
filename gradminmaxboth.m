function [gradtheta,gradlambda,minrate] = gradminmaxboth(nUsers,aki,bki,ck,mytheta,mylambda)
gradtheta=zeros(nUsers,1);
gradlambda=zeros(nUsers,1);
for kUser=1:nUsers
    bterm=bki(kUser,:).*exp(mytheta-mytheta(kUser))';
    d=bterm.';
    TN=ck(kUser)*exp(-mytheta(kUser));
    IUI=sum(bterm);
    d(kUser)=d(kUser)-TN-IUI;
    gradlambda(kUser)=IUI+TN;
    gradtheta = gradtheta + d*mylambda(kUser);
end
SINR = 1./gradlambda;
minrate = log2(1+min(SINR));
