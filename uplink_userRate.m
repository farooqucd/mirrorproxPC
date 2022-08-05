function [out,uRates,uSINR] = uplink_userRate(nTx,nUsers,myNu,myBeta,myPsi,myEta,RxCoeff,myomega)
uSINR=zeros(nUsers,1);
for kUser=1:nUsers
    sig=(RxCoeff(:,kUser)'*myNu(:,kUser))^2*myEta(kUser);
    Dki=myNu(:,kUser).*myBeta(:,1:nUsers);
    IUI=sum((1/nTx)*(RxCoeff(:,kUser).^2)'*Dki*myEta);
    TN = myomega*(1/nTx)*(RxCoeff(:,kUser).^2)'*myNu(:,kUser);
    interference=IUI+TN;
    uSINR(kUser)=sig/interference;
end
uRates=log2(1+uSINR);
out=min(uRates);
