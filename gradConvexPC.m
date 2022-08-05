function [out] = gradConvexPC(nUsers,mytheta,aki,bki,ck,mytau)
ematrix=eye(nUsers);
invSINR=zeros(1,nUsers);
gradinvSINR=zeros(nUsers,nUsers);
for kUser=1:nUsers
    jUser=1:nUsers~=kUser;
    aterm=aki(kUser,jUser).*exp(mytheta(jUser)-mytheta(kUser))';
    BU=sum(aterm);
    dBU=sum(aterm.*(ematrix(:,jUser)-ematrix(:,kUser)),2);
    bterm=bki(kUser,:).*exp(mytheta-mytheta(kUser))';
    IUI=sum(bterm);
    dIUI=sum(bterm.*(ematrix-ematrix(:,kUser)),2);
    TN=ck(kUser)*exp(-mytheta(kUser));
    dTN=-TN*ematrix(:,kUser);

    % calculate inverse SINR and its gradient
    invSINR(kUser)=BU+IUI+TN;
    gradinvSINR(:,kUser)=dBU+dIUI+dTN;
end
out = sum(exp(mytau*invSINR).*gradinvSINR/sum(exp(mytau*invSINR)),2);

