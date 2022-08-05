function [out,invSINR] = uplinkConvexPC(nUsers,mytheta,aki,bki,ck)
invSINR=zeros(nUsers,1);
for kUser=1:nUsers
    jUser=1:nUsers~=kUser;
    BU=aki(kUser,jUser)*exp(mytheta(jUser)-mytheta(kUser));
    IUI=bki(kUser,:)*exp(mytheta-mytheta(kUser));
    TN=ck(kUser)*exp(-mytheta(kUser));
    % calculate inverse SINR
    invSINR(kUser)=BU+IUI+TN;
end
out=max(invSINR);