function out = proj_hyperplane(x,lb,ub)
% out=zeros(length(x),1);
% for iUser=1:length(x)
%     if(x(iUser)<lb(iUser))
%         out(iUser)=lb(iUser);
%     elseif(x(iUser)>ub(iUser))
%         out(iUser)=ub(iUser);
%     else
%         out(iUser)=x(iUser);
%     end
% end 
if lb==-Inf && ub==Inf
    out=x;
elseif lb==-Inf
    out=double(x>ub).*ub+double(x>lb&x<ub).*x;
elseif ub==Inf
    out=double(x<lb).*lb+double(x>lb&x<ub).*x;
else
    out=double(x<lb).*lb+double(x>ub).*ub+double(x>lb&x<ub).*x;
end