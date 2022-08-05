function y = proj_theta(x,thetamax)
    y=double(x>thetamax)*thetamax+double(x<thetamax).*x;
end