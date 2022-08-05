function x = projsplx(y)
x = max(y-max((cumsum(sort(y,1,'descend'),1)-1)./(1:size(y,1))'),0);
return;