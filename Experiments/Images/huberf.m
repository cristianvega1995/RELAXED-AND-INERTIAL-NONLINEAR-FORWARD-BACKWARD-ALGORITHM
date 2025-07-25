function f = huberf(x,lambda,mu)

f = lambda*sum(sum( (abs(x)-mu/2).*(abs(x)>mu)+(abs(x).^2)/(2*mu).*(abs(x)<=mu)));

