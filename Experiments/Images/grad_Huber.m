function df = grad_Huber(x,mu)



%f = sum(sum( (abs(x)-mu/2).*(abs(x)>mu)+(abs(x).^2)/(2*mu).*(abs(x)<=mu)));

df = ((abs(x)>mu & (x)>0)-(abs(x)>mu & (x)<0)+(x./mu).*(abs(x)<=mu));



