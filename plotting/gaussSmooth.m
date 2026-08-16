function g = gaussSmooth(v,sigma)
 
v = v(:);

len = 2*round(sigma*4)+1;

f = gausswin(len,len/sigma);

if len > length(v)
    if mod(length(v),2)
        len = length(v);
    else
        len = length(v)-1;
    end
    
    f = f((length(f)+1)/2-len+1:(length(f)+1)/2+len-1);
end

f = f/sum(f);

g = conv(v,f);

c = cumsum(f);

g(1:length(c)) = g(1:length(c))./c;
g(end-length(c)+1:end) = g(end-length(c)+1:end) ./ flipud(c);

g = g((len-1)/2+1:end-(len-1)/2);
