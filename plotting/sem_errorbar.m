function [mn,sem,yu,yl] = sem_errorbar(X)

mn = mean(X,1); 
sem = std(X)/sqrt(size(X,1));
yu = mn+sem; yl = mn-sem;

end