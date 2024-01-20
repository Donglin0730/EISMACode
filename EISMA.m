function [Best_pos,Best_score,curve]=EISMA(pop,Max_iter,lb,ub,dim,fobj)
	zz=0.03; 
if(max(size(ub)) == 1)
    ub = ub.*ones(1,dim);
    lb = lb.*ones(1,dim);
end

X0=initialization(pop,dim,ub,lb);

fitness = zeros(1,pop);
for i = 1:pop
    fitness(i)=fobj(X0(i,:));
end
[fitness, index]= sort(fitness);
GBestF = fitness(1);
for i = 1:pop
    X(i,:) = X0(index(i),:);
end
GBestX = X(1,:);
curve=zeros(1,Max_iter);
W = zeros(pop,dim);
for t = 1: Max_iter
    worstFitness = fitness(end);
    bestFitness = fitness(1);
    S=bestFitness-worstFitness+eps;  
    for i = 1: pop
        if i<pop/2 
            W(i,:)= 1+rand(1,dim).*log10((bestFitness-fitness(i))/(S)+1);
        else 
            W(i,:)= 1-rand(1,dim).* log10((bestFitness-fitness(i))/(S)+1);
        end
    end
    a = atanh(-(t/Max_iter)+1);
    b = 1-t/Max_iter;
    for i=1:pop
        Xnew1 =GBestX +a*Levy(dim).*(mean(X) - X(i,:));
        if fobj(Xnew1)< fobj(X(i,:))
            flag =1;   
        else
             flag =-1; 
        end
        if rand<zz
            X(i,:) = (ub-lb)*rand+lb;
        else
            p =tanh(abs(fitness(i)-GBestF));
            vb = unifrnd(-a,a,1,dim);
            vc = unifrnd(-b,b,1,dim);
            for j=1:dim
                r = rand();
                A = randi([1,pop]); 
                B = randi([1,pop]);
                if r<p
                    if flag == 1
                        X(i,j) =Xnew1(j);
                    else
                        X(i,j) = GBestX(j)+ vb(j)*(W(i,j)*X(A,j)-X(B,j));
                    end
                else
                    X(i,j) = vc(j)* X(i,j);
                end
            end
        end
    end
    for j = 1:pop
        for oo = 1: dim
            if(X(j,oo)>ub(oo))
                X(j,oo) =ub(oo);
            end
            if(X(j,oo)<lb(oo))
                X(j,oo) =lb(oo);
            end
        end
    end
    for j=1:pop
        fitness(j)=fobj(X(j,:));
    end
    fmax = 1; fmin = 0.25;
    w = (fmax - fmin)*sin(((t/Max_iter)^3)*pi/2) + fmin;
    for it = 1:pop
        if rand<0.5
            index = randperm(pop,3);
            Stemp = X(it,:) + rand.*(X(index(1),:)-X(it,:)) + rand.*(X(index(2),:)-X(index(3),:));
        else
            index = randperm(pop,3);
            Stemp = GBestX +w.*(X(index(1),:)-X(it,:)) +w.*((X(index(2),:)-X(index(3),:)));
        end
        Stemp(Stemp>ub) = ub(Stemp>ub);
        Stemp(Stemp<lb) = lb(Stemp<lb);
        fitvnew =fobj(Stemp);
        fitv = fobj(X(it,:));
        if(fitvnew <fitv )
            X(it,:) =  Stemp;
        end
    end
    for j=1:pop
        fitness(j)=fobj(X(j,:));
    end
    [fitness, index]= sort(fitness);
    for j = 1:pop
        X(j,:) = X(index(j),:);
    end
    XMb=(X(1,:)+X(2,:)+X(3,:))/3;
    S2 = XMb+XMb.* (randn(1,dim)/norm(randn(1,dim)));
    S2 = rand.*X(1,:)+(1-rand).*S2;
    S2(S2>ub) = ub(S2>ub);
    S2(S2<lb) = lb(S2<lb);
    if  fobj(S2)<fobj(X(1,:))
        X(1,:) = S2;
    end    
      if(fobj(X(1,:))<GBestF)
        GBestF=fobj(X(1,:))
        GBestX = X(1,:);
    end
    curve(t) = GBestF;
end
Best_pos = GBestX;
Best_score = curve(end);
end

function o=Levy(d)
beta=1.5;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
u=randn(1,d)*sigma;v=randn(1,d);step=u./abs(v).^(1/beta);
o=step;
end
