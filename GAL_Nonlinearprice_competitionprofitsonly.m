%%Genetic Algorith Loop--Nonlinear in prices
clear
n=60; %number trials array=(n,N,t)
N=12;  %number of agents or length of string
One=ones(1,N);%creates a 1xN matrix of ones
t=100;    %number of generations
y=repmat(N+1, n,1);
probmut=.001;%probability mutation will occur
a=5; %multiplier of transport costs
b=1; %power price is raised to
c=.5; %power transport costs raised to
d=1;
e=repmat(N^d,n,1); %haven't added d&e to solutions yet!!!
solution=0; %if solution is 1, run genetic algorithm plugging in the solution, if zero, don't run algorithm on solution


%Generating First Generation values

 %firstgen=binornd(ones(n,N),.5); %chose binomial because this gives us a 1 or 0 each with probability of .5.
 bob=rand(n,N);
 firstgen=zeros(n,N);
 firstgen(find(bob<=.5))=1;
 
    %%Fitness value first generation
mktnum=sum(firstgen,2); %sums number of agents in market in each trial. 2=sums across columns. this is an array: (n,1,t)
price=(y-(mktnum.^d)+e).^b;%calculates price in market. did a non linear price calculation=(1+number of agents-the number in the market for each string)^2.
% sizef=size(firstgen);
% [rowf,columnf]=ind2sub(sizef,find(firstgen==1));   
% indfitness=zeros(n,N);
% indfitness(rowf,columnf)=price(rowf,:)-a*(columnf).^c;


    for i=1:n
    for j=1:N  
        if firstgen(i,j)==0; %if the agent is not in the market
            indfitness(i,j)=0;%profits when not in the market
        else
            indfitness(i,j)=price(i,1)-a*j^c; %profit from staying in market
        end
    end  %if statement describes what fitness is used for each individual. if they're in the market we calculate their profit if out of the market we calculate the opportunity cost of not being in the market. 
    end
    fitness=sum(indfitness,2)+(mktnum.^2)./2; %summing across columns to get trial fitness, nx1xt
    minfit=min(fitness(:,:,1), [], 1); %finds minimum value for fitness.  The one denotes finding the min of all of the row values. min value may be negative so we must map to positive numbers
    if minfit<=0
        posfitness=fitness+(-1*minfit)*ones(n,1)+ones(n,1); %if the minimum fitness is negative, add the absolute value plus one to all fitness values.  (adding one gets rid of zeros)
    else
        posfitness=fitness;
    end
 

    %Calculating fitness proportion
    totalfitness=sum(posfitness);  %to get the total we must sum all string fitness values together. 
    fitnessprop=posfitness./totalfitness; %finds proportion of fitness for each trial by dividing each string fitness value by the total fitness.

    %New Parent Generation
    Parents=randsample(1:n,n,true,fitnessprop'); %randomly generates a number 1-n depending on the corresponding probabilities
    %randsample(1:n,n,true,fitnessprop') chose to you this because its the
    %only way to sample with replacement from a probability distribution of
    %your choice.
    newparents=zeros(n,N);
    for k=1:n
        newparents(k,:)=firstgen(Parents(1,k),:);  %taking the randomly generated numbers and puting together a matrix with the selected parents
    end

    %Matching In this process I reorganize the new parent matrix randomly
    %and then group strings together in pairs of two.
    shuffledparents=newparents(randperm(size(newparents,1)),:); %shuffles the newparent matrix around in preparation for matching/crossover
    %randperm selects rows and randomly reorganizes them.

    %Crossover
    for r=1:2:n %gives values for r so that we only do crossover in each pair once. r goes up by two so that the pair r,r+1 is unique. 
        stop=randsample(N,1);%generates randomly a number in 1 through N in which crossover will take place for the r,r+1 pair
       
        newgen(r:r+1,:)=[shuffledparents(r,1:stop) shuffledparents(r+1,stop+1:N); shuffledparents(r+1,1:stop) shuffledparents(r,stop+1:N)];
        %generates a new matrix in which for each pair there is crossover
        %at the randomly chosen stopping point for that pair.  if stop=2,
        %the first two columns of row r and r+1 remain the same.  columns 3
        %through N swap.
    end


%Generations 2-t


newgen_g(:,:,1)=newgen;
%newgen_g(:,:,1)=[ones(n,9) zeros(n,N-9)];

for g=2:t
    
 %%Fitness value
    mktnum_g=sum(newgen_g,2); %sums number of agents in market in each trial this is an array: (n,1,t)
    price_g=(y-(mktnum_g(:,:,g-1).^d)+e).^b;%calculates price in market. did a linear price calculation=1+number of agents-the number in the market for each string.
    priceplus_g=(y-(mktnum_g(:,:,g-1).^d)+e-(d*mktnum_g(:,:,g-1)+ones(n,1))).^b;%price calculated for oppt cost calculations
    %used repmat(N+1,[n 1 t]) because I needed an array of the number N+1
    %added ones(n,1,t) to price because we want to get the price if the agent
    %had entered the market.
    
 
    indfitness_g=zeros(n,N,t-1);
    for i=1:n
    for j=1:N
        if newgen_g(i,j,g-1)==0;%if the agent is not in the market
            indfitness_g(i,j,g-1)=0;%profit when not in market
        else
            indfitness_g(i,j,g-1)=price_g(i,1)-a*j^c; %profit from entering market
        end
    end
    end
    
    indfitness_g;
    fitness_g=sum(indfitness_g(:,:,g-1),2)+(mktnum_g(:,:,g-1).^2)./2; %summing across columns to get trial fitness, nx1xt
    minfit_g=min(fitness_g, [], 1) ; %calculates min fitness down each generation's fitness column
     if minfit_g<=0;  %fitness may be less than zero.  if so we add 1 to the absolute value of min fit and then add this to all fitness values for each trial.  If min fit isn't negative we leave everything the same.
        posfitness_g=fitness_g+(-1*minfit_g)*ones(n,1)+ones(n,1);
     else
        posfitness_g=fitness_g;
     end
    posfitness_g;
    %Calculating fitness proportion
    totalfitness_g=sum(posfitness_g,1);  %1x1
    fitnessprop_g=zeros(n,1,t-1);
    for i=1:n
    fitnessprop_g(i,1,g-1)=posfitness_g(i,:)./totalfitness_g;  %nx1xt
    end
    %finds proportion of fitness for each trial
    
    %New Parent Generation
    Parents_g=randsample(1:n,n,true,fitnessprop_g(:,:,g-1)'); %randomly generates a number 1-n depending on the corresponding probabilities
    %randsample(1:n,n,true,fitnessprop') chose to you this because its the only
    %way to sample with replacement from a probability distribution of your
    %choice.
    newparents_g=zeros(n,N);
    for k=1:n
        newparents_g(k,:)=newgen_g(Parents_g(1,k),:,g-1);  %taking the randomly generated numbers and puting together a matrix with the selected parents
    end

    %Matching
    %In this process I reorganize the new parent matrix randomly and then
    %group strings together in pairs of two.
    shuffledparents_g=newparents_g(randperm(size(newparents_g,1)),:); %shuffles the newparent matrix around in preparation for matching/crossover
    %randperm selects rows and randomly reorganizes them. 
    

    %Crossover
    for r=1:2:n %gives values for r so that we only do crossover in each pair once. r goes up by two so that the pair r,r+1 is unique. 
        stop=randsample(N,1);%generates randomly a number in 1 through N in which crossover will take place for the r,r+1 pair
        newgen_g(r:r+1,:,g)=[shuffledparents_g(r,1:stop) shuffledparents_g(r+1,stop+1:N); shuffledparents_g(r+1,1:stop) shuffledparents_g(r,stop+1:N)];
        
        %generates a new matrix in which for each pair there is crossover at the
        %randomly chosen stopping point for that pair.  if stop=2, the first two
        %columns of row r and r+1 remain the same.  columns 3 through N swap.  
    end
   
    
  %Mutation
    bob=rand(n,N);
    mutate_g=zeros(n,N);
    mutate_g(find(bob<=probmut))=1;
  
  for i=1:n
  for j=1:N
      if mutate_g(i,j)==1 & newgen_g(i,j,g)==0;  %if statement so that if there is a mutation of the agent, that agent newgen_g is changed to the opposite of what it was.  
          newgen_g(i,j,g)=mutate_g(i,j);  %if statement says, if mutation agent=1 and newgen agent=0 then newgen agent=1, if mutation agent =1 and newgen agent=1 then newgen agent=0, else (mutation agent=0) newgen stays the same.
      elseif  mutate_g(i,j)==1 & newgen_g(i,j,g)==1;
          newgen_g(i,j,g)=0;
      end
  end
  end
  
  
      
end
%Getting every permutation dealing with N=length of string
%generating a matrix with each row representing a possibility for the
%number of ones in each string

perm_ops=zeros(N,N);
fact=ones(N,1);
for i=1:N
    perm_ops(i,:)=[ones(1,i) zeros(1,N-i)];
    fact(i,1)=factorial(N)/(factorial(i)*factorial(N-i));
end
Perm=zeros(1,N);
for i=1:N
  C=nextperms(perm_ops(i,:), fact(i,1));
  Perm=[Perm; C'];
end

C=cell(N,1);

% %generating the unique permutations
% for i=1:N
%   [num x p]=uperms(perm_ops(i,:));
%   C{i,1}=p;
% end
% 
% Perm=cell2mat(C);
 m=size(Perm,1);
yp=repmat(N+1, m,1);
ep=repmat(N^d,m,1);
 mktnump=sum(Perm,2); %sums number of agents in market in each trial. 2=sums across columns. this is an array: (n,1,t)
pricep=(yp-(mktnump.^d)+ep).^b;%calculates price in market. did a non linear price calculation=(1+number of agents-the number in the market for each string)^2.priceplus=(y-mktnum-ones(n,1)).^3;
priceplusp=(yp-(mktnump.^d)+ep-(d.*mktnump+ones(m,1))).^b;

indfitnessp=zeros(m,N);
for i=1:m
    for j=1:N
        if Perm(i,j)==0; %if the agent is not in the market
            indfitnessp(i,j)=0;%profit when not in market
        else
            indfitnessp(i,j)=pricep(i,1)-a*j^c; %profit from entering market
        end
    end  %if statement describes what fitness is used for each individual. if they're in the market we calculate their profit if out of the market we calculate the opportunity cost of not being in the market. 
end
    fitnessp=sum(indfitnessp,2)+(mktnump.^2)./2; %summing across columns to get trial fitness, nx1xt
    minfitp=min(fitnessp(:,:,1), [], 1); %finds minimum value for fitness.  The one denotes finding the min of all of the row values. min value may be negative so we must map to positive numbers
    if minfitp<=0
        posfitnessp=fitnessp+(-1*minfitp)*ones(m,1)+ones(m,1); %if the minimum fitness is negative, add the absolute value plus one to all fitness values.  (adding one gets rid of zeros)
    else
        posfitnessp=fitnessp;
    end
[max_val, in]=max(posfitnessp);
Best=Perm(in,:)
plot(fitnessp, '.')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if solution==1;
firstgen=repmat(Best,n,1);


  %%Fitness value first generation
    mktnum=sum(firstgen,2); %sums number of agents in market in each trial. 2=sums across columns. this is an array: (n,1,t)
    price=(y-2*mktnum).^b;%calculates price in market. did a non linear price calculation=(1+number of agents-the number in the market for each string)^2.
    priceplus=(y-2*mktnum-ones(n,1)).^b;%price calculated for oppt cost calculations
    %used repmat(N+1,[n 1 t]) because I needed an array of the number N+1
    %subtracted ones(n,1) to price because we want to get the price if the agent
    %had entered the market.
    
    for i=1:n
    for j=1:N
        if firstgen(i,j)==0; %if the agent is not in the market
            indfitness(i,j)=0;%profits when not in market
           
        else
            indfitness(i,j)=price(i,1)-a*j^c; %profit from staying in market
        end
    end  %if statement describes what fitness is used for each individual. if they're in the market we calculate their profit if out of the market we calculate the opportunity cost of not being in the market. 
    end
    fitness=sum(indfitness,2)+(mktnum.^2)./2; %summing across columns to get trial fitness, nx1xt
    minfit=min(fitness(:,:,1), [], 1); %finds minimum value for fitness.  The one denotes finding the min of all of the row values. min value may be negative so we must map to positive numbers
    if minfit<=0
        posfitness=fitness+(-1*minfit)*ones(n,1)+ones(n,1); %if the minimum fitness is negative, add the absolute value plus one to all fitness values.  (adding one gets rid of zeros)
    else
        posfitness=fitness;
    end
 

    %Calculating fitness proportion
    totalfitness=sum(posfitness);  %to get the total we must sum all string fitness values together. 
    fitnessprop=posfitness./totalfitness; %finds proportion of fitness for each trial by dividing each string fitness value by the total fitness.

    %New Parent Generation
    Parents=randsample(1:n,n,true,fitnessprop'); %randomly generates a number 1-n depending on the corresponding probabilities
    %randsample(1:n,n,true,fitnessprop') chose to you this because its the
    %only way to sample with replacement from a probability distribution of
    %your choice.
    for k=1:n
        newparents(k,:)=firstgen(Parents(1,k),:);  %taking the randomly generated numbers and puting together a matrix with the selected parents
    end

    %Matching In this process I reorganize the new parent matrix randomly
    %and then group strings together in pairs of two.
    shuffledparents=newparents(randperm(size(newparents,1)),:); %shuffles the newparent matrix around in preparation for matching/crossover
    %randperm selects rows and randomly reorganizes them.

    %Crossover
    for r=1:2:n %gives values for r so that we only do crossover in each pair once. r goes up by two so that the pair r,r+1 is unique. 
        stop=randsample(N,1);%generates randomly a number in 1 through N in which crossover will take place for the r,r+1 pair
       
        newgen(r:r+1,:)=[shuffledparents(r,1:stop) shuffledparents(r+1,stop+1:N); shuffledparents(r+1,1:stop) shuffledparents(r,stop+1:N)];
        %generates a new matrix in which for each pair there is crossover
        %at the randomly chosen stopping point for that pair.  if stop=2,
        %the first two columns of row r and r+1 remain the same.  columns 3
        %through N swap.
    end


%Generations 2-t


newgen_g(:,:,1)=newgen;
%newgen_g(:,:,1)=[ones(n,9) zeros(n,N-9)];

for g=2:t
    
 %%Fitness value
    mktnum_g=sum(newgen_g,2); %sums number of agents in market in each trial this is an array: (n,1,t)
    price_g=(y-2*mktnum_g(:,:,g-1)).^b;%calculates price in market. did a linear price calculation=1+number of agents-the number in the market for each string.
    priceplus_g=(y-2*mktnum_g(:,:,g-1)-ones(n,1)).^b;%price calculated for oppt cost calculations
    %used repmat(N+1,[n 1 t]) because I needed an array of the number N+1
    %added ones(n,1,t) to price because we want to get the price if the agent
    %had entered the market.
    
 
    
    for i=1:n
    for j=1:N
        if newgen_g(i,j,g-1)==0;%if the agent is not in the market
            indfitness_g(i,j,g-1)=0;%profit when not in market
        else
            indfitness_g(i,j,g-1)=price_g(i,1)-a*j^c; %profit from entering market
        end
    end
    end
    
    indfitness_g;
    fitness_g=sum(indfitness_g(:,:,g-1),2)+(mktnum_g(:,:,g-1).^2)./2; %summing across columns to get trial fitness, nx1xt
    minfit_g=min(fitness_g, [], 1) ; %calculates min fitness down each generation's fitness column
     if minfit_g<=0;  %fitness may be less than zero.  if so we add 1 to the absolute value of min fit and then add this to all fitness values for each trial.  If min fit isn't negative we leave everything the same.
        posfitness_g=fitness_g+(-1*minfit_g)*ones(n,1)+ones(n,1);
     else
        posfitness_g=fitness_g;
     end
    posfitness_g;
    %Calculating fitness proportion
    totalfitness_g=sum(posfitness_g,1);  %1x1
    for i=1:n
    fitnessprop_g(i,1,g-1)=posfitness_g(i,:)./totalfitness_g;  %nx1xt
    end
    %finds proportion of fitness for each trial
    
    %New Parent Generation
    Parents_g=randsample(1:n,n,true,fitnessprop_g(:,:,g-1)'); %randomly generates a number 1-n depending on the corresponding probabilities
    %randsample(1:n,n,true,fitnessprop') chose to you this because its the only
    %way to sample with replacement from a probability distribution of your
    %choice.
    for k=1:n
        newparents_g(k,:)=newgen_g(Parents_g(1,k),:,g-1);  %taking the randomly generated numbers and puting together a matrix with the selected parents
    end

    %Matching
    %In this process I reorganize the new parent matrix randomly and then
    %group strings together in pairs of two.
    shuffledparents_g=newparents_g(randperm(size(newparents_g,1)),:); %shuffles the newparent matrix around in preparation for matching/crossover
    %randperm selects rows and randomly reorganizes them. 
    

    %Crossover
    for r=1:2:n %gives values for r so that we only do crossover in each pair once. r goes up by two so that the pair r,r+1 is unique. 
        stop=randsample(N,1);%generates randomly a number in 1 through N in which crossover will take place for the r,r+1 pair
        newgen_g(r:r+1,:,g)=[shuffledparents_g(r,1:stop) shuffledparents_g(r+1,stop+1:N); shuffledparents_g(r+1,1:stop) shuffledparents_g(r,stop+1:N)];
        
        %generates a new matrix in which for each pair there is crossover at the
        %randomly chosen stopping point for that pair.  if stop=2, the first two
        %columns of row r and r+1 remain the same.  columns 3 through N swap.  
    end
    
  %Mutation
    bob=rand(n,N);
    mutate_g=zeros(n,N);
    mutate_g(find(bob<=probmut))=1;
  
  for i=1:n
  for j=1:N
      if mutate_g(i,j)==1 & newgen_g(i,j,g)==0;  %if statement so that if there is a mutation of the agent, that agent newgen_g is changed to the opposite of what it was.  
          newgen_g(i,j,g)=mutate_g(i,j);  %if statement says, if mutation agent=1 and newgen agent=0 then newgen agent=1, if mutation agent =1 and newgen agent=1 then newgen agent=0, else (mutation agent=0) newgen stays the same.
      elseif  mutate_g(i,j)==1 & newgen_g(i,j,g)==1;
          newgen_g(i,j,g)=0;
      end
  end
  end
  
  
      
end
else
    newgen_g=newgen_g;
end

