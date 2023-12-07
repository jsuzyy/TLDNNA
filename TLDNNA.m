function [BestCost,BestValue,XTarget]=TLDNNA(CostFunction,nPop,nVar,VarMin,VarMax,MaxIt,X)

% CostFunction:Objective function which you wish to minimize or maximize
% VarMin:Lower bound of a problem
% VarMax:Upper bound of a problem
% nVar:Number of design variables
% nPop:Population size
% MaxIt:Maximum number of iterations

% BestCost:Convergenc curve
% BestValue: The optimal objective function value
% XTarget:The optimal solution
%% --------------------Initialization----------------------------------------
X_LB=VarMin;
X_UB=VarMax;

beta=1;

% Creat random initial population
for i=1:nPop
    x_pattern(:,:,i)=X(:,:,i);
    cost(i)=CostFunction(x_pattern(:,:,i));
end

[COST,index]=max(cost);
%--------------------------------------------------------------------------
% Creat random initial weights with constraint of Summation each column = 1
ww=ones(1,nPop)*0.5;
w=diag(ww);
for i=1:nPop
    t=rand(1,nPop-1)*0.5;
    t=(t./sum(t))*0.5;
    w(w(:,i)==0,i)=t;
end
XTarget=[x_pattern(:,1,index);x_pattern(:,2,index)]';
Target=COST;                  % Best obtained objetive function value
wtarget=w(:,index);           % Best obtained weight (weight target)
BestCost(1)=Target;
%% -------------------- Main Loop for the NNA -------------------------------
tic

for i=1:nPop
    x_patternt(i,:)=[x_pattern(:,1,i);x_pattern(:,2,i)]';
end
x_pattern=x_patternt;
x_patternnew=x_pattern;
PBEST=x_pattern;
GBEST(1,:)=XTarget;
GBEST_COST(1)=Target;
for ii=2:MaxIt
    x_patternold=x_pattern;

    costold=cost;
    con1=0;
    con2=0;
    
    for i=1:nPop
        if rand<rand
            con1=con1+1;
            WININDEX(con1)=i;
        else
            con2=con2+1;
            LOSINDEX(con2)=i;
        end
        x_patternx(i,:)=1/3.*(x_pattern(i,:)+PBEST(i,:)+GBEST(randperm(length(GBEST_COST),1),:));
    end
    x_new=w(LOSINDEX)*x_patternx(LOSINDEX,:);
    x_pattern(LOSINDEX,:)=x_new+x_pattern(LOSINDEX,:);
    
    %------------------- Updating the weights -----------------------------
    for i=1:con2
        w(:,LOSINDEX(i))=abs(w(:,LOSINDEX(i))+((wtarget-w(:,LOSINDEX(i)))*2.*rand(nPop,1)));
    end
    
    for i=1:con2
        w(:,LOSINDEX(i))=w(:,LOSINDEX(i))./sum(w(:,LOSINDEX(i)));    % Summation of each column = 1
    end
 
    %----------------------- Creat new input solutions --------------------
    for i=1:nPop
        if rand<beta
            
            %------------- Bias for input solutions -----------------------
            N_Rotate=ceil(beta*nVar);
            
            xx=VarMin+(VarMax-VarMin).*rand(1,nVar);
            rotate_postion=randperm(nVar);
            rotate_postion=rotate_postion(1:N_Rotate);
            
            for m=1:N_Rotate
                x_pattern(i,rotate_postion(m))=xx(m);
            end
            %---------- Bias for weights ----------------------------------
            N_wRotate=ceil(beta*nPop);
            
            w_new=rand(N_wRotate,nPop);
            rotate_position=randperm(nPop);rotate_position=rotate_position(1:N_wRotate);
            
            for j=1:N_wRotate
                w(rotate_position(j),:)=w_new(j,:);
            end
            
            for iii=1:nPop
                w(:,iii)=w(:,iii)./sum(w(:,iii));   % Summation of each column = 1
            end
        else
            %------------ Transfer Function Operator
            %----------------------zz
            a=randperm(nPop,1);
            b=randperm(nPop,1);
        
            while a==b 
                a=randperm(nPop,1);
                b=randperm(nPop,1);
      
            end
            if costold(a)<costold(b)
                xtemp1=x_patternold(b,:)-x_patternold(a,:);
            else
                xtemp1=x_patternold(a,:)-x_patternold(b,:);
            end
   
            
            xtemp1=2.*rand(1,nVar).*xtemp1;
            
            fi=rand;
            
            if fi<=1/3
                x_pattern(i,:)=x_pattern(i,:)+2.*(XTarget-x_pattern(i,:)).*(rand(1,nVar))+xtemp1;
            elseif fi>=2/3
                x_pattern(i,:)=x_pattern(i,:)+2.*(PBEST(i,:)-x_pattern(i,:)).*(rand(1,nVar))+xtemp1;
            else
                x_pattern(i,:)=x_pattern(i,:)+2.*(GBEST(randperm(length(GBEST_COST),1),:)-x_pattern(i,:)).*(rand(1,nVar))+xtemp1;
            end
     
           
            
            x_pattern(i,:)=max(x_pattern(i,:),X_LB);    x_pattern(i,:)=min(x_pattern(i,:),X_UB);
        end
    end
    % ---------------------- Bias Reduction -------------------------------
    
    beta=beta*0.99;
    
    %----------------------------------------------------------------------
    %  x_pattern=max(x_pattern,X_LB);    x_pattern=min(x_pattern,X_UB);     % Check the side constraints
    %-------------- Calculating objective function values -----------------
    for i=1:nPop
        xx=x_pattern(i,:);
        xt=[xx(1:20)' xx(21:40)'];
        cost(i)=CostFunction(xt);
        if cost(i)<costold(i)
            cost(i)=costold(i);
            x_pattern(i,:)=x_patternold(i,:);
        else
            PBEST(i,:)=xx;
        end
    end
    %% ------ Selection ---------------------------------------------------
    [FF,Index]=max(cost);

    if FF>Target
        Target=FF;
        XTarget=x_pattern(Index,:);
        wtarget=w(:,Index);
        ml=size(GBEST,1);
        if ml<=nPop
            GBEST(ml,:)=XTarget;
            GBEST_COST(ml)=FF;
        else
            [~,lo]=min(GBEST_COST);
            GBEST(lo,:)=x_pattern(Index,:);
            GBEST_COST(lo)=FF;
        end
    else
        [~,Indexx]=min(cost);
        x_pattern(Indexx,:)=XTarget;
        w(:,Indexx)=wtarget;
    end
    
    x_patternnew=x_pattern;
    BestCost(ii)=Target;
   % disp(['Iteration ' num2str(ii) ': Best Cost = ' num2str(Target)]);
end

%% -------------------------------- NNA Finishes ----------------------------
BestValue=Target;%Output the optimal solution
end