clc;
clear;
close all;
MaxIt=1000;          
nPop=100;          
model = CreateModel_Case;
CostFunction=@(x) MyCost(x,model);   
nVar = model.n;      
VarSize=[nVar 2];  
VarMin=-model.MRANGE;          
VarMax = model.MRANGE;         
for i=1:30
    for k=1:nPop
        X(:,:,k)=CreateRandomSolution(model);
    end
    [BestCost(i,:,1),BestValue(1,i),Best(i,:,1)]=TLDNNA(CostFunction,nPop,2.*nVar,VarMin.*ones(1,2*nVar),VarMax.*ones(1,2*nVar),MaxIt,X);
    disp(['MNNA     Iteration ' num2str(i) ': Mean Best Cost = ' num2str(mean(BestValue15_1(1,i)),15) ': Best Cost = ' num2str(max(BestValue15_1(1,:)),15)]);
    disp(['**************************************************************************************************************************************************']);
end
