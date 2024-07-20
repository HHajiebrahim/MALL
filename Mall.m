clc;
clear;
close all;
%% Problem definition
tic
CostFunction = @(x) f_3(x);          % Objective Function

FunctionType = 'MIN';                % Define function type as MIN of MAX


nVar = 30;                           % Number of decision variables
VarSize = [nVar,1];                  % Decision variable matrix size

VarMin = -100;                       % Lower bound of decision variables
VarMax = 100;                        % Upper bound of decision variables
%% MallAlgorithm Parameters

MaxIt = 300;                         % Maximum Number of Iterations

SellerPop = 4;                       % Sellers Population Size includ:
                                        % 1. Global Best Seller
                                        % 2. Global Worst Seller
                                        % 3. Yesterday Best Seller
                                        % 4. New Investor

CustomerPop = 12;                     % Customers Population Size

Opportunities = 4;                   % Number of opportunities
RewardPunishment = 1;                % Reward and punishment value

alfa = 0.9;                          % Inertia to target
alfaCoef=0.998;                      % Inertia update coefficint in each iteration

%% 1 the first selected sellers

% Empty Seller Structure
EmptySeller.Position = [];
EmptySeller.Objective = [];
EmptySeller.Attraction = [];
EmptySeller.NearestIndex=[];

Seller = repmat(EmptySeller,SellerPop,1);   %Initial Sellers population arraySeller(i).Objective

for i = 1:SellerPop
    %Initial Seller Position
    Seller(i).Position = unifrnd(VarMin,VarMax,VarSize);
    %Evaluation
    Seller(i).Objective = CostFunction(Seller(i).Position);
end

% Seller Nearest and Attraction
NearestIndex=NearestFunction(Seller,SellerPop);
AttractionIndex=Attractionfunction(Seller,FunctionType);
for j=1:SellerPop
    Seller(j).NearestIndex=NearestIndex(j);
    Seller(j).Attraction=AttractionIndex(j);
end

%%  Structures
% Empty Customer Structure
EmptyCustomer.Position = [];
EmptyCustomer.Objective = [];
EmptyCustomer.Characteristics.Gender = [];
EmptyCustomer.Characteristics.Age = [];
EmptyCustomer.Characteristics.FinancialSituation = [];
EmptyCustomer.Characteristics.FreeTime = [];
EmptyCustomer.Motivation = [];
EmptyCustomer.Inquiries = [];
EmptyCustomer.TargetList = [];
EmptyCustomer.TargetIndex = [];
EmptyCustomer.Destination = [];
EmptyCustomer.Distance = [];
EmptyCustomer.Sigma = [];
EmptyCustomer.NBS = [];
EmptyCustomer.NGS = [];
EmptyCustomer.Direction = [];
EmptyCustomer.MyOpportunities = [];


EmptyDirection.Position=[];
EmptyDirection.StepSize=[];
EmptyDirection.Objective=[];

EmptySol.Position=[];%nan(VarSize)';
EmptySol.Objective=[];%nan(1);

Sol = repmat(EmptySol,1,1);
BestSols = repmat(EmptySol,MaxIt,1);

GlobalSol=[];


%% Main loop
for it=1:MaxIt
    %Initial Customer population array
    clear Customer
    Customer = repmat(EmptyCustomer,CustomerPop,1);

    % 2 customer entrance to the mall
    for c=1:CustomerPop

        %Characteristics
        [Customer(c).Characteristics.Gender,...
         Customer(c).Characteristics.Age,...
         Customer(c).Characteristics.FinancialSituation,...
         Customer(c).Characteristics.FreeTime,...
         Customer(c).Motivation]=...
                                 CharacteristicsFunction(c,CustomerPop);

        %Number of price inquiries  
        Customer(c).Inquiries=...
                              InquiryFunction( ...
                                       Customer(c).Characteristics.Age,...
                                       Customer(c).Characteristics.Gender,...
                                       Customer(c).Characteristics.FreeTime,...
                                       Customer(c).Characteristics.FinancialSituation);

        %Customer Destination
        Customer(c).TargetList = 1:SellerPop;
        Customer(c).TargetIndex=...
                                TargetFunction(Customer(c), Seller, FunctionType);

        if  Customer(c).TargetIndex ~= 0
            Customer(c).Destination = Seller(Customer(c).TargetIndex).Position;
        else
            Customer(c).Destination = unifrnd(VarMin,VarMax,VarSize);
        end
            %Update Target List
            Customer(c).TargetList(Customer(c).TargetList==Customer(c).TargetIndex)=[];

        % customer entrance
        if Customer(c).Motivation~=0
            % Standard deviation for new position
            Customer(c).Sigma=...
                abs(Customer(c).Destination - ...
                Seller(Seller(Customer(c).TargetIndex).NearestIndex).Position)./1;
            % Customer Positon
            Customer(c).Position=...
                normrnd(Customer(c).Destination,Customer(c).Sigma);    
        elseif Customer(c).Motivation==0
            Customer(c).Position=...
                unifrnd(VarMin,VarMax,VarSize);
        end
        Customer(c).Position= max(min(Customer(c).Position,VarMax),VarMin);
        Customer(c).Objective=CostFunction(Customer(c).Position);

        % Customer Distance to Destination
        Customer(c).Distance= Customer(c).Destination - Customer(c).Position ;

        %Initial Customer Direction array
        Customer(c).Direction=repmat(EmptyDirection,Customer(c).Inquiries,1);

    % 3 Move
        for e=1:Customer(c).Inquiries %inquiring for price

                %Step Controler
                if e<= Opportunities
                    Customer(c).NGS=0;
                    Customer(c).NBS=0;
                %count Number of bad/good Steps
                elseif strcmp(FunctionType ,'MIN')
                    % in MIN function count Number of bad/good Steps
                    if  Customer(c).Direction(e-1).Objective <= Customer(c).Direction(e-2).Objective
                        Customer(c).NGS = Customer(c).NGS + 1;
                    else
                        Customer(c).NBS = Customer(c).NBS + 1;
                    end
                elseif strcmp(FunctionType ,'MAX')
                    % in Max function count Number of bad/good Steps
                    if  Customer(c).Direction(e-1).Objective >= Customer(c).Direction(e-2).Objective
                        Customer(c).NGS = Customer(c).NGS + 1;
                    else
                        Customer(c).NBS = Customer(c).NBS + 1;
                    end
                end

                StepControl=1+(Customer(c).NBS ./ Opportunities)-(Customer(c).NGS ./ (Customer(c).Inquiries+1));

                % Determining the Step Size
                Customer(c).Direction(e).StepSize =StepControl * (...
                            exp((Customer(c).Inquiries*(e-Customer(c).NBS-1))/(21.4-3.529*Customer(c).Inquiries))-...
                            exp((Customer(c).Inquiries*(e-Customer(c).NBS))/(21.4-3.529*Customer(c).Inquiries))); 

                %current position controller (just for step 1)
                if e == 1
                    CurrentPosition = Customer(c).Position; 
                else
                    CurrentPosition = Customer(c).Direction(e-1).Position;
                end

                % Position Attraction
                PA=0;
                if numel(Customer(c).TargetList)>0
                    for jj=Customer(c).TargetList
                        PA=PA+(Seller(jj).Position-CurrentPosition)*Seller(jj).Attraction;
                    end
                end
                % Step
                Step = Customer(c).Direction(e).StepSize * unifrnd(0,2,VarSize).* ...
                        (alfa*(Customer(c).Distance))+...
                        (1-alfa).*(PA);

                %new Position
                Customer(c).Direction(e).Position= CurrentPosition + Step;
                Customer(c).Direction(e).Position = min(max(Customer(c).Direction(e).Position,VarMin),VarMax);

                % Objective calcuation
                Customer(c).Direction(e).Objective...
                            =CostFunction(Customer(c).Direction(e).Position);

    % 4 Opportunities and Actions
                % opportunities
                if e <= Opportunities
                    Customer(c).MyOpportunities=Opportunities;
                elseif strcmp(FunctionType ,'MIN')
                    % in min function calculate My Opportunities
                    if  Customer(c).Direction(e).Objective > Customer(c).Direction(e-1).Objective
                        Customer(c).MyOpportunities = Customer(c).MyOpportunities - RewardPunishment;
                    else
                        Customer(c).MyOpportunities = Customer(c).MyOpportunities + RewardPunishment;
                    end
                elseif strcmp(FunctionType ,'MAX')
                    % in Max function count Number of bad/good Steps
                    if  Customer(c).Direction(e).Objective < Customer(c).Direction(e-1).Objective
                        Customer(c).MyOpportunities = Customer(c).MyOpportunities - RewardPunishment;
                    else
                        Customer(c).MyOpportunities = Customer(c).MyOpportunities + RewardPunishment;
                    end
                end

                %Action
                if Customer(c).MyOpportunities <= 0
                    
                    %NEW Customer Destination
                    if numel(Customer(c).TargetList) == 0 %any target can't get better
                        Customer(c).Motivation = 0;
                    end
                    Customer(c).TargetIndex=...
                                            TargetFunction(Customer(c), Seller, FunctionType);
                    if  Customer(c).TargetIndex ~= 0
                        Customer(c).Destination = Seller(Customer(c).TargetIndex).Position;
                    else
                        Customer(c).Destination = unifrnd(VarMin,VarMax,VarSize);
                    end
                    %Update Target List
                    Customer(c).TargetList(Customer(c).TargetList==Customer(c).TargetIndex)=[];

                    % NEW Distance
                    Customer(c).Distance= Customer(c).Destination - Customer(c).Position;

                    % NEW Opportunities
                    Customer(c).MyOpportunities = Opportunities;
                end
        end

        %add today sol to list
        Sol.Position=[Sol.Position;[Customer(c).Position]';[Customer(c).Direction(:).Position]'];
        Sol.Objective=[Sol.Objective;[Customer(c).Objective]';[Customer(c).Direction(:).Objective]'];
    end

        % Sort today sol
        TodaySol=[Sol.Objective,Sol.Position];
        Sol.Objective=[];
        Sol.Position=[];
        TodaySol=unique(TodaySol,"rows");%Unique and sort a->z
        if strcmp(FunctionType,'MAX')
            TodaySol=sort(TodaySol,'descend');
        end

    % 5 Determine selected sellers
        %New investor
        NewPosition=zeros(VarSize);
        theta=0;

        NewSol=[];
        nSeller=[];


    for k=1:SellerPop
        Beta=unifrnd(-0.1,0.01);

                % Seller Distance to near seller
        Seller(k).Distance=sqrt(sumsqr(Seller(k).Position -...
                            Seller(Seller(k).NearestIndex).Position)) ;
            ii=0;
        while theta <= 360*(3+floor(it*5/(MaxIt)))
             ii=ii+1;
            theta =randi([-5,25])+theta;    %comulative theta angel


            NewPosition(1:2:end)=Seller(k).Distance*exp(Beta*deg2rad(theta))*cos(deg2rad(theta));
            if nVar > 1
            NewPosition(2:2:end)=Seller(k).Distance*exp(Beta*deg2rad(theta))*sin(deg2rad(theta));
            end

            nSeller(ii).Position=Seller(k).Position+NewPosition;%#ok
            nSeller(ii).Position= max(min(nSeller(ii).Position,VarMax),VarMin);%#ok

            nSeller(ii).Objective = CostFunction(nSeller(ii).Position);%#ok
        end
            NewSol=[NewSol;[nSeller.Objective]',[nSeller.Position]'];      %#ok
            nSeller=[];
            theta = 0;
    end

        GlobalSol=[GlobalSol;TodaySol;NewSol];             %#ok
        GlobalSol=unique(GlobalSol,"rows");%Unique and sort a->z
        if strcmp(FunctionType,'MAX')
            GlobalSol=sort(GlobalSol,'descend');
        end
        
% MallCheck;

        %Selected Sellers
        Seller(1).Position = GlobalSol(1,2:end)';    %Global Best Seller

        if GlobalSol(1,1)~=TodaySol(1,1)
        Seller(2).Position = TodaySol(1,2:end)';     %Today Best Seller
        else
        Seller(2).Position = TodaySol(2,2:end)';     %Today Best Seller
        end
        TodaySol=[];                                 % Clear Today list
        %Zero and One Checker
            if it == floor(MaxIt*1/3) && VarMax >= 0 && VarMin <= 0
                Seller(2).Position=zeros(VarSize)+unifrnd(-1e-5,1e-5);
            elseif it == floor(MaxIt*3/4) && VarMax >= 0 && VarMin <= 0
                Seller(2).Position=ones(VarSize)+unifrnd(-1e-25,1e-25);
            end

        NewSol=unique(NewSol,"rows");%Unique and sort a->z
        if strcmp(FunctionType,'MAX')
            GlobalSol=sort(GlobalSol,'descend');
        end

        Seller(3).Position = NewSol(2,2:end)';       % New Investor
        NewSol=[];                                   % Clear New Sol list
        Seller(4).Position = GlobalSol(end,2:end)';  %Global Worst Seller

        for kk=1:SellerPop
            Seller(kk).Position=min(max(Seller(kk).Position,VarMin),VarMax);
        end


        for i = 1:SellerPop
            %Evaluation
            Seller(i).Objective = CostFunction(Seller(i).Position);
        end

        % update Seller Nearest and Attraction
        NearestIndex=NearestFunction(Seller,SellerPop);
        AttractionIndex=Attractionfunction(Seller,FunctionType);
        for j=1:SellerPop
            Seller(j).NearestIndex=NearestIndex(j);
            Seller(j).Attraction=AttractionIndex(j);
        end

        %Best Sol
        BestSols(it).Position = GlobalSol(1,2:end)';
        BestSols(it).Objective  = GlobalSol(1,1)';

        disp(['Iteration: ' num2str(it)  '  Best Objective = ' num2str(BestSols(it).Objective) ]);
    alfa = alfa * alfaCoef;


end


toc
%% results
        BestSol=BestSols(end);
        BestSol.Position=BestSol.Position';

        figure(1);
%         plot([BestSols.Objective],'LineWidth',2);
        hold on;
        semilogy([BestSols.Objective],'LineWidth',2);
        xlabel('Iteration');
        ylabel('Best Cost');

%% Functions
%Characteristics
function [Gender,Age,FinancialSituation,FreeTime,Motivation]=CharacteristicsFunction(i,CustomerPop)
    %Customer Gender 1--> Male , 0--> Female
    Gender = randi([0,1]);
    %Customer Age 1--> Old , 0--> Young
    Age = randi([0,1]);
    %Customer Financial Situation 2--> Poor , 1--> Medium , 0--> Rich
    FinancialSituation = randi([0,2]);
    %Customer Free Time 1--> Busy , 0--> Free
    FreeTime = randi([0,1]);
        if i <= floor(CustomerPop * 0.8)
        Motivation = 1;
        else
        Motivation = 0;
        end
end
%inquiry
function Inquiries = InquiryFunction(Age,Gender,FreeTime,FinancialSituation)
    %Customer Visit Duration
    VisitDuration = normrnd(170,30) ...
        - 26 * Age ...
        - 10 * FreeTime ...
        - 29 * FinancialSituation ...
        - 45 * Gender;
    %Customer Number of calls for price. at least 5 time inquiry
    Inquiries = max([7, ceil(VisitDuration./(5+randi([0,3])))]);
end
% Target with Roulett Wheel
function TargetIndex=TargetFunction(Customer, Seller , FunctionType)
    if Customer.Motivation == 0 || numel(Customer.TargetList) < 1
        TargetIndex = 0;
    else
        if strcmp(FunctionType ,'MIN')
            Cost=1./([Seller.Objective]+eps);
        elseif strcmp(FunctionType ,'MAX')
            Cost=Seller.Objective;
        end
        C=cumsum(Cost)./sum(Cost);
        TargetIndex = find(rand<=C,1,"first");
    end

end
% Near Seller
function NearestIndex=NearestFunction(Seller,SellerPop)
   route=nan(SellerPop);% should use nan instead of zerose or once
   for i=1:SellerPop-1
       for j=i+1:SellerPop
           route(i,j)=sqrt(sumsqr(Seller(i).Position-Seller(j).Position));
           route(j,i)=route(i,j);
       end
   end
[~,NearestIndex]=min(route,[],2,"omitnan");
end
%Sellers Attraction
function AttractionIndex=Attractionfunction(Seller,FunctionType)
    %Best / Worst Objective
    if strcmp(FunctionType ,'MIN')
        BestObj=min([Seller.Objective]);
        WorstObj=max([Seller.Objective]);
    elseif strcmp(FunctionType ,'MAX')
        BestObj=max([Seller.Objective]);
        WorstObj=min([Seller.Objective]);
    end
    AttractionIndex= abs([Seller.Objective] - WorstObj)./abs(WorstObj-BestObj);
end