clear all

Year=[2005,2010,2015];
n_p=200;
F_HH_3=struct();
  
load EXIOBASE_HHFD_split_by_age

EXIOBASE_HH_FD=load([pwd '\EX_v3_7_mat_current_price\IOT_'  '2005' '_pxp.mat']);
 name=EXIOBASE_HH_FD.('meta').('raw_countries')    ;
 name1=[{'BR'},{'WA'},{'WL'},{'WE'},{'WF'},{'WM'}];
 max=[1:49];
 CO2_total=struct();

 GHG=EXIOBASE_HH_FD.IO.char(215,:);
 CO2_F_P_1=zeros(1,6);
 
for i=1:size(Year,2)
        
    % Deal with Final demand (assuming import structure identical among income groups)
    EXIOBASE_HH_FD=load([pwd '\EX_v3_7_mat_current_price\IOT_'  num2str(Year(i)) '_pxp.mat']);
        
    Country=fieldnames(Household_Final_Demand_split.(['IO_',num2str(Year(i))]));
    
       % Direct emissions
        for idd=1:49*7
          Y_direct(:,idd)=transpose(GHG*EXIOBASE_HH_FD.IO.hhld_emis(:,:,idd));% CO2eq
          
        end
          Y_direct1=Y_direct;   
                
          
     %Conversion to constant price
     Constant_F=readtable('macro_desire_const_euro.xls', 'FileType', 'spreadsheet','Sheet',['Y',num2str(Year(i))]);      
     Constant_F_H=Constant_F.('HouseholdConsumptionExpenditure_includingNon_profitInstitutions');  
    
     %%% Calucating E L
    A=EXIOBASE_HH_FD.('IO').('A');
    I=eye(size(A,1));
    L=(I-A)^-1;
    
    %%%
    %input=(EXIOBASE_HH_FD.('IO').('x'));
    %input(find(input<1))=0;
     
    Invent=(EXIOBASE_HH_FD.('IO').('S'));%.*input_delfator1');
    %Intensity=(Invent)./input';
    Invent(isnan(Invent))=0;
      Invent(isinf(Invent))=0;
    %Intensity(isinf(Intensity))=0;
     Intensity= GHG* Invent;

    %%% 
    CO2multiplier=diag(Intensity)*L;
     
    
    for k=1:size(Country,1)%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Country
    
    FD_quintile=Household_Final_Demand_split.(['IO_',num2str(Year(i))]).(Country{k}).Household_FD_BasicPrice_pq;  
    
    k1=sum(strcmp(name,Country{k}).*[1:49]');
      
    % final demand
    F_HH=EXIOBASE_HH_FD.IO.Y(:,(k1-1)*7+1);
    Y_direct2=Y_direct1(:,(k1-1)*7+1);
    
    deflator=(Constant_F_H(k1,:)./1000000)./sum(F_HH,1);
      
    % constant price input
      
    share_exp_pq = FD_quintile ./ repmat(sum(FD_quintile, 2), 1, size(FD_quintile,2));
    share_exp_pq(isnan(share_exp_pq)) = 0;
    FD_per_quintile = repmat(share_exp_pq, 49, 1) .* repmat(F_HH, 1, size(FD_quintile,2));
    
    FD_HH_Quintile=FD_per_quintile.*repmat(deflator,9800,1);%
    
         
    % Direct emissions by quintile
    Share_direct=(FD_quintile./(sum(FD_quintile,2)));
    Share_direct(isnan(Share_direct))=0;
    Y_direct3(:,(k-1)*size(FD_quintile,2)+1:(k-1)*size(FD_quintile,2)+size(FD_quintile,2))=Y_direct2.*Share_direct;
    
    CO2_total.(['IO_',num2str(Year(i))]).Direct=Y_direct3;
      
  
 for c1=1:size(FD_HH_Quintile,2)
    CO2_F_M=CO2multiplier*diag(FD_HH_Quintile(:,c1));%(FD_HH_Quintile(:,cl))
    CO2_F(:,c1)=sum(CO2_F_M,2);
    CO2_F_c(c1,:)=sum(CO2_F_M,1);
    
    %%%%% Backward linkage   
    CO2_F_M1=zeros(9800,200);
    for i1=1:49
    CO2_F_M1=CO2_F_M(:,(i1-1)*200+1:(i1-1)*200+200)+CO2_F_M1;
    end
    CO2_F_M2(:,(c1-1)*200+1:(c1-1)*200+200)=CO2_F_M1;
 end
   
    for i1=1:49
    CO2_F_M3(i1,:)=sum(CO2_F_M2((i1-1)*200+1:(i1-1)*200+200,:),1);
    end

 %%%%%%%%%%%%    
    if k==29 %%% US
    CO2_F_M31(:,(1-1)*200+1:(1-1)*200+200)=CO2_F_M3(:,(1-1)*200+1:(1-1)*200+200)+CO2_F_M3(:,(2-1)*200+1:(2-1)*200+200)./2;
    CO2_F_M31(:,(2-1)*200+1:(2-1)*200+200)=CO2_F_M3(:,(2-1)*200+1:(2-1)*200+200)./2+CO2_F_M3(:,(3-1)*200+1:(3-1)*200+200);
    
    CO2_F_M31(:,(3-1)*200+1:(3-1)*200+200)=CO2_F_M3(:,(4-1)*200+1:(4-1)*200+200)+CO2_F_M3(:,(5-1)*200+1:(5-1)*200+200)./2;
    CO2_F_M31(:,(4-1)*200+1:(4-1)*200+200)=CO2_F_M3(:,(5-1)*200+1:(5-1)*200+200)./2+CO2_F_M3(:,(6-1)*200+1:(6-1)*200+200);
        
    CO2_F_M32=CO2_F_M31;
    elseif k==30 %% JP
    CO2_F_M31(:,(1-1)*200+1:(1-1)*200+200)=CO2_F_M3(:,(1-1)*200+1:(1-1)*200+200);
    CO2_F_M31(:,(2-1)*200+1:(2-1)*200+200)=CO2_F_M3(:,(2-1)*200+1:(2-1)*200+200)+CO2_F_M3(:,(3-1)*200+1:(3-1)*200+200)./2;
    CO2_F_M31(:,(3-1)*200+1:(3-1)*200+200)=(CO2_F_M3(:,(3-1)*200+1:(3-1)*200+200))./2+CO2_F_M3(:,(4-1)*200+1:(4-1)*200+200);
    CO2_F_M31(:,(4-1)*200+1:(4-1)*200+200)=CO2_F_M3(:,(5-1)*200+1:(5-1)*200+200)+CO2_F_M3(:,(6-1)*200+1:(6-1)*200+200);
            
    CO2_F_M32=CO2_F_M31;  
    
    elseif k==31 %% AUS
    CO2_F_M31(:,(1-1)*200+1:(1-1)*200+200)=CO2_F_M3(:,(1-1)*200+1:(1-1)*200+200)+CO2_F_M3(:,(2-1)*200+1:(2-1)*200+200)./2;
    CO2_F_M31(:,(2-1)*200+1:(2-1)*200+200)=CO2_F_M3(:,(2-1)*200+1:(2-1)*200+200)./2+CO2_F_M3(:,(3-1)*200+1:(3-1)*200+200);
    CO2_F_M31(:,(3-1)*200+1:(3-1)*200+200)=CO2_F_M3(:,(4-1)*200+1:(4-1)*200+200)+CO2_F_M3(:,(5-1)*200+1:(5-1)*200+200)./2;
    CO2_F_M31(:,(4-1)*200+1:(4-1)*200+200)=CO2_F_M3(:,(5-1)*200+1:(5-1)*200+200)./2+CO2_F_M3(:,(6-1)*200+1:(6-1)*200+200);
        
    CO2_F_M32=CO2_F_M31;
        
    else
    CO2_F_M32=CO2_F_M3;
    end
    
    CO2_Trade(:,(k-1)*size(CO2_F,2)+1:(k-1)*size(CO2_F,2)+size(CO2_F,2))=CO2_F;%%% Forward linkage
    CO2_Trade_back(:,(k-1)*4*200+1:(k-1)*4*200+4*200)=CO2_F_M32;%%% Backward linkage
        
    CO2_F_P_1(k,1:size(FD_HH_Quintile,2))=sum(CO2_F_c,2)'+sum(Y_direct2.*Share_direct,1);
       
    for j=1:size(FD_HH_Quintile,2)
    for i1=1:49
    CO2_F_c1(i1,:)=CO2_F_c(j,(i1-1)*200+1:(i1-1)*200+200);
    end
    CO2_F_c2(j,:)=sum(CO2_F_c1,1);
    end    

    CO2_F_P(:,(k-1)*size(CO2_F_c2,1)+1:(k-1)*size(CO2_F_c2,1)+size(CO2_F_c2,1))=CO2_F_c2';
    
    
    clearvars CO2_F_c CO2_F_c2 CO2_F CO2_F_M3 CO2_F_M2
    end
  
    CO2_total.(['IO_',num2str(Year(i))]).Total=CO2_F_P+Y_direct3;%%% Direct+Indirect
    CO2_total.(['IO_',num2str(Year(i))]).by_products=CO2_F_P;;%%% Indirect
    CO2_total.(['IO_',num2str(Year(i))]).national=CO2_F_P_1;%%% Direct+Indirect
    CO2_total.(['IO_',num2str(Year(i))]).Trade=CO2_Trade;
    CO2_total.(['IO_',num2str(Year(i))]).Trade_back=CO2_Trade_back;
    clearvars FD_HH_Quintile 
end

clearvars -except CO2_total Year Country Household_Final_Demand_split name

%%%% JP  US AUS
%%% 
for i=1:size(Year,2)
for k=1:size(Country,1)
    CO2_total_2=CO2_total.(['IO_',num2str(Year(i))]).national(k,:);
    
    xxx=strcmp(Country{k},'US')| strcmp(Country{k},'JP')| strcmp(Country{k},'AU');
    
    if xxx==0;
    CO2_total_new.(['IO_',num2str(Year(i))]).national(k,:)=CO2_total_2;
    elseif ([Country{k}])=='US'
    CO2_total_1=CO2_total.(['IO_',num2str(Year(i))]).national(k,:);
    CO2_total_2(:,1)=CO2_total_1(:,1)+CO2_total_1(:,2)./2;
    CO2_total_2(:,2)=CO2_total_1(:,2)./2+CO2_total_1(:,3);
    CO2_total_2(:,3)=CO2_total_1(:,4)+CO2_total_1(:,5)./2;
    CO2_total_2(:,4)=(CO2_total_1(:,5))./2+CO2_total_1(:,6);
    CO2_total_new.(['IO_',num2str(Year(i))]).national(k,:)=CO2_total_2;
       
    elseif ([Country{k}])=='JP'
    CO2_total_1=CO2_total.(['IO_',num2str(Year(i))]).national(k,:);    
    CO2_total_2(:,1)=CO2_total_1(:,1);
    CO2_total_2(:,2)=CO2_total_1(:,2)+CO2_total_1(:,3)./2;
    CO2_total_2(:,3)=CO2_total_1(:,3)./2+CO2_total_1(:,4);
    CO2_total_2(:,4)=CO2_total_1(:,5)+CO2_total_1(:,6);
    CO2_total_new.(['IO_',num2str(Year(i))]).national(k,:)=CO2_total_2;
    
    elseif ([Country{k}])=='AU'
    CO2_total_1=CO2_total.(['IO_',num2str(Year(i))]).national(k,:);
    CO2_total_2(:,1)=CO2_total_1(:,1)+CO2_total_1(:,2)./2;
    CO2_total_2(:,2)=CO2_total_1(:,2)./2+CO2_total_1(:,3);
    CO2_total_2(:,3)=CO2_total_1(:,4)+CO2_total_1(:,5)./2;
    CO2_total_2(:,4)=(CO2_total_1(:,5))./2+CO2_total_1(:,6);
    CO2_total_new.(['IO_',num2str(Year(i))]).national(k,:)=CO2_total_2;
       
    end
    
end
CO2_total_new.(['IO_',num2str(Year(i))]).national(:,5:6)=[];
end

%%%%% Product
for  i=1:size(Year,2)
for k=1:size(Country,1)
      xxx=strcmp(Country{k},'US')| strcmp(Country{k},'JP')| strcmp(Country{k},'AU');
    if xxx==0;
    CO2_total_2P=CO2_total.(['IO_',num2str(Year(i))]).Total(:,(k-1)*4+1:(k-1)*4+4);
    CO2_total_new.(['IO_',num2str(Year(i))]).Total(:,(k-1)*4+1:(k-1)*4+4)=CO2_total_2P;
    
    elseif ([Country{k}])=='US'
     Bridge=xlsread('Bridge', ([Country{k}]));
     Bridge(isnan(Bridge))=0;
    CO2_total_2P=CO2_total.(['IO_',num2str(Year(i))]).Total(:,(k-1)*6+1:(k-1)*6+6);    
    CO2_total_21=CO2_total_2P*Bridge;
   
    CO2_total_new.(['IO_',num2str(Year(i))]).Total(:,(k-1)*4+1:(k-1)*4+4)=CO2_total_21;
       
    elseif ([Country{k}])=='JP'
      Bridge=xlsread('Bridge', ([Country{k}]));
     Bridge(isnan(Bridge))=0;
    CO2_total_2P=CO2_total.(['IO_',num2str(Year(i))]).Total(:,(k-1)*6+1:(k-1)*6+6);    
    CO2_total_21=CO2_total_2P*Bridge;
    CO2_total_new.(['IO_',num2str(Year(i))]).Total(:,(k-1)*4+1:(k-1)*4+4)=CO2_total_21;
    
    elseif ([Country{k}])=='AU'
     Bridge=xlsread('Bridge', ([Country{k}]));
     Bridge(isnan(Bridge))=0;
    CO2_total_2P=CO2_total.(['IO_',num2str(Year(i))]).Total(:,(k-1)*6+1:(k-1)*6+6);    
    CO2_total_21=CO2_total_2P*Bridge;
   
    CO2_total_new.(['IO_',num2str(Year(i))]).Total(:,(k-1)*4+1:(k-1)*4+4)=CO2_total_21;
    
    
    end
end
end



%%%%%%%%-------------------------Trade_product forward-----------
[~,~,Bridge_Product]=xlsread('Bridge_product.xlsx');
Name_prod=Bridge_Product(1,2:10)';
Bridge_Product=cell2mat(Bridge_Product(2:201,2:10));
Bridge_Product(isnan(Bridge_Product))=0;
%%%% Trade_product forward
for i=1:size(Year,2)
for k=1:size(Country,1)
    xxx=([Country{k}])~='US'&([Country{k}])~='JP';
    if xxx(:,1)==1;
    CO2_total_trade=CO2_total.(['IO_',num2str(Year(i))]).Trade(:,(k-1)*4+1:(k-1)*4+4);
    
    for j=1:49
    CO2_total_trade1((j-1)*size(Bridge_Product,2)+1:(j-1)*size(Bridge_Product,2)+size(Bridge_Product,2),(k-1)*4+1:(k-1)*4+4)=Bridge_Product'*CO2_total_trade((j-1)*200+1:(j-1)*200+200,:)./1e+6;% Million
    end
    
    elseif ([Country{k}])=='US'
     Bridge=xlsread('Bridge', ([Country{k}]));
     Bridge(isnan(Bridge))=0;   
      
     CO2_total_trade=CO2_total.(['IO_',num2str(Year(i))]).Trade(:,(k-1)*6+1:(k-1)*6+6)*Bridge;
    
    for j=1:49
    CO2_total_trade1((j-1)*size(Bridge_Product,2)+1:(j-1)*size(Bridge_Product,2)+size(Bridge_Product,2),(k-1)*4+1:(k-1)*4+4)=Bridge_Product'*CO2_total_trade((j-1)*200+1:(j-1)*200+200,:)./1e+6;% Million
    end 
          
    elseif ([Country{k}])=='JP'  
     Bridge=xlsread('Bridge', ([Country{k}]));
     Bridge(isnan(Bridge))=0;
         
     CO2_total_trade=CO2_total.(['IO_',num2str(Year(i))]).Trade(:,(k-1)*6+1:(k-1)*6+6)*Bridge;
    
    for j=1:49
    CO2_total_trade1((j-1)*size(Bridge_Product,2)+1:(j-1)*size(Bridge_Product,2)+size(Bridge_Product,2),(k-1)*4+1:(k-1)*4+4)=Bridge_Product'*CO2_total_trade((j-1)*200+1:(j-1)*200+200,:)./1e+6;% Million
    end 
     
    end
    CO2_total_new.(['IO_',num2str(Year(i))]).Trade=CO2_total_trade1;
end
end


%%%% SORT
Age_n=[1,2,3,4];
for i=1:size(Year,2)
CO2_total_trade1=CO2_total_new.(['IO_',num2str(Year(i))]).Trade;

CO2_total_trade_Domestic=zeros(10*49,32*4);
%%% Extract Domestic
for j=1:size(Country,1)

    j1=sum(strcmp(name,Country{j}).*[1:49]');
CO2_total_trade_Domestic((j1-1)*9+1:(j1-1)*9+9,(j-1)*4+1:(j-1)*4+4)=CO2_total_trade1((j1-1)*9+1:(j1-1)*9+9,(j-1)*4+1:(j-1)*4+4);
end


%%%Year
Year1=repmat(Year(i),size(CO2_total_trade1(:),1),1);

%%%Consumer name
for j=1:size(Country,1)
Consumer((j-1)*9*49*4+1:(j-1)*4*9*49+4*9*49,1)=repmat(Country(j),4*9*49,1);
end

%%%Age
for k=1:4
Age_n1((k-1)*49*9+1:(k-1)*49*9+49*9,1)=repmat(Age_n(k),49*9,1);
end

%%% Producer name
for j1=1:size(name,1)
Producer((j1-1)*9+1:(j1-1)*9+9,1)=repmat(name(j1),9,1);
end

%%%% Product

Product=repmat(Name_prod,49,1);


[~,~,Typc]=xlsread('Bridge_product.xlsx','TypeC');
%%%% Type C
for  j=1:size(Typc,1)
Typc1((j-1)*9*49*4+1:(j-1)*9*49*4+4*9*49,1)=repmat(Typc(j,2),4*9*49,1);
end

[~,~,Mature]=xlsread('Bridge_product.xlsx','MatureP');
%%%% Mature
for  j=1:size(Mature,1)
Mature1((j-1)*9+1:(j-1)*9+9,1)=repmat(Mature(j,2),9,1);
end


Producer1= repmat(repmat(Producer,4,1),32,1);
Product1= repmat(repmat(Product,4,1),32,1);
Mature2=repmat(repmat(Mature1,4,1),32,1);

CO2_FLOW=table(Year1,Consumer,repmat(Age_n1,32,1),Producer1,Product1,Typc1,Mature2,CO2_total_trade1(:),CO2_total_trade_Domestic(:));
CO2_FLOW1((i-1)*size(CO2_FLOW,1)+1:(i-1)*size(CO2_FLOW,1)+size(CO2_FLOW,1),:)=CO2_FLOW;
end

%---------------------------------------Backward leakage--------------
% Trade_product bakcward
[~,~,Bridge_Product]=xlsread('Bridge_product.xlsx');
Name_prod=Bridge_Product(1,2:11)';
Bridge_Product=cell2mat(Bridge_Product(2:201,2:end));
Bridge_Product(isnan(Bridge_Product))=0;

for i=1:size(Year,2)
for k=1:size(Country,1)
 
    CO2_total_trade=CO2_total.(['IO_',num2str(Year(i))]).Trade_back(:,(k-1)*4*200+1:(k-1)*4*200+4*200);
    
    for j=1:4    
        CO2_total_trade1(:,(j-1)*size(Bridge_Product,2)+1:(j-1)*size(Bridge_Product,2)+size(Bridge_Product,2))=CO2_total_trade(:,(j-1)*200+1:(j-1)*200+200)*Bridge_Product./1e+6;% Million
    end
    
    CO2_total_tradeb(:,(k-1)*size(CO2_total_trade1,2)+1:(k-1)*size(CO2_total_trade1,2)+size(CO2_total_trade1,2))=CO2_total_trade1;
  
end
    CO2_total_new.(['IO_',num2str(Year(i))]).Trade_back=CO2_total_tradeb;
end


%%%% SORT
Age_n=[1,2,3,4];
for i=1:size(Year,2)
CO2_total_trade1=CO2_total_new.(['IO_',num2str(Year(i))]).Trade_back;

CO2_total_trade_Domestic=zeros(49,32*4*10);
%%% Extract Domestic
for j=1:size(Country,1)
    j1=sum(strcmp(name,Country{j}).*[1:49]');
CO2_total_trade_Domestic(j1,(j-1)*4*10+1:(j-1)*4*10+4*10)=CO2_total_trade1(j1,(j-1)*4*10+1:(j-1)*4*10+4*10);
end
sum(CO2_total_trade_Domestic,1);

%%%Year
Year1=repmat(Year(i),size(CO2_total_trade1(:),1),1);

%%%Age
for k=1:4
Age_n1((k-1)*49*10+1:(k-1)*49*10+49*10,1)=repmat(Age_n(k),49*10,1);
end

%%%Consumer name
for j=1:size(Country,1)
Consumer((j-1)*10*49*4+1:(j-1)*4*10*49+4*10*49,1)=repmat(Country(j),4*10*49,1);
end


%%% Producer name
Producer=repmat(name,10*4*32,1);

%%%% Product
for i1=1:size(Name_prod)
Product((i1-1)*49+1:(i1-1)*49+49,1)=repmat(Name_prod(i1),49,1);
end
Product1= repmat(repmat(Product,4,1),32,1);

%%%% Type C
[~,~,Typc]=xlsread('Bridge_product.xlsx','TypeC');
for  j=1:size(Typc,1)
Typc1((j-1)*10*49*4+1:(j-1)*10*49*4+4*10*49,1)=repmat(Typc(j,2),4*10*49,1);
end
%%%% Mature
[~,~,Mature]=xlsread('Bridge_product.xlsx','MatureP');
Mature1=repmat(Mature(:,2),10*4*32,1);

CO2_FLOW=table(Year1,repmat(Age_n1,32,1),Consumer,Product1,Producer,Typc1,Mature1,CO2_total_trade1(:),CO2_total_trade_Domestic(:));
CO2_FLOW1((i-1)*size(CO2_FLOW,1)+1:(i-1)*size(CO2_FLOW,1)+size(CO2_FLOW,1),:)=CO2_FLOW;
end

%%%%%% Layout backward leakage
for  j=1:size(Typc,1)
Typc11((j-1)*10*4+1:(j-1)*10*4+4*10,1)=repmat(Typc(j,2),4*10,1);
end

%%%Age
for k=1:4
Age_n11((k-1)*10+1:(k-1)*10+10,1)=repmat(Age_n(k),10,1);
end
repmat(Age_n11,32,1);

for i=1:size(Country,1)
   Country1((i-1)*40+1:(i-1)*40+40,:)=repmat(Country(i),40,1);
end




%---------------------------------Aggregated Trade Matrix---------------------------
%%%Aggregated 
for i=1:size(Year,2)
   CO2_total_trade1=CO2_total_new.(['IO_',num2str(Year(i))]).Trade;
for k=1:49
   CO2_total_trade2(k,:)=sum(CO2_total_trade1((k-1)*9+1:(k-1)*9+9,:),1);
end
for k=1:32
    k1=sum(strcmp(name,Country{k}).*[1:49]');
   CO2_total_dome(1,(k-1)*4+1:(k-1)*4+4)=CO2_total_trade2(k1,(k-1)*4+1:(k-1)*4+4);
end   
   CO2_total_new.(['IO_',num2str(Year(i))]).Trade_c2c=CO2_total_trade2;
   Domestic{i}=CO2_total_dome;
end

for i=1:size(Country,1)
   Country1(:,(i-1)*4+1:(i-1)*4+4)=repmat(Country(i),1,4);
end


%
%------------------------------------LMDI-----------------------------------------%
LMDI=struct();
for i=1:size(Year,2)
   CO2_1=CO2_total_new.(['IO_',num2str(Year(i))]).national;
   CO2_p=CO2_total_new.(['IO_',num2str(Year(i))]).Total;
for k=1:size(Country,1)
    
   Info=(Household_Final_Demand_split.(['IO_',num2str(Year(i))]).([Country{k}]));
   xxx=strcmp(Country{k},'US')| strcmp(Country{k},'JP')| strcmp(Country{k},'AU');
   
   if xxx==0;
   CO2_2= CO2_1(k,:);
   
   HH=Info.Households_q;
   HH_pp=Info.Population_q;
   HH_size=HH_pp./HH;
   total_EX=sum(Info.Household_FD_BasicPrice_pq,1);
   Per_C_EX=total_EX./HH_pp;% 1000 Euro£»
   Structure_EX=Info.Household_FD_BasicPrice_pq./total_EX;
   co2_multiplier=CO2_p(:,(k-1)*4+1:(k-1)*4+4)./Info.Household_FD_BasicPrice_pq;
    Expenditure1=Info.Household_FD_BasicPrice_pq;
    
   elseif ([Country{k}])=='US'
   CO2_2= CO2_1(k,:);    
   Bridge=xlsread('Bridge', ([Country{k}]));
   Bridge(isnan(Bridge))=0;
   
   HH=Info.Consumer_Units_q*Bridge;
   HH_pp=Info.Population_q*Bridge;
   HH_size=HH_pp./HH;
   total_EX=sum(Info.Household_FD_BasicPrice_pq*Bridge,1);
   Per_C_EX=total_EX./HH_pp;% 1000 Euro£»
   Structure_EX=Info.Household_FD_BasicPrice_pq*Bridge./total_EX;
   co2_multiplier=(CO2_p(:,(k-1)*4+1:(k-1)*4+4))./(Info.Household_FD_BasicPrice_pq*Bridge);    
    Expenditure1=Info.Household_FD_BasicPrice_pq*Bridge;; 
    
     elseif ([Country{k}])=='AU'
   CO2_2= CO2_1(k,:);    
   Bridge=xlsread('Bridge', ([Country{k}]));
   Bridge(isnan(Bridge))=0;
   
   HH=Info.Households_q*Bridge;
   HH_pp=Info.Population_q*Bridge;
   HH_size=HH_pp./HH;
   total_EX=sum(Info.Household_FD_BasicPrice_pq*Bridge,1);
   Per_C_EX=total_EX./HH_pp;% 1000 Euro£»
   Structure_EX=Info.Household_FD_BasicPrice_pq*Bridge./total_EX;
   co2_multiplier=(CO2_p(:,(k-1)*4+1:(k-1)*4+4))./(Info.Household_FD_BasicPrice_pq*Bridge);    
   Expenditure1=Info.Household_FD_BasicPrice_pq*Bridge;; 
      
    
    
   elseif ([Country{k}])=='JP'
   CO2_2= CO2_1(k,:);    
   Bridge=xlsread('Bridge', ([Country{k}]));
   Bridge(isnan(Bridge))=0;
   
   HH=Info.Households_q*Bridge;
   HH_pp=Info.Population_q*Bridge;
   HH_size=HH_pp./HH;
   total_EX=sum(Info.Household_FD_BasicPrice_pq*Bridge,1);
   Per_C_EX=total_EX./HH_pp;% 1000 Euro£»
   Structure_EX=Info.Household_FD_BasicPrice_pq*Bridge./total_EX;
   co2_multiplier=(CO2_p(:,(k-1)*4+1:(k-1)*4+4))./(Info.Household_FD_BasicPrice_pq*Bridge);
    Expenditure1=Info.Household_FD_BasicPrice_pq*Bridge;;  
      end
   
   LMDI.(['IO_',num2str(Year(i))]).([Country{k}]).CO2=CO2_2;
   LMDI.(['IO_',num2str(Year(i))]).([Country{k}]).CO2_p=CO2_p(:,(k-1)*4+1:(k-1)*4+4); 
   LMDI.(['IO_',num2str(Year(i))]).([Country{k}]).HH=HH;
   LMDI.(['IO_',num2str(Year(i))]).([Country{k}]).Pop=HH_pp;
   LMDI.(['IO_',num2str(Year(i))]).([Country{k}]).HZ=HH_size;
   LMDI.(['IO_',num2str(Year(i))]).([Country{k}]).TE=total_EX;
   LMDI.(['IO_',num2str(Year(i))]).([Country{k}]).PE=Per_C_EX;
   LMDI.(['IO_',num2str(Year(i))]).([Country{k}]).SE=Structure_EX;
   LMDI.(['IO_',num2str(Year(i))]).([Country{k}]).int=co2_multiplier;
   LMDI.(['IO_',num2str(Year(i))]).([Country{k}]).Ex=Expenditure1;
   end
end

%%%%%%%% LMDI-I additive model
for i=1:size(Year,2)-1
    for k=1:size(Country,1)
   CO2_0=LMDI.(['IO_',num2str(Year(i))]).([Country{k}]).CO2;
   CO2_p_0= LMDI.(['IO_',num2str(Year(i))]).([Country{k}]).CO2_p;
   HH_0=LMDI.(['IO_',num2str(Year(i))]).([Country{k}]).HH;
   Pop_0=LMDI.(['IO_',num2str(Year(i))]).([Country{k}]).Pop;
   HZ_0=LMDI.(['IO_',num2str(Year(i))]).([Country{k}]).HZ;
   TE_0=LMDI.(['IO_',num2str(Year(i))]).([Country{k}]).TE;
   PE_0=LMDI.(['IO_',num2str(Year(i))]).([Country{k}]).PE;
   SE_0=LMDI.(['IO_',num2str(Year(i))]).([Country{k}]).SE;
   int_0=LMDI.(['IO_',num2str(Year(i))]).([Country{k}]).int;
   
   
   CO2_t=LMDI.(['IO_',num2str(Year(i+1))]).([Country{k}]).CO2;
   CO2_p_t= LMDI.(['IO_',num2str(Year(i+1))]).([Country{k}]).CO2_p;
   HH_t=LMDI.(['IO_',num2str(Year(i+1))]).([Country{k}]).HH;
   Pop_t=LMDI.(['IO_',num2str(Year(i+1))]).([Country{k}]).Pop;
   HZ_t=LMDI.(['IO_',num2str(Year(i+1))]).([Country{k}]).HZ;
   TE_t=LMDI.(['IO_',num2str(Year(i+1))]).([Country{k}]).TE;
   PE_t=LMDI.(['IO_',num2str(Year(i+1))]).([Country{k}]).PE;
   SE_t=LMDI.(['IO_',num2str(Year(i+1))]).([Country{k}]).SE;
   int_t=LMDI.(['IO_',num2str(Year(i+1))]).([Country{k}]).int;

%CO2_p_t(find(CO2_p_t==0))=1e-500;
%CO2_p_0(find(CO2_p_0==0))=1e-500;


PE_t(isnan(PE_t))=1e-300;
SE_t(isnan(SE_t))=1e-300;
int_t(isnan(int_t))=1e-300;
%SE_t(find(SE_t==0))=1e-500;
%int_t(find(int_t==0))=1e-500;

Dtot1=sum(sum(CO2_p_t))-sum(sum(CO2_p_0));

%%% total Household number
for q=1:4
for pi=1:200    
Dact(pi,q)=(CO2_p_t(pi,q)-CO2_p_0(pi,q))./(log(CO2_p_t(pi,q))-log(CO2_p_0(pi,q))).*log(sum(HH_t,2)./sum(HH_0,2));
end
end
Dact(isnan(Dact))=0;
Dact(isinf(Dact))=0;
Factor_2005(1,:)=((sum(Dact,1)));

%%% household structure
for q=1:4
S=[HH_0(q)./sum(HH_0,2),HH_t(q)./sum(HH_t,2)];%%
for pi=1:200    
L=(CO2_p_t(pi,q)-CO2_p_0(pi,q))./(log(CO2_p_t(pi,q))-log(CO2_p_0(pi,q))); 
Dstr(pi,q)=L.*log(S(:,2)./S(:,1));
end
end
Dstr(isnan(Dstr))=0;
Dstr(isinf(Dstr))=0;
Factor_2005(2,:)=((sum(Dstr,1)));

%%% Household size
for q=1:4
for pi=1:200    
L=(CO2_p_t(pi,q)-CO2_p_0(pi,q))./(log(CO2_p_t(pi,q))-log(CO2_p_0(pi,q))); 
DHZ(pi,q)=L.*log(HZ_t(q)./HZ_0(q));
end
end
DHZ(isnan(DHZ))=0;
DHZ(isinf(DHZ))=0;
Factor_2005(3,:)=((sum(DHZ,1)));

%%% per capita expenditure
for q=1:4
for pi=1:200    
L=(CO2_p_t(pi,q)-CO2_p_0(pi,q))./(log(CO2_p_t(pi,q))-log(CO2_p_0(pi,q))); 
DPEX(pi,q)=L.*log(PE_t(q)./PE_0(q));
end
end
DPEX(isnan(DPEX))=0;
DPEX(isinf(DPEX))=0;
Factor_2005(4,:)=((sum(DPEX,1)));

%%% consumption strucuture
for q=1:4
for pi=1:200    
L=(CO2_p_t(pi,q)-CO2_p_0(pi,q))./(log(CO2_p_t(pi,q))-log(CO2_p_0(pi,q))); 
DEXSTR(pi,q)=L.*log(SE_t(pi,q)./SE_0(pi,q));
end
end
DEXSTR(isnan(DEXSTR))=0;
DEXSTR(isinf(DEXSTR))=0;
Factor_2005(5,:)=((sum(DEXSTR,1)));

%%%intensity
for q=1:4
for pi=1:200    
L=(CO2_p_t(pi,q)-CO2_p_0(pi,q))./(log(CO2_p_t(pi,q))-log(CO2_p_0(pi,q))); 
Dint(pi,q)=L.*log(int_t(pi,q)./int_0(pi,q));
end
end
Dint(isnan(Dint))=0;
Dint(isinf(Dint))=0;
Factor_2005(6,:)=((sum(Dint,1)));

Delta=Dtot1-sum(sum(Factor_2005,2));

LMDI_C((k-1)*6+1:(k-1)*6+6,(i-1)*4+1:(i-1)*4+4)=Factor_2005./1000000000;

Check(k,i)=Delta./sum(sum(CO2_p_0));
end  
end


%%% Combaing Household number and strucuture
for i=1:32
LMDI_D((i-1)*5+1,:)=sum(LMDI_C((i-1)*6+1:(i-1)*6+2,:),1);
LMDI_D((i-1)*5+2:(i-1)*5+5,:)=LMDI_C((i-1)*6+3:(i-1)*6+6,:);
end


%%%% AGE 1+2

LMDI_D1(:,1)=LMDI_D(:,1)+LMDI_D(:,2);
LMDI_D1(:,2:3)=LMDI_D(:,3:4);
LMDI_D1(:,4)=LMDI_D(:,5)+LMDI_D(:,6);
LMDI_D1(:,5:6)=LMDI_D(:,7:8);

LMDI_D1(:);

emppp=[];

for i=1:32*4*3
    emppp1((i-1)*10+1:(i-1)*10+10,1)=repmat(emppp(i),10,1); 
end

  
%%%%% CO2 by products
for i=1:size(Year,2)
      for k=1:32
CO2_p_t= LMDI.(['IO_',num2str(Year(i))]).([Country{k}]).CO2_p;

LMDI_T(k,i)=sum(sum(CO2_p_t,1))./1000000000;
      end
end
LMDI_T(:);

% Household size
for i=1:size(Year,2)
    for k=1:32
      HZ_Y(k,:) =LMDI.(['IO_',num2str(Year(i))]).([Country{k}]).HZ;
      HZ_1.(['IO_',num2str(Year(i))])=HZ_Y;
    end
end

% Population
for i=1:size(Year,2)
    for k=1:32
      Pop1(k,:) =LMDI.(['IO_',num2str(Year(i))]).([Country{k}]).Pop;
      Pop12.(['IO_',num2str(Year(i))])=Pop1;
    end
end


%%%%%%%%%%%%% Extract the expenditure strucuture

for k=1:32
 for i=1:3   
  SE_P=LMDI.(['IO_',num2str(Year(i))]).([Country{k}]).SE;
%for q=1:4
  SE_P1(:,i)=SE_P(:,4);
%end
end
  Age_EXst.([Country{k}])=SE_P1;
end



%%%%% Layout
[~,~,LMDI_n]=xlsread('Bridge_product.xlsx','LMDI');

repmat(repmat(repmat(LMDI_n,32,1),4,1),2,1);

Age=[1;2;3;4];
    for j=1:size(Age,1)
       Age1(:,j)=repmat(repmat(Age(j),32,1),6,1);
    end
Age2=repmat(Age1,1,2);
Age2(:);


for i=1:size(Country,1)
Country_1((i-1)*6+1:(i-1)*6+6,:)=repmat(Country(i),6,1);
end
Country_2=repmat(Country_1,1,8);
Country_2(:);

%%%% Type
[~,~,Name_type]=xlsread('Bridge_product.xlsx','Typec');;%%% Copy from EXCEL
for i=1:32
   Name_type2((i-1)*6+1:(i-1)*6+6,1)=repmat(Name_type(i,2),6,1); 
end
Name_type3=repmat(Name_type2,1,8);
Name_type3(:);




%%%%%%%%%%%%%%%-------------------Per capita GHG ----------------------------
[~,~,Bridge_Product]=xlsread('Bridge_product.xlsx');
Name_prod=Bridge_Product(1,2:11)';
Bridge_Product=cell2mat(Bridge_Product(2:201,2:11));
Bridge_Product(isnan(Bridge_Product))=0;
%%%%% Per capita GHG
for i=1:size(Year,2)
    for k=1:size(Country,1)
   CO2_pc=LMDI.(['IO_',num2str(Year(i))]).([Country{k}]).CO2;
   CO2_p_pc= LMDI.(['IO_',num2str(Year(i))]).([Country{k}]).CO2_p;
   HH_pc=LMDI.(['IO_',num2str(Year(i))]).([Country{k}]).HH;
   Pop_pc=LMDI.(['IO_',num2str(Year(i))]).([Country{k}]).Pop;
   %HZ_pc=LMDI.(['IO_',num2str(Year(i))]).([Country{k}]).HZ;
   %TE_pc=LMDI.(['IO_',num2str(Year(i))]).([Country{k}]).TE;
   PE_pc=LMDI.(['IO_',num2str(Year(i))]).([Country{k}]).PE;
   %SE_pc=LMDI.(['IO_',num2str(Year(i))]).([Country{k}]).SE;
   %int_pc=LMDI.(['IO_',num2str(Year(i))]).([Country{k}]).int;
   Expenditure2=LMDI.(['IO_',num2str(Year(i))]).([Country{k}]).Ex;
   
   CO2_PC0(k,(i-1)*4+1:(i-1)*4+4)=CO2_pc./1e+9;
   CO2_PC1(k,(i-1)*4+1:(i-1)*4+4)=CO2_pc./Pop_pc./1e+3; %%% ton/person
   CO2_PC2(k,(i-1)*4+1:(i-1)*4+4)=PE_pc*1e+6;%%% Million Euro/person
   CO2_PC3(k,(i-1)*4+1:(i-1)*4+4)=Pop_pc;
   CO2_PC4((k-1)*10+1:(k-1)*10+10,(i-1)*4+1:(i-1)*4+4)=(CO2_p_pc'*Bridge_Product)'./1e+9;%%% Emissions by products
   CO2_PC5((k-1)*10+1:(k-1)*10+10,(i-1)*4+1:(i-1)*4+4)=(Expenditure2'*Bridge_Product)';%%% Expenditure by products
   CO2_PC6(k,(i-1)*4+1:(i-1)*4+4)=HH_pc;
   
   %CO2_PC7(k,(i-1)*4+1:(i-1)*4+4)=sum(Expenditure2,1);
    
    end
end

CO2_PC0(:); % GHG total
CO2_PC1(:); % Per capita GHG total
CO2_PC2(:); % Expenditure per capita total Euro
CO2_PC3(:); % Pop total
CO2_PC4(:); % GHG by product 
CO2_PC5(:); % Expenditure by product
CO2_PC6(:); %household number


for  j=1:32 
    CO2_PC3_1((j-1)*10+1:(j-1)*10+10,:)=repmat(CO2_PC3(j,:),10,1);
end

CO2_PC3_2=CO2_PC5./CO2_PC3_1;%%%%%  CO2_PC3_2   Expenditure/pop by sectors
CO2_PC_in=CO2_PC4./CO2_PC5; %%%%  CO2_PC_in   GHG multipliers


save 'AAA' CO2_PC_in CO2_PC5 CO2_PC4 CO2_PC3_2;
%%%%% Adult equivalent
Adult=xlsread('Result.xlsx','adult-equevalent');
Adult1(:,1:4)=Adult(:,1:4);
Adult1(:,5:8)=Adult(:,6:9);
Adult1(:,9:12)=Adult(:,11:14);

Total_people=Adult1.*CO2_PC6;

Total_people(:);

%%%%%%%%%%%%% Flow characteristics
for i=1:32
Pop_p((i-1)*49*10+1:(i-1)*49*10+49*10,:)=repmat(CO2_PC3(i,:),49*10,1); % Pop by product by age
end
Pop_p(:);


%%% Flow pop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Backward leakage pop
for j=1:3
    for k=1:4
    CO2_ppp=CO2_PC3(:,(j-1)*4+1:(j-1)*4+4)';
    CO3_ppp((k-1)*10*49+1:(k-1)*10*49+10*49,:)=repmat(CO2_ppp(k,:),10*49,1);
    CO4_ppp((j-1)*size(CO3_ppp(:),1)+1:(j-1)*size(CO3_ppp(:),1)+size(CO3_ppp(:),1),1)=CO3_ppp(:);
end
end



ppp=CO2_PC3(:,1:4);
for i=1:32
    for j=1:4
Pop_p1((j-1)*10+1:(j-1)*10+10,i)=repmat(ppp(i,j),10,1); % Pop by product by age
    end
end
Pop_p1(:);

%%%%%%Layout
Country1=repmat(Country,4,1);
Country2=repmat(Country1,3,1);

for i=1:size(Year,2)
Year1((i-1)*size(Country,1)*4+1:(i-1)*size(Country,1)*4+4*size(Country,1),1)=repmat(Year(i),4*size(Country,1),1);
end

Age=[1;2;3;4];
    for j=1:size(Age,1)
       Age1((j-1)*32+1:(j-1)*32+32,1)=repmat(Age(j),32,1);
    end
repmat(Age1,3,1);

%%%%%%Layout products
for i=1:size(Country,1)
Country1((i-1)*10+1:(i-1)*10+10,1)=repmat(Country(i),10,1);
end
Country2=repmat(Country1,4,1);
Country3=repmat(Country2,3,1);

for i=1:size(Year,2)
Year1((i-1)*10*32*4+1:(i-1)*10*32*4+10*32*4,1)=repmat(Year(i),10*32*4,1);
end

for j=1:size(Age,1)
       Age1((j-1)*32*10+1:(j-1)*32*10+32*10,1)=repmat(Age(j),32*10,1);
end
repmat(Age1,3,1);

%%%product layout
repmat(repmat(repmat(Name_prod,32,1),4,1),3,1);

%%%% Type
[~,~,Name_type]=xlsread('Bridge_product.xlsx','Typec');;%%% Copy from EXCEL
for i=1:32
   Name_type2((i-1)*10+1:(i-1)*10+10,1)=repmat(Name_type(i,2),10,1); 
end
repmat(repmat(Name_type2,4,1),3,1);




%%%%%%%%%%%%%%%%%%% Regression

%%% per capita asset
%%% per capita expenditure

Asset=xlsread('regression','Read');

for i=1:3
for j=1:4
   Asset_1((i-1)*32+1:(i-1)*32+32,j)=Asset(:,(i-1)*4+j);
end
end

EX_Per0=CO2_PC5./CO2_PC3_1;


%%%Total expenditure elasticity
for i=1:3   

EX_Pert((i-1)*32+1:(i-1)*32+32,:)=CO2_PC2(:,(i-1)*4+1:(i-1)*4+4);

end

 for k=1:10
   for p=1:32
       EX_Per(p,:)=EX_Per0((p-1)*10+k,:);
   end
       
for i=1:3
for j=1:4
   EX_Per_2((i-1)*32+1:(i-1)*32+32,j)=EX_Per(:,(i-1)*4+j);
end
end
   
   
    
    EX_Per1.([Name_prod{k}])=EX_Per_2;
    EX_Per1.(['total'])=EX_Pert;
 end
 Name_prod{11}='total';
EX_TOTOAL=EX_Per1.(['total']);
for k=1:11
     y=EX_Per1.([Name_prod{k}]);
     
for i=1:4
       
    x1=EX_TOTOAL(:,i);
    x1(find(x1<0))=0.00000000000000000001;
    
     x2=Asset_1(:,i);
    x2(find(x1<0))=0.000000000000001;
    
    
    X = [ones(size(x1)),log(x1),log(x2)];

   [b1((k-1)*3+1:(k-1)*3+3,i),bint((k-1)*3+1:(k-1)*3+3,(i-1)*3+1:(i-1)*3+2),~,~,stats((k-1)*5+i,:)] = regress(log(y(:,i)),X) ;

   rel=regstats(log(y(:,i)),X(:,2:end),'linear',{'adjrsquare','tstat','fstat'});
   
   T_reg((k-1)*3+1:(k-1)*3+3,i)=rel.tstat.pval;
   
end

end



%%%%%%% Asset Elasticity
for i=1:3   

EX_Pert((i-1)*32+1:(i-1)*32+32,:)=CO2_PC2(:,(i-1)*4+1:(i-1)*4+4);

end

 for k=1:10
   for p=1:32
       EX_Per(p,:)=EX_Per0((p-1)*10+k,:);
   end
       
for i=1:3
for j=1:4
   EX_Per_2((i-1)*32+1:(i-1)*32+32,j)=EX_Per(:,(i-1)*4+j);
end
end
   
   
    
    EX_Per1.([Name_prod{k}])=EX_Per_2;
    EX_Per1.(['total'])=EX_Pert;
 end
 Name_prod{11}='total';

for k=1:11
     y=EX_Per1.([Name_prod{k}]);
     
for i=1:4
       
    x1=Asset_1(:,i);
    x1(find(x1<0))=0.00000000000000000000000000000000001;
    
    X = [ones(size(x1)),log(x1)];

   [b1((k-1)*3+1:(k-1)*3+2,i),bint((k-1)*3+1:(k-1)*3+2,(i-1)*3+1:(i-1)*3+2),~,~,stats((k-1)*5+i,:)] = regress(log(y(:,i)),X) ;

   rel=regstats(log(y(:,i)),X(:,2:end),'linear',{'adjrsquare','tstat','fstat'});
   
   T_reg((k-1)*3+1:(k-1)*3+2,i)=rel.tstat.pval;
   
end

end
