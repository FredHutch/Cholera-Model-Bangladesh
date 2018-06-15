% run vaccination campaigns with calibrated parameter sets

clear all;
global a b d de ga id im 
global k la m nu om p q th ve

%initiate rand and randn
rng('shuffle')

col='r r:b b:g g:m m:k k:';

titleComp={'Susceptibles','Vaccinated','Symptomatic ','Asymptomatic','Recovered', ...
        'Cholera in water'};
titleVac={'Proportional vaccination: all ages','Proportional vaccination: 1-14', ...
    'Proportional vaccination: 1-4 ','Fixed number of vaccine doses: all ages',...
    'Fixed number of vaccine doses: 1-14'};

fname='M66';
load (['best' fname '_97']);
b(:,4:6)=[];
load(['vac' fname 'a_100camp'],'storePar');
a=[1/2 1/3 1/10 0]/360; % exchange between age groups
ga=zeros(4,3);
ga(:,1)=par(1); % 5 days, andrews11
ga(:,2)=par(2); % 5 days, andrews11
ga(:,3)=par(3); % 3 years, dennis

m=par(4:7);  % mortality rates by age (matlab 1997)
d=ones(4,1)*par(8); % mortality rates due to Cholera
                      % 2% (mortality) / 5 days from table A2 andrews11 
k=par(9); % kL from andrews 11
om=zeros(4,2);
om(:,1)=par(10); %based on 15l per person in Matlab, andrews11
om(:,2)=par(11); %1/100 andrews11
th=par(12); % 30 days, andrews11
nu=zeros(4,4);
p=ones(4,1)*par(15); %andrews11
q=ones(4,3)*par(16); % guess
r=par(17);
de=par(18:21);
InitP=par(22);
InitI=par(23);
la=[par(24) 0 0 0]; % entry rates: guess

%ve=[0.38 0.85 0.69]; % (vacM6_61) vaccine efficacy in reducing infections (VEs) 1-preschoolers, 2-school age, 3-adults
ve=[0.42 0.68 0.74]; % (vacM6_62 from Bhattacharya et al 2013) vaccine efficacy in reducing infections (VEs) 1-preschoolers, 2-school age, 3-adults
b(1,1:3)=ageSusc(1)*b(4,1:3);
b(2,1:3)=ageSusc(2)*b(4,1:3);
b(3,1:3)=ageSusc(3)*b(4,1:3);


ageDist=[4.2 6.7 21.5 67.6]; % population fractions by age groups 
infDist=[13.9 16.1 23.3 46.7]; % infections distribution by age groups 

simN=100; % number of simulations
yrN=10; %number of years to simulate
period=30; % days, monthly accounting
revacTime=3; % period for revaccination (years)
vacCover=[70 70 70 55]; % percentage coverage
imRate=[0 0.1 0.25]/360; % immigration rate
vacStrategy=3; % testing 3 strategies (vacc 1-4, 1-14 or 1+)
scN=vacStrategy*length(imRate); % number of scenarios 

vacFrac=[0.5 1 0 0; % initially vaccinated 1-4y old
           0 0 1 0; % initially vaccinated 5-14y old
           0 0 0 1]; % initially vaccinated 15+y old
vacType=1; % 1-proportional campaign, 2-fixed doses campaign, 3-continuous

popDist=zeros(simN*(scN+1)*(yrN+1),57); % store population disctribution
 

for sim=1:simN
    sim
    idM=zeros(12*yrN,16); %monthly epidemic distribution  
    im=0; %immigration (exchange) rate

    sc=1; %current scenario
    %run without vaccination
    y1=zeros(1,57);
    y1([1 2 5:7 12 13 16:18 23 24 27:29 34 35 38:40 45])=y0(1:21);
    popDist((sim-1)*(scN+1)*(yrN+1)+(sc-1)*(yrN+1)+1,:)=y1;
    for yr=1:yrN
        % from file
        de([1 3 4])=storePar(sim,3*yr-2:3*yr);
        
        
        for mn=1:12
            id=reshape(y1(1:44),11,4); %initial epidemic distribution
            id([2:4 8:11],:)=[];
            id=id./repmat(sum(id),4,1);
            id=id';
            idM(12*(yr-1)+mn,:)=reshape(id,1,16);
            [T,Y] = ode45(@M57f,[mn-1 mn]*period,y1); %frequency dependence
            n=size(Y,1);
            y1=Y(n,:);
        end
        popDist((sim-1)*(scN+1)*(yrN+1)+(sc-1)*(yrN+1)+yr+1,:)=y1;
         
    end
   
    %vaccination scenarios
    cov=vacCover;
    InitVF=vacFrac.*repmat(cov,3,1)/100; % initially vaccinated fraction
            
    for vs=1:vacStrategy
        for im=imRate %immigration (exchange) rate
            sc=sc+1;
            y1=zeros(1,57);
            y1([1 2 5:7 12 13 16:18 23 24 27:29 34 35 38:40 45])=y0(1:21);
            
            %initial vaccination
            if (vs==1) %vaccinate only 1-4
                y1(2:11:44)=InitVF(1,:).*(y1(1:11:44)+y1(7:11:44));
                y1(7:11:44)=(1-InitVF(1,:)).*y1(7:11:44);
                y1(1:11:44)=(1-InitVF(1,:)).*y1(1:11:44);
            elseif (vs==2) %vaccinate only 1-14
                y1(2:11:44)=InitVF(1,:).*(y1(1:11:44)+y1(7:11:44));
                y1(3:11:44)=InitVF(2,:).*(y1(1:11:44)+y1(7:11:44));
                y1(7:11:44)=(1-sum(InitVF(1:2,:))).*y1(7:11:44);
                y1(1:11:44)=(1-sum(InitVF(1:2,:))).*y1(1:11:44);
            else % vaccinate all 1+
                y1(2:11:44)=InitVF(1,:).*(y1(1:11:44)+y1(7:11:44));
                y1(3:11:44)=InitVF(2,:).*(y1(1:11:44)+y1(7:11:44));
                y1(4:11:44)=InitVF(3,:).*(y1(1:11:44)+y1(7:11:44));
                y1(7:11:44)=(1-sum(InitVF)).*y1(7:11:44);
                y1(1:11:44)=(1-sum(InitVF)).*y1(1:11:44);
            end
            popDist((sim-1)*(scN+1)*(yrN+1)+(sc-1)*(yrN+1)+1,:)=y1;
            
            for yr=1:yrN
                
                de([1 3 4])=storePar(sim,3*yr-2:3*yr);
                
                for mn=1:12
                    id=reshape(idM(12*(yr-1)+mn,:),4,4);
                    [T,Y] = ode45(@M57f,[mn-1 mn]*period,y1); %frequency dependence
                    n=size(Y,1);
                    y1=Y(n,:);
                end
                
                %campaign revaccination
                if (rem(yr,3)==0)
                    if (vs==1) %vaccinate only 1-4
                        %eligible for revaccination
                        revacc=y1(1:11:44)+y1(2:11:44)+y1(3:11:44)+y1(4:11:44)+y1(7:11:44)+y1(10:11:44)+y1(11:11:44);
                                               
                        y1(1:11:44)=(1-InitVF(1,:)).*(y1(1:11:44)+y1(2:11:44)+y1(3:11:44)+y1(4:11:44)+y1(11:11:44));
                        y1(7:11:44)=(1-InitVF(1,:)).*(y1(7:11:44)+y1(10:11:44));
                        y1(10:11:44)=0;
                        y1(11:11:44)=0;
                        
                        y1(2:11:44)=InitVF(1,:).*revacc;
                        y1(3:11:44)=0;
                        y1(4:11:44)=0;
                        
                    elseif (vs==2) %vaccinate only 1-14
                                           
                        revacc=y1(1:11:44)+y1(2:11:44)+y1(3:11:44)+y1(4:11:44)+y1(7:11:44)+y1(10:11:44)+y1(11:11:44);
                        y1(1:11:44)=(1-sum(InitVF(1:2,:))).*(y1(1:11:44)+y1(2:11:44)+y1(3:11:44)+y1(4:11:44)+y1(11:11:44));
                        y1(7:11:44)=(1-sum(InitVF(1:2,:))).*(y1(7:11:44)+y1(10:11:44));
                        y1(10:11:44)=0;
                        y1(11:11:44)=0;
                        
                        y1(2:11:44)=InitVF(1,:).*revacc;
                        y1(3:11:44)=InitVF(2,:).*revacc;
                        y1(4:11:44)=0;
                    else % vaccinate all 1+
                        revacc=y1(1:11:44)+y1(2:11:44)+y1(3:11:44)+y1(4:11:44)+y1(7:11:44)+y1(10:11:44)+y1(11:11:44);
                        y1(1:11:44)=(1-sum(InitVF)).*(y1(1:11:44)+y1(2:11:44)+y1(3:11:44)+y1(4:11:44)+y1(11:11:44));
                        y1(7:11:44)=(1-sum(InitVF)).*(y1(7:11:44)+y1(10:11:44));
                        y1(10:11:44)=0;
                        y1(11:11:44)=0;
                        
                        y1(2:11:44)=InitVF(1,:).*revacc;
                        y1(3:11:44)=InitVF(2,:).*revacc;
                        y1(4:11:44)=InitVF(3,:).*revacc;
                    end
                else
                     %revaccinate infants
                    revacc=y1(1)+y1(2)+y1(3)+y1(4)+y1(7)+y1(10)+y1(11);
                    
                    y1(1)=(1-InitVF(1,1))*(y1(1)+y1(2)+y1(3)+y1(4)+y1(11));
                    y1(7)=(1-InitVF(1,1))*(y1(7)+y1(10));
                    y1(10)=0;
                    y1(11)=0;
                    
                    y1(2)=InitVF(1,1)*revacc;
                    y1(3)=0;
                    y1(4)=0;
                    
                end

                popDist((sim-1)*(scN+1)*(yrN+1)+(sc-1)*(yrN+1)+yr+1,:)=y1;
                
                
            end
            
        end
        
    end
    
end

save(['vac' fname '_62'],'par','b','ageSusc','storePar','popDist','idM','imRate');


    
    
    
