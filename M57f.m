% cholera model with asymptomatic infections and leaky vaccine,
% differentiated by age and seasonal forcing, frequency-dependent transmission 
% with completely separate vaccinated and 3 different vaccine regimens 
% to keep track if vaccinated as infant, child or adult
% added migration (in->unvaccinated, out->all)

function z=M57f(t,x)
global a b d de ga id im k la m nu om p q th ve

z=zeros(57,1);
N=sum(reshape(x(1:44),11,4)); % total population
bb=b(:,3)*(1+de(1)*(rem(t,360)<de(2))...             
   +de(3)*((rem(t,360)>de(4)+150) & (rem(t,360)<de(4)+210))); % stepwise seasonal forcing (double dip)

for ag=0:3
% susceptible unvaccinated ages (0-1y, 2-4y,5-14y, 15+y)
z(11*ag+1)=la(ag+1)*sum(N)+im*id(ag+1,1)*N(ag+1)+ga(ag+1,3)*x(11*ag+7)-(sum(nu(ag+1,1:3))+m(ag+1)+a(ag+1)+im)*x(11*ag+1)...
         -(b(ag+1,1)*sum(x([5:11:44 8:11:44]))/sum(N)+b(ag+1,2)*sum(x([6:11:44 9:11:44]))/sum(N)+bb(ag+1)*x(45)/(k+x(45)))*x(11*ag+1);
     
% vaccinated as infants(1-4y) - NO REVACCINATION !!!
z(11*ag+2)=nu(ag+1,1)*(x(11*ag+1)+x(11*ag+7))-(m(ag+1)+nu(ag+1,4)+a(ag+1)+im)*x(11*ag+2)...
         -(1-ve(1))*(b(ag+1,1)*sum(x([5:11:44 8:11:44]))/sum(N)+b(ag+1,2)*sum(x([6:11:44 9:11:44]))/sum(N)+bb(ag+1)*x(45)/(k+x(45)))*x(11*ag+2);

% vaccinated as children(5-14y) - NO REVACCINATION !!!
z(11*ag+3)=nu(ag+1,2)*(x(11*ag+1)+x(11*ag+7))-(m(ag+1)+nu(ag+1,4)+a(ag+1)+im)*x(11*ag+3)...
         -(1-ve(2))*(b(ag+1,1)*sum(x([5:11:44 8:11:44]))/sum(N)+b(ag+1,2)*sum(x([6:11:44 9:11:44]))/sum(N)+bb(ag+1)*x(45)/(k+x(45)))*x(11*ag+3);

% vaccinated as adult(15+y) - NO REVACCINATION !!!
z(11*ag+4)=nu(ag+1,3)*(x(11*ag+1)+x(11*ag+7))-(m(ag+1)+nu(ag+1,4)+a(ag+1)+im)*x(11*ag+4)...
         -(1-ve(3))*(b(ag+1,1)*sum(x([5:11:44 8:11:44]))/sum(N)+b(ag+1,2)*sum(x([6:11:44 9:11:44]))/sum(N)+bb(ag+1)*x(45)/(k+x(45)))*x(11*ag+4);

% unvaccinated symptomatic infected ages (0-1y, 2-4y,5-14y, 15+y)
z(11*ag+5)=im*id(ag+1,2)*N(ag+1)+p(ag+1)*(b(ag+1,1)*sum(x([5:11:44 8:11:44]))/sum(N)+b(ag+1,2)*sum(x([6:11:44 9:11:44]))/sum(N)+bb(ag+1)*x(45)/(k+x(45)))*x(11*ag+1)...
        -(m(ag+1)+d(ag+1)+ga(ag+1,1)+a(ag+1)+im)*x(11*ag+5);

% unvaccinated asymptomatic infected ages (0-1y, 2-4y,5-14y, 15+y)
z(11*ag+6)=im*id(ag+1,3)*N(ag+1)+(1-p(ag+1))*(b(ag+1,1)*sum(x([5:11:44 8:11:44]))/sum(N)+b(ag+1,2)*sum(x([6:11:44 9:11:44]))/sum(N)+bb(ag+1)*x(45)/(k+x(45)))*x(11*ag+1)...
         -(m(ag+1)+ga(ag+1,2)+a(ag+1)+im)*x(11*ag+6);

% unvaccinated recovered ages (0-1y, 2-4y,5-14y, 15+y)
z(11*ag+7)=im*id(ag+1,4)*N(ag+1)+ga(ag+1,1)*x(11*ag+5)+ga(ag+1,2)*x(11*ag+6)-(sum(nu(ag+1,1:3))+m(ag+1)+ga(ag+1,3)+a(ag+1)+im)*x(11*ag+7);


% vaccinated symptomatic infected ages (0-1y, 2-4y,5-14y, 15+y)
z(11*ag+8)=p(ag+1)*(b(ag+1,1)*sum(x([5:11:44 8:11:44]))/sum(N)+b(ag+1,2)*sum(x([6:11:44 9:11:44]))/sum(N)+bb(ag+1)*x(45)/(k+x(45)))*x(11*ag+11)...
    +q(ag+1,1)*(1-ve(1))*(b(ag+1,1)*sum(x([5:11:44 8:11:44]))/sum(N)+b(ag+1,2)*sum(x([6:11:44 9:11:44]))/sum(N)+bb(ag+1)*x(45)/(k+x(45)))*x(11*ag+2)...
    +q(ag+1,2)*(1-ve(2))*(b(ag+1,1)*sum(x([5:11:44 8:11:44]))/sum(N)+b(ag+1,2)*sum(x([6:11:44 9:11:44]))/sum(N)+bb(ag+1)*x(45)/(k+x(45)))*x(11*ag+3)...
    +q(ag+1,3)*(1-ve(3))*(b(ag+1,1)*sum(x([5:11:44 8:11:44]))/sum(N)+b(ag+1,2)*sum(x([6:11:44 9:11:44]))/sum(N)+bb(ag+1)*x(45)/(k+x(45)))*x(11*ag+4)...
 -(m(ag+1)+d(ag+1)+ga(ag+1,1)+a(ag+1)+im)*x(11*ag+8);

% vaccinated asymptomatic infected ages (0-1y, 2-4y,5-14y, 15+y)
z(11*ag+9)=(1-p(ag+1))*(b(ag+1,1)*sum(x([5:11:44 8:11:44]))/sum(N)+b(ag+1,2)*sum(x([6:11:44 9:11:44]))/sum(N)+bb(ag+1)*x(45)/(k+x(45)))*x(11*ag+11)...
    +(1-q(ag+1,1))*(1-ve(1))*(b(ag+1,1)*sum(x([5:11:44 8:11:44]))/sum(N)+b(ag+1,2)*sum(x([6:11:44 9:11:44]))/sum(N)+bb(ag+1)*x(45)/(k+x(45)))*x(11*ag+2)...
    +(1-q(ag+1,2))*(1-ve(2))*(b(ag+1,1)*sum(x([5:11:44 8:11:44]))/sum(N)+b(ag+1,2)*sum(x([6:11:44 9:11:44]))/sum(N)+bb(ag+1)*x(45)/(k+x(45)))*x(11*ag+3)...
    +(1-q(ag+1,3))*(1-ve(3))*(b(ag+1,1)*sum(x([5:11:44 8:11:44]))/sum(N)+b(ag+1,2)*sum(x([6:11:44 9:11:44]))/sum(N)+bb(ag+1)*x(45)/(k+x(45)))*x(11*ag+4)...
 -(m(ag+1)+ga(ag+1,2)+a(ag+1)+im)*x(11*ag+9);

% vaccinated recovered ages (0-1y, 2-4y,5-14y, 15+y)
z(11*ag+10)=ga(ag+1,1)*x(11*ag+8)+ga(ag+1,2)*x(11*ag+9)-(nu(ag+1,1)+m(ag+1)+ga(ag+1,3)+a(ag+1)+im)*x(11*ag+10);

% vaccinated who lost protection ages (0-1y, 2-4y,5-14y, 15+y)
z(11*ag+11)=nu(ag+1,4)*sum(x(11*ag+2:11*ag+4))+ga(ag+1,3)*x(11*ag+10)-(nu(ag+1,1)+m(ag+1)+a(ag+1)+im)*x(11*ag+11)...
        -(b(ag+1,1)*sum(x([5:11:44 8:11:44]))/sum(N)+b(ag+1,2)*sum(x([6:11:44 9:11:44]))/sum(N)+bb(ag+1)*x(45)/(k+x(45)))*x(11*ag+11);

%aging
if (ag>0)
  z(11*ag+1:11*ag+11)=z(11*ag+1:11*ag+11)+a(ag)*x(11*ag-10:11*ag);
    
end

% cumulative new infections unvaccinated (symptomatic) ages (0-1y, 2-4y,5-14y, 15+y)
z(3*ag+46)=p(ag+1)*(b(ag+1,1)*sum(x([5:11:44 8:11:44]))/sum(N)+b(ag+1,2)*sum(x([6:11:44 9:11:44]))/sum(N)+bb(ag+1)*x(45)/(k+x(45)))*x(11*ag+1);

% cumulative new infections (asymptomatic) ages (0-1y, 2-4y,5-14y, 15+y)
z(3*ag+47)=(1-p(ag+1))*(b(ag+1,1)*sum(x([5:11:44 8:11:44]))/sum(N)+b(ag+1,2)*sum(x([6:11:44 9:11:44]))/sum(N)+bb(ag+1)*x(45)/(k+x(45)))*(x(11*ag+1)+x(11*ag+11))...
    +(1-q(ag+1,1))*((1-ve(1))*b(ag+1,1)*sum(x([5:11:44 8:11:44]))/sum(N)+(1-ve(1))*b(ag+1,2)*sum(x([6:11:44 9:11:44]))/sum(N)+(1-ve(1))*bb(ag+1)*x(45)/(k+x(45)))*x(11*ag+2)...
    +(1-q(ag+1,2))*((1-ve(2))*b(ag+1,1)*sum(x([5:11:44 8:11:44]))/sum(N)+(1-ve(2))*b(ag+1,2)*sum(x([6:11:44 9:11:44]))/sum(N)+(1-ve(2))*bb(ag+1)*x(45)/(k+x(45)))*x(11*ag+3)...
    +(1-q(ag+1,3))*((1-ve(3))*b(ag+1,1)*sum(x([5:11:44 8:11:44]))/sum(N)+(1-ve(3))*b(ag+1,2)*sum(x([6:11:44 9:11:44]))/sum(N)+(1-ve(3))*bb(ag+1)*x(45)/(k+x(45)))*x(11*ag+4);
 
% cumulative new infections vaccinated (symptomatic) ages (0-1y, 2-4y,5-14y, 15+y)
z(3*ag+48)=p(ag+1)*(b(ag+1,1)*sum(x([5:11:44 8:11:44]))/sum(N)+b(ag+1,2)*sum(x([6:11:44 9:11:44]))/sum(N)+bb(ag+1)*x(45)/(k+x(45)))*x(11*ag+11)...
    +q(ag+1,1)*((1-ve(1))*b(ag+1,1)*sum(x([5:11:44 8:11:44]))/sum(N)+(1-ve(1))*b(ag+1,2)*sum(x([6:11:44 9:11:44]))/sum(N)+(1-ve(1))*bb(ag+1)*x(45)/(k+x(45)))*x(11*ag+2)...
    +q(ag+1,2)*((1-ve(2))*b(ag+1,1)*sum(x([5:11:44 8:11:44]))/sum(N)+(1-ve(2))*b(ag+1,2)*sum(x([6:11:44 9:11:44]))/sum(N)+(1-ve(2))*bb(ag+1)*x(45)/(k+x(45)))*x(11*ag+3)...
    +q(ag+1,3)*((1-ve(3))*b(ag+1,1)*sum(x([5:11:44 8:11:44]))/sum(N)+(1-ve(3))*b(ag+1,2)*sum(x([6:11:44 9:11:44]))/sum(N)+(1-ve(3))*bb(ag+1)*x(45)/(k+x(45)))*x(11*ag+4);

end


% bacteria in environment
z(45)=om(:,1)'*(x(5:11:44)+x(8:11:44))+om(:,2)'*(x(6:11:44)+x(9:11:44))-th*x(45);





