%% function of getting all radical reactions
%Added KIE for both sites of beta scission
%Added degree-specific beta-scission rates
function [BetaSciList, BetaSciRate, CapRate, CapRate_w, TerminationList, TerminationRate]=getRadicalReactions_bondM_BSKIE(ConnM, BondM, radicalI, Cmarker, Dmarker, molmass, atomorder, T)
            BetaSciList=[];        %setting capping reaction rates%Calculate rates of beta scission row1:graph # row2:radical # row3:begin atom #row4: end atom
            
            CapRate=zeros(length(radicalI),1);  %capping reaction rates; capping by organic
            CapRate_w=CapRate;  %capping by water hydrogen

%            CapRate_w(:)=100;         
            radicaldeg=sum(BondM(radicalI,:)~=0,2);
            
            ConcOrgH=3*10^6*0.015*(4246/55202)/2/2; %concentration of organicH in 2mol/m3; 1.5% is EFK80:C3849H4246O298
%             These are for attacking methylene groups
%             CapRate(radicaldeg==0)=ConcOrgH*1.6E-8*T^4.34*exp(-27656/8.314/T);
%             CapRate(radicaldeg==1)=ConcOrgH*2.16E-9*T^4.29*exp(-31532/8.314/T);
%             CapRate(radicaldeg==2)=ConcOrgH*2.57E-9*T^4.29*exp(-32269/8.314/T);
%             CapRate(radicaldeg==3)=ConcOrgH*1.43E-10*T^4.56*exp(-28507/8.314/T);
% %             These are for attacking methyl groups
%             CapRate(radicaldeg==0)=ConcOrgH*2E-6*T^3.57*exp(-32283/8.314/T);
%             CapRate(radicaldeg==1)=ConcOrgH*2.16E-9*T^4.29*exp(-31532/8.314/T);
%             CapRate(radicaldeg==2)=ConcOrgH*2.57E-9*T^4.29*exp(-32269/8.314/T);
%             CapRate(radicaldeg==3)=ConcOrgH*6.36E-9*T^4.34*exp(-40584/8.314/T);
%            These are for attacking methane
            CapRate(radicaldeg==0)=ConcOrgH*1.55E-8*T^4.34*exp(-43095/8.314/T);
            CapRate(radicaldeg==1)=ConcOrgH*8.64E-8*T^4.14*exp(-52551/8.314/T);
            CapRate(radicaldeg==2)=ConcOrgH*7.24E-10*T^4.40*exp(-45145/8.314/T);
            CapRate(radicaldeg==3)=ConcOrgH*6.16E-9*T^4.34*exp(-52718/8.314/T);
            
            ConcH2O=3*10^6/18*0.16;  %concentration of water in mol/m3; 0.16 is water content from Eagle Ford Shale
            %June-23-2021 
            CapRate_w(radicaldeg==0)=ConcH2O*6.4E-6*T^3.31*exp(-57250/8.314/T);
            CapRate_w(radicaldeg==1)=ConcH2O*3.4E0*T^1.44*exp(-84810/8.314/T);
            CapRate_w(radicaldeg==2)=ConcH2O*8.52E-9*T^3.99*exp(-84955/8.314/T);
            CapRate_w(radicaldeg==3)=ConcH2O*9.26E-6*T^3.06*exp(-93148/8.314/T);
            
            for i=1:length(radicalI)
%                     [edges, neigh]=outedges(BoG_mega{1},radicalI(i));
                    neigh=ConnM(radicalI(i), BondM(radicalI(i),:)==1);
                    for n=1:length(neigh)
                        if  atomorder(neigh(n))=='C'        %only for carbons
%                             [secondaryedges, secondaryneigh]=outedges(BoG_mega{1},neigh(n));
                            secondaryneigh=ConnM(neigh(n), BondM(neigh(n),:)==1);
                            for p=1:length(secondaryneigh)
                                if secondaryneigh(p)~=radicalI(i) && atomorder(secondaryneigh(p))=='C'
                                    BetaSciList=[BetaSciList; 1, radicalI(i), neigh(n), secondaryneigh(p)];
                                end
                            end
                        end
                    end

            end
            if BetaSciList

            BetaSciRate=zeros(size(BetaSciList,1),1);
            BetaSciRate(:)=1.22e13*T^0.31*exp(-117916/(8.314*T));
            
            CKIE=exp((-1.51613E06/T^2 - 8.73339E03/T - 1.50983E+01)/1000); %Xiao et al., for ethane
            DKIE=1.2.*exp(-327*4.184/8.314/T);    %Ni et al for homolytic cleavage
            
            BetaSciRate(Cmarker(BetaSciList(:,4)))=BetaSciRate(Cmarker(BetaSciList(:,4)))*CKIE;
            BetaSciRate(Cmarker(BetaSciList(:,3)))=BetaSciRate(Cmarker(BetaSciList(:,3)))*CKIE;
            BetaSciRate(Dmarker(BetaSciList(:,4)))=BetaSciRate(Dmarker(BetaSciList(:,4)))*DKIE;
            BetaSciRate(Dmarker(BetaSciList(:,3)))=BetaSciRate(Dmarker(BetaSciList(:,3)))*DKIE;
            else
                BetaSciRate=0;
            end
            TerminationList=[];%List all of the termination reactions row1:graph1 row2:radical1 row3:graph3 row4:radical4
            for i=1:length(CapRate_w)
                for j=i:length(CapRate_w)
                    if atomorder(radicalI(i))=='C'
                        neigh=ConnM(radicalI(i),ConnM(radicalI(i),:)~=0);
                        if not(ismember(radicalI(j),neigh)) && radicalI(j)~=radicalI(i)
                            TerminationList=[TerminationList; 1, radicalI(i), 1, radicalI(j)];
                        end
                    end
                end
            end
            TerminationRate=zeros(size(TerminationList,1),1);
            TerminationRate(:)=3e6*2/(molmass/10^6)*0.015;
end
