%% function of getting all radical reactions
function [BetaSciList, BetaSciRate, CapRate, TerminationList, TerminationRate]=getRadicalReactions_bondM(ConnM, BondM, radicalI, Cmarker, Dmarker, molmass, atomorder, T)
            BetaSciList=[];        %setting capping reaction rates%Calculate rates of beta scission row1:graph # row2:radical # row3:begin atom #row4: end atom
            CapRate=zeros(length(radicalI),1);  %capping reaction rates
            CapRate(:)=100;         
            radicaldeg=sum(BondM(radicalI,:)~=0,2);
            CapRate(radicaldeg==0)=CapRate(radicaldeg==0)*4;
            CapRate(radicaldeg==1)=CapRate(radicaldeg==1)*2;
            
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
            
            BetaSciRate=zeros(size(BetaSciList,1),1);
            BetaSciRate(:)=1.22e13*T^0.31*exp(-117916/(8.314*T));  %setting beta Scission rates
            if BetaSciList
            CKIE=exp((-1.51613E06/T^2 - 8.73339E03/T - 1.50983E+01)/1000); %Xiao et al., for ethane
            DKIE=1.2.*exp(-327*4.184/8.314/T);    %Ni et al for homolytic cleavage
            
            BetaSciRate(Cmarker(BetaSciList(:,4)))=BetaSciRate(Cmarker(BetaSciList(:,4)))*CKIE;
            BetaSciRate(Dmarker(BetaSciList(:,4)))=BetaSciRate(Dmarker(BetaSciList(:,4)))*DKIE;
            end
            TerminationList=[];%List all of the termination reactions row1:graph1 row2:radical1 row3:graph3 row4:radical4
            for i=1:length(CapRate)
                for j=i:length(CapRate)
                    if atomorder(radicalI(i))=='C'
                        neigh=ConnM(radicalI(i),ConnM(radicalI(i),:)~=0);
                        if not(ismember(radicalI(j),neigh)) && radicalI(j)~=radicalI(i)
                            TerminationList=[TerminationList; 1, radicalI(i), 1, radicalI(j)];
                        end
                    end
                end
            end
            TerminationRate=zeros(size(TerminationList,1),1);
            TerminationRate(:)=3e6*2/(molmass/10^6)*0.05;
end
