%% function of gettting the radical Isomerization reaction rates
function [IsomeList, Isomerate]=getIsoRates_bondM(ConnM, BondM, radicalI, T)

singlebondneigh=[]; %n*3 matrix; first column: radical position second column:distance; third column:neutral position;

for r=1:length(radicalI)
    currentatom=radicalI(r);
    sgbn=[];
    sgbn_loop=[];
    for dist=1:6           %find the neighbours
        neighs=[];
        count=0;
        for c=1:length(currentatom)
            newneigh=ConnM(currentatom(c), BondM(currentatom(c),:)==1);
            newneigh=newneigh(not(ismember(newneigh, [sgbn_loop; radicalI(r);neighs])));
            neighs=[neighs; newneigh'];
%             for nei=1:length(newneigh)
%                 neighs(count+nei)=newneigh(nei);
%             end
%             count=count+nei+1;
        end
%         neighs(neighs==0)=[];
%         neighs=unique(neighs);
%         neighs(neighs==radicalI(r))=[];
        
%         neighs=neighs(not(ismember(neighs, sgbn_loop)));
        sgbn_loop=[sgbn_loop; neighs];
        if neighs
            if dist>2
            for n=neighs'
                if sum(BondM(n,:))<4 && sum(radicalI==n)==0       %if there's still hydrogen
                     sgbn=[sgbn; radicalI(r), dist, n];
                end
            end
            end
        end
        currentatom=neighs;
    end
    singlebondneigh=[singlebondneigh;sgbn];
end


if singlebondneigh
    IsomeList=singlebondneigh; %list of radical isomerization;
    Isomerate=zeros(size(IsomeList,1),1);  %list of radical isomerization rate;
    isorate3=7.54e5*T^1.75*exp(-74893.6/8.314/T);
    isorate4=2.62e5*T^1.62*exp(-46442.4/8.314/T);
    isorate5=2.58e4*T^1.67*exp(-42677.4/8.314/T);
    isorate6=1.062e3*T^1.67*exp(-42677.4/8.314/T);
    Isomerate(IsomeList(:,2)==3)=isorate3;
    Isomerate(IsomeList(:,2)==4)=isorate4;
    Isomerate(IsomeList(:,2)==5)=isorate5;
    Isomerate(IsomeList(:,2)==6)=isorate6;
else
    IsomeList=[];
    Isomerate=0;
end

end
