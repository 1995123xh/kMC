%% Making a kinetic MonteCarlo for simulating propane generation from
% Dec-2020 add beta scission
%%% kerogen molecules
%%%
clear all
Filename='oil';

disp('loading...')
load(Filename)

nmol=100;         %number of molecules in the graph. 1000-2000 atoms have highest efficiency
sz=1;         %number of graphs in the cluster
enrichfactor=50; %the factor of deuterium concentration
enrichfactor_13C=5;
crkratio=100;    %fraction of bonds going broken
stopat=0.5;      %percentage of completion
sims=500;        %number of MC simulations
rep=8;       %number of reps in the parfor loop
updatestepsize=100;  %step size (HomoC+BetaC+Termination) for refreshing the reaction rates due to temeprature change and structural change.
CapIso=30; %maximum isomerization step number
LowN=10; %minimum isomerization step number

ThermalProgram=[0 0.2 0.5 0.7 1; 500 500 500 500 500];  % A thermal program where first row is fraction to completion and second row is temeprature
inittemp=ThermalProgram(2,1);
stopat=stopat/100;

BoM_initial=BoM;    %transfering the original molecule somewhere else...
atoms_initial=atoms;

graphstable=cell(sims,rep);
gcount=0;

[RowI,ColumnI]=find(BoM~=0);
RowI0=RowI;
ColumnI0=ColumnI;
Ele=[BoM(BoM~=0)];
Ele0=Ele;
for i=1:nmol-1
    RowI=[RowI;RowI0+length(BoM)*i];
    ColumnI=[ColumnI;ColumnI0+length(BoM)*i];
    Ele=[Ele;Ele0];
    atoms=strcat(atoms,atoms_initial);
end
BoM=sparse(RowI,ColumnI, Ele);
BoG=graph(BoM, 'upper','omitselfloops');    %Construct a graph of bond orders
BoG0=BoG;           %the original molecule

ConnM=zeros(length(atoms),4);  %connection matrix: n*4 matrix, where each row is an atom and each column is a neighbor atom;...
%Note that this is a full-loop (not upper triangle), a bond shows
%up twice in the matrix
BondM=ConnM;                   %bonding matrix: n*4 matrix, where each row is an atom and each column is the bond weight
NodesList=BoG.Edges.EndNodes;
EdgeWeights=BoG.Edges.Weight;
pointer=zeros(length(atoms),1);
pointer(:)=1;
for i=1:length(NodesList)
    ConnM(NodesList(i,1),pointer(NodesList(i,1)))=NodesList(i,2);
    ConnM(NodesList(i,2),pointer(NodesList(i,2)))=NodesList(i,1);
    BondM(NodesList(i,1),pointer(NodesList(i,1)))=EdgeWeights(i);
    BondM(NodesList(i,2),pointer(NodesList(i,2)))=EdgeWeights(i);
    pointer(NodesList(i,1))=pointer(NodesList(i,1))+1;
    pointer(NodesList(i,2))=pointer(NodesList(i,2))+1;
end

ConnM0=ConnM;
BondM0=BondM;

% BoG_mega=cell(1,sz);
% BoG_mega(:)={BoG};
% BoG_mega0=BoG_mega;
%
% BoM_mega=cell(1,sz);
% BoM_mega(:)={BoM};
% BoM0=BoM;
% BoM_mega0=BoM_mega;
% for i=1:length(BoM_mega0)
%     BoM_mega0{i}=sparse(BoM_mega0{i});
% end

deg=degree(BoG0)';
deg_mega=zeros(sz, length(BoM));
for i=1:sz
    deg_mega(i,:)=deg;
end
deg0=deg;
deg_mega0=deg_mega;

ethane_13Cn=0;
ethane_D1=0;
ethane_D2=0;
ethane_Dn=0;

ethane_13C1=0;
ethane_13C2=0;

methane_n=0;
methane_D=0;
methane_13C=0;
methane_13CD=0;

for i=1:6                %creating output data matching structure for C3+ straight chain n alkanes
    cnout(i).len=i+2;
    singlelabel=cnout(i).len;   %singly-substituted species
    %   doublelabel=cnout(i).length*(cnout(i).length-1)/2;
    labelmatrix=zeros(singlelabel,cnout(i).len);
    for j=1:cnout(i).len
        labelmatrix(j,j)=1;
    end
    labelmatrix=[labelmatrix;zeros(1,cnout(i).len)];
    cnout(i).label=labelmatrix;       %ways of isotopic substitution
    cnm=zeros(cnout(i).len,cnout(i).len);
    for j=1:cnout(i).len-1
        cnm(j,j+1)=1;
    end
    cnout(i).BoM=cnm;
    cnout(i).BoG=graph(cnm, 'upper','omitselfloops');
    cnout(i).labelgraph=cell(size(labelmatrix,1),1);
    for j=1:length(cnout(i).labelgraph)
        cnout(i).labelgraph{j}=cnout(i).BoG;
        cnout(i).labelgraph{j}.Nodes.CLabel=labelmatrix(j,:)';
        cnout(i).labelgraph{j}.Nodes.DLabel=labelmatrix(j,:)';
    end
end
cnoutiso_D=zeros(sims,length(cnout),8+1);
cnoutiso_C=zeros(sims,length(cnout),8+1);


%% Constants
R0=155.76e-6*(1-100/1000)*enrichfactor;
R0_CAP=155.76e-6*(1-20/1000)*enrichfactor;   %water-sourced
R0_13C=0.01118*enrichfactor_13C;
F0=R0/(R0+1);
F0_13C=R0_13C/(R0_13C+1);
% T=273.15+180; %kelvin
R=8.314;
RH2=(1-836/1000)*R0;             %hydrogen isotopes of the hydrogen gas
aH=3;        %fractionation factor for capping hydrogen
molmass=12*sum(atoms=='C')+14*sum(atoms=='N')+16*sum(atoms=='O')+32*sum(atoms=='S')+sum(4-deg);
molmass=molmass*nmol;
disp('Initializing...')

%% Mark the aromaticity
arobonds=find(abs(BoG.Edges.Weight-1.5)<0.001);
for b=arobonds
    [beginA, endA] = findedge(BoG, b);
    atoms(beginA)='A';
    atoms(endA)='A';
end

%% Monte-Carlo
higherbonds=0;

ethane_D1_log=zeros(sims,1);
ethane_D2_log=zeros(sims,1);
ethane_Dn_log=zeros(sims,1);
ethane_13C1_log=zeros(sims,1);
ethane_13C2_log=zeros(sims,1);
ethane_13Cn_log=zeros(sims,1);

methane_D_log=zeros(sims,1);
methane_13C_log=zeros(sims,1);
methane_13CD_log=zeros(sims,1);
methane_n_log=zeros(sims,1);

hnumber=4-sum(BondM,2)';
dabundance=hnumber.*F0.*(1-F0).^(hnumber-1); % fixed a bug of no exponent
% parpool('local')
tic
% BoG_mega0=parallel.pool.Constant(BoG_mega0);
disp('Starts Monte-Carlo...');
propaneN=[];
% timeM=zeros(sims, rep); for calculating time
for sim=1:sims
    cnoutisotemp_D=cell(rep,1);   %pre-allocate these counters
    cnoutisotemp_C=cell(rep,1);
    %           mypool=gcp;
    %           addAttachedFiles(mypool,{'getRateMatrix.m','kinPara.m'})
    %          ticBytes(gcp);
    
    parfor iteration=1:rep
        %random number shuffler
        rng('shuffle');
        %data transfer
        ConnM=ConnM0;
        BondM=BondM0;
        %         BoG_mega=cell(1,sz);
        %         BoG_mega(:)={BoG};
        %                 BoM_mega=full(BoM_mega0);
        %         [beginN, numc]=find(BondM~=0);
        %         endN=beginN;
        %         for i=1:length(endN)
        %             endN(i)=ConnM(beginN(i), numc(i));
        %         end
        %         edgeweights=BondM(BondM~=0);
        BoM_begin=BoM;
        BoG_mega=cell(1,sz);
        BoG_mega{1}=graph(BoM_begin, 'upper', "omitselfloops");
        
        deg_mega=deg_mega0;
        atomorder=atoms;
        degorder=deg;
        %        BoMorder=BoM;
        ThermalProgramPar=ThermalProgram;
        %         cutedgelist=[];
        %         singlebonds=singlebonds0;
        
        b=0;
        cnoutiso_Dpar=zeros(length(cnout),8+1);
        cnoutiso_Cpar=zeros(length(cnout),8+1);
        cnoutorder=cnout;
        %         BoG=BoG0;
        hnumberPar=hnumber;
        
        %generating random isotope distribution
        cabundance=F0_13C; %
        randomarray=rand(size(hnumberPar));
        Dmarker=randomarray<dabundance;
        randomarray=rand(size(hnumberPar));
        Cmarker=randomarray<F0_13C;
        Cmarker(:,atoms~='C')=0;                   % 2019-06-15 fixed the error of wrong indexing here
        
        %radical=zeros(size(Cmarker));       %2020-12-04 added radical label;
        radicalI=[];
        T=ThermalProgramPar(2,1);
        %calculating rate
        [rate_mega, singlebonds]=getRateMatrix_bondM(ConnM, BondM, Dmarker, Cmarker, atomorder, T+273);
        totalbonds=length(singlebonds);
        homostep=updatestepsize;
        isostep=0;
        isotreshold=randi([LowN, CapIso])+1;
        for step=1:round(totalbonds*crkratio*stopat) % core loop for reactions
            
            T=interp1(ThermalProgramPar(1,:), ThermalProgramPar(2,:), step/totalbonds/crkratio);
            %[radicalG,radicalI]=find(radical>0);    %locate the radicals
            
            [BetaSciList, BetaSciRate, CapRate, CapRate_w, TerminationList, TerminationRate]=getRadicalReactions_bondM_BSKIEonly(ConnM, BondM, radicalI, Cmarker, Dmarker, molmass, atomorder, T);
            cumBetaSciRate=cumsum(BetaSciRate);
            cumCapRate=cumsum(CapRate);
            cumCapRate_w=cumsum(CapRate_w);
            cumTerminationRate=cumsum(TerminationRate);
            
            if isostep==isotreshold
                Isomerate=0;
                isostep=0;
                isotreshold=randi([LowN, CapIso])+1;
            else
                [IsomeList, Isomerate]=getIsoRates_bondM(ConnM, BondM, radicalI, T);
            end
            
            if homostep==updatestepsize         %when counter meets treshold, update homolytic cleavage rates
                
                %  [predic, KIE, KIE_13C]=kinPara(T);  %update kinetics based on temperature
                [rate_mega, singlebonds]=getRateMatrix_bondM(ConnM, BondM, Dmarker, Cmarker,atomorder, T);
                %  rate_mega(rate_mega~=0)=new_rate_mega(rate_mega~=0);
                sumrate=sum(sum(rate_mega));
                cumrate=cumsum(rate_mega(:)); % Generate a cumulative sum of rates
                homostep=1;
                
            end
            
            sumrate_total=sumrate+sum(BetaSciRate)+sum(TerminationRate)+sum(CapRate)+sum(CapRate_w)+sum(Isomerate);
            if sumrate_total==0
                % disp('No Reaction Found');
                break
            end
            dice=rand()*sumrate_total;   %selection: 4 types of reactions
            %timeM(sim, iteration)=timeM(sim,iteration)+1/sumrate_total*log(1/rand());
            if dice<sumrate              %homolytic cleavage
                dice=dice*sumrate/sumrate_total;
                i = find(cumrate > dice,1);
                index=ceil(i/sz); %the graph
                g=i-(index-1)*sz;
                
                %cut the bond
                %                 [A1, A2] = findedge(BoG,singlebonds(index,1));
                A1=singlebonds(index,2);
                A2=singlebonds(index,3);
                if BondM(A1, ConnM(A1,:)==A2)==1
                    %                     deg_mega(g,A1)=deg_mega(g,A1)-1;
                    %                     deg_mega(g,A2)=deg_mega(g,A2)-1;
                    %                     radical(g,A1)=radical(g,A1)+1;
                    %                     radical(g,A2)=radical(g,A2)+1;
                    
                    radicalI=[radicalI, A1, A2];
                    rate_mega(g, index)=0;         %clean that from the rate calculations
                    BondM(A1, ConnM(A1,:)==A2)=0;
                    BondM(A2, ConnM(A2,:)==A1)=0;
                    ConnM(A1, ConnM(A1,:)==A2)=0;
                    ConnM(A2, ConnM(A2,:)==A1)=0;
                    %                     BoG_mega{g}=rmedge(BoG_mega{g}, A1,A2);
                    % disp(['homocleavage' ,' ', num2str(A1),' ', num2str(A2)]);
                    homostep=homostep+1;
                end
            elseif dice<sumrate+sum(BetaSciRate)   %beta Scission
                dice=(dice-sumrate);
                i = find(cumBetaSciRate > dice,1);
                %                 BoM_mega{BetaSciList(i,1)}(BetaSciList(i,3),BetaSciList(i,4))=0;
                %                 BoM_mega{BetaSciList(i,1)}(BetaSciList(i,4),BetaSciList(i,3))=0;
                B2=BetaSciList(i,2);
                B3=BetaSciList(i,3);
                B4=BetaSciList(i,4);
                BondM(B3, ConnM(B3,:)==B4)=0;
                BondM(B4, ConnM(B4,:)==B3)=0;
                ConnM(B3, ConnM(B3,:)==B4)=0;
                ConnM(B4, ConnM(B4,:)==B3)=0;
                
                %                BoG_mega{BetaSciList(i,1)}=rmedge(BoG_mega{BetaSciList(i,1)}, BetaSciList(i,3),BetaSciList(i,4));
                %                radical(BetaSciList(i,1), BetaSciList(i,2))=radical(BetaSciList(i,1), BetaSciList(i,2))-1;
                %                radical(BetaSciList(i,1), BetaSciList(i,4))=radical(BetaSciList(i,1), BetaSciList(i,4))+1;
                radicalI(radicalI==BetaSciList(i,2))=[];
                radicalI=[radicalI, BetaSciList(i,4)];
                %                 BoM_mega{BetaSciList(i,1)}(BetaSciList(i,2),BetaSciList(i,3))=2;
                %                 BoM_mega{BetaSciList(i,1)}(BetaSciList(i,3),BetaSciList(i,2))=2;
                %                 newDouble=findedge(BoG_mega{BetaSciList(i,1)}, BetaSciList(i,2),BetaSciList(i,3));
                %                 BoG_mega{BetaSciList(i,1)}.Edges.Weight(newDouble)=2;
                BondM(B2, ConnM(B2,:)==B3)=2;
                BondM(B3, ConnM(B3,:)==B2)=2;
                % disp(['betaSci' ,' ', num2str(BetaSciList(i,3)),' ', num2str(BetaSciList(i,4))]);
                homostep=homostep+1;
            elseif dice<sumrate+sum(BetaSciRate)+sum(TerminationRate) %termination
                dice=dice-sumrate-sum(BetaSciRate);
                i = find(cumTerminationRate > dice,1);
                T2=TerminationList(i,2);
                T4=TerminationList(i,4);
                %      radical(TerminationList(i,1),TerminationList(i,2))=radical(TerminationList(i,1),TerminationList(i,2))-1;
                %      radical(TerminationList(i,3),TerminationList(i,4))=radical(TerminationList(i,3),TerminationList(i,4))-1;
                radicalI(radicalI==TerminationList(i,2))=[];
                radicalI(radicalI==TerminationList(i,4))=[];
                %                 BoM_mega{1}(TerminationList(i,2),TerminationList(i,4))=1;
                %                 BoM_mega{1}(TerminationList(i,4),TerminationList(i,2))=1;
                T2numc=find(BondM(T2,:)==0,1);
                T4numc=find(BondM(T4,:)==0,1);
                BondM(T2, T2numc)=1;
                BondM(T4, T4numc)=1;
                ConnM(T2, T2numc)=T4;
                ConnM(T4, T4numc)=T2;
                %   BoG_mega{1}=addedge(BoG_mega{1}, TerminationList(i,2),TerminationList(i,4),1);
                % disp(['Termination' ,' ', num2str(TerminationList(i,2)),' ', num2str(TerminationList(i,4))]);
                homostep=homostep+1;
            elseif dice<sumrate+sum(BetaSciRate)+sum(TerminationRate)+sum(CapRate)  %capping with organic hydrogen
                dice=(dice-sumrate-sum(BetaSciRate)-sum(TerminationRate)); %capping
                i = find(cumCapRate > dice,1);
                %   radical(radicalG(i),radicalI(i))=radical(radicalG(i),radicalI(i))-1;
                rIi=radicalI(i);
                
                TransferDice=randi(length(atomorder));   %determine the transfer reaction (where H comes from)
                
                Hatoms=4-sum(BondM(TransferDice,:));
                if sum(radicalI==TransferDice)==0 && Hatoms>0
                    if Dmarker(TransferDice)
                        newdice=rand();
                        if newdice<(Hatoms-1+1/aH)/Hatoms    %the last bracket does nothing in order to account for hydrogen isotope KIE
                            radicalI(i)=[];
                            radicalI=[radicalI, TransferDice];
                            % homostep=homostep+1;
                            %H migration
                            if newdice>(Hatoms-1)/(Hatoms)
                                Dmarker(rIi)=1;%D migration
                                Dmarker(TransferDice)=0;
                            end
                        end
                    else    %protium transfer
                        radicalI(i)=[];
                        radicalI=[radicalI, TransferDice];
                    end
                end
                
                % disp(['Capping' ,' ', num2str(rIi)]);
            elseif dice<sumrate+sum(BetaSciRate)+sum(TerminationRate)+sum(CapRate)+sum(CapRate_w) %capping with H from water
                dice=dice-sumrate-sum(BetaSciRate)-sum(TerminationRate)-sum(CapRate); %capping
                i = find(cumCapRate_w > dice,1);
                rIi=radicalI(i);
                radicalI(i)=[];
                % homostep=homostep+1;
                newdice=rand();
                if newdice<R0_CAP/aH/(R0_CAP/aH+1)
                    Dmarker(rIi)=1;
                end
                % disp(['Capping_w' ,' ', num2str(rIi)]);
            else
                dice=dice-sumrate-sum(BetaSciRate)-sum(TerminationRate)-sum(CapRate)-sum(CapRate_w);
                i=find(cumsum(Isomerate)>dice,1);
                % radical(1,IsomeList(i,1))=radical(1,IsomeList(i,1))-1;
                % radical(1,IsomeList(i,3))=radical(1,IsomeList(i,3))+1;
                if Dmarker(IsomeList(i,3))
                    Hatoms=4-sum(BondM(IsomeList(i,3),:));  %number of H atoms
                    newdice=rand();
                    if newdice<(Hatoms-1+1/aH)/Hatoms    %the last bracket does nothing in order to account for hydrogen isotope KIE
                        radicalI(radicalI==IsomeList(i,1))=[];
                        radicalI=[radicalI, IsomeList(i,3)];
                        % homostep=homostep+1;
                        isostep=isostep+1;
                        %H migration
                        if newdice>(Hatoms-1)/(Hatoms)
                            Dmarker(IsomeList(i,1))=1;%D migration
                            Dmarker(IsomeList(i,3))=0;
                        end
                    end
                else
                    radicalI(radicalI==IsomeList(i,1))=[];
                    radicalI=[radicalI, IsomeList(i,3)];
                    %homostep=homostep+1;
                    isostep=isostep+1;
                end
                
                %              disp(['Radical isomerization',' ', num2str(IsomeList(i,1)), ' ',num2str(IsomeList(i,3))])
                
            end
        end
        %Count hydrocarbons
        %         tables=cell(sz,1);
        [beginN, numc]=find(BondM~=0);
        endN=beginN;
        for i0=1:length(endN)
            endN(i0)=ConnM(beginN(i0), numc(i0));
        end
        edgeweights=BondM(BondM~=0);
        BoM_end=sparse(beginN, endN, edgeweights);
        BoG_mega{1}=graph(BoM_end, 'upper', "omitselfloops");
        
        for g0=1:sz
            %             tables{g}=BoG_mega{g}.Edges;
            splitg=conncomp(BoG_mega{g0},'OutputForm','cell'); %decompose eachgraph into sepearte parts
            for i1=1:length(splitg)
                if length(splitg{i1})<9          %ignore larger fragments
                    if length(splitg{i1})==2  %Find the ethanes
                        pnodes=splitg{i1};
                        if atomorder(pnodes)=='CC'
                            if BoM_end(pnodes(1), pnodes(2))==1 %single C-C bonds
                                Dlabel=[(Dmarker(g0, pnodes(1))), (Dmarker(g0, pnodes(2)))];
                                if or(isequal(Dlabel, [1, 0]), isequal(Dlabel,[0, 1]))
                                    ethane_D1=ethane_D1+1;
                                elseif isequal(Dlabel, [1, 1])
                                    ethane_D2=ethane_D2+1;
                                else
                                    ethane_Dn=ethane_Dn+1;
                                end
                                Clabel=[Cmarker(g0, pnodes(1)), Cmarker(g0, pnodes(2))];
                                if or(isequal(Clabel, [1, 0]), isequal(Clabel,[0, 1]))
                                    ethane_13C1=ethane_13C1+1;
                                elseif isequal(Clabel, [1, 1])
                                    ethane_13C2=ethane_13C2+1;
                                else
                                    ethane_13Cn=ethane_13Cn+1;
                                end
                            end
                        end
                    elseif length(splitg{i1})==1  %Find the methanes
                        pnodes=splitg{i1};
                        if atomorder(pnodes)=='C'
                            Dlabel=Dmarker(g0, pnodes);
                            Clabel=Cmarker(g0, pnodes);
                            if Dlabel
                                methane_D=methane_D+1;
                                if Clabel
                                    methane_13CD=methane_13CD+1;
                                end
                            end
                            if Clabel
                                methane_13C=methane_13C+1;
                            elseif not(Dlabel)
                                methane_n=methane_n+1;
                            end
                        end
                    else               %creating output for C3+ straight chain n alkanes
                        cn=length(splitg{i1})-2;   %
                        pnodes=splitg{i1};                           %
                        if all(atomorder(pnodes)=='C')                   %no hetero atoms
                            splitg_i=subgraph(BoG_mega{g0},splitg{i1});   %construct a seperate graph
                            Dlabel=[Dmarker(g0, pnodes)];
                            Clabel=[Cmarker(g0, pnodes)];
                            splitg_i.Nodes.DLabel=double(Dlabel');
                            splitg_i.Nodes.CLabel=double(Clabel');
                            if isisomorphic(cnoutorder(cn).BoG,splitg_i, 'EdgeVariables','Weight')
                                if cn==1
                                    propaneN=[propaneN; pnodes];
                                end
                                for isotopomer=1:size(cnoutorder(cn).label,1)
                                    if isequal(cnoutorder(cn).label(isotopomer,:),Dlabel)
                                        cnoutiso_Dpar(cn,isotopomer)=cnoutiso_Dpar(cn,isotopomer)+1;
                                    end
                                    if isequal(cnoutorder(cn).label(isotopomer,:),Clabel)
                                        cnoutiso_Cpar(cn,isotopomer)=cnoutiso_Cpar(cn,isotopomer)+1;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        cnoutisotemp_D{iteration}=cnoutiso_Dpar;
        cnoutisotemp_C{iteration}=cnoutiso_Cpar;
        %          graphstable{sim,iteration}=tables;
    end
    %     tocBytes(gcp);
    for iteration=1:rep
        m_D=cnoutisotemp_D{iteration};
        m_C=cnoutisotemp_C{iteration};
        for cn=1:length(cnout)
            for isotopomer=1:size(cnout(cn).label,1)
                cnoutiso_D(sim,cn,isotopomer)=cnoutiso_D(sim,cn,isotopomer)+m_D(cn,isotopomer);
                cnoutiso_C(sim,cn,isotopomer)=cnoutiso_C(sim,cn,isotopomer)+m_C(cn,isotopomer);
            end
        end
    end
    
    ethane_D1_log(sim)=ethane_D1;
    ethane_D2_log(sim)=ethane_D2;
    ethane_Dn_log(sim)=ethane_Dn;
    ethane_13C1_log(sim)=ethane_13C1;
    ethane_13C2_log(sim)=ethane_13C2;
    ethane_13Cn_log(sim)=ethane_13Cn;
    ethane_log=ethane_13Cn+ethane_13C1+ethane_13C2;  %total amount of ethane
    
    methane_D_log(sim)=methane_D;
    methane_13C_log(sim)=methane_13C;
    methane_13CD_log(sim)=methane_13CD;
    methane_n_log(sim)=methane_n;
    
    
    disp(sim)
    toc
end
%% Parsing outcome
ethane_13Cclumplog=ethane_13C2_log.*ethane_13Cn_log./(ethane_13C1_log).^2*4;
if sims>1
    cnoutiso_log_D=cumsum(cnoutiso_D);
    cnoutiso_log_C=cumsum(cnoutiso_C);
else
    cnoutiso_log_D=cnoutiso_D;
    cnoutiso_log_C=cnoutiso_C;
end

for sim=1:sims
    for cn=1:6
        for pos=1:cn+2  %1-21-2021 added for new isotopomer collection
            if pos<(cn+3)/2
                cnoutiso_D(sim,cn, pos)=cnoutiso_D(sim,cn, pos)+cnoutiso_D(sim,cn, cn+3-pos);
                cnoutiso_C(sim,cn, pos)=cnoutiso_C(sim,cn, pos)+cnoutiso_C(sim,cn, cn+3-pos);
                cnoutiso_log_D(sim,cn, pos)=cnoutiso_log_D(sim,cn, pos)+cnoutiso_log_D(sim,cn, cn+3-pos);
                cnoutiso_log_C(sim,cn, pos)=cnoutiso_log_C(sim,cn, pos)+cnoutiso_log_C(sim,cn, cn+3-pos);
            elseif pos>(cn+3)/2
                cnoutiso_D(sim,cn, pos)=cnoutiso_D(sim,cn, cn+3-pos);
                cnoutiso_C(sim,cn, pos)=cnoutiso_C(sim,cn, cn+3-pos);
                cnoutiso_log_D(sim,cn, pos)=cnoutiso_log_D(sim,cn, cn+3-pos);
                cnoutiso_log_C(sim,cn, pos)=cnoutiso_log_C(sim,cn, cn+3-pos);
            end
        end
        %         if mod(cn+2,2)==1
        %             midpos=round((cn+1+2)/2);
        %             cnoutiso_D(sim,cn, midpos)=cnoutiso_D(sim,cn, midpos)/2;
        %             cnoutiso_C(sim,cn, midpos)=cnoutiso_C(sim,cn, midpos)/2;
        %             cnoutiso_log_D(sim,cn, midpos)=cnoutiso_log_D(sim,cn, midpos)/2;
        %             cnoutiso_log_C(sim,cn, midpos)=cnoutiso_log_C(sim,cn, midpos)/2;
        %         end
    end
end

D_clean=cnoutiso_log_D;
C_clean=cnoutiso_log_C;

for sim=1:sims
    for cn=1:6
        for pos=1:cn+2
            if or(pos==1,pos==cn+2)
                D_clean(sim,cn,pos)=D_clean(sim,cn,pos)/3;
            else
                D_clean(sim,cn,pos)=D_clean(sim,cn,pos)/2;
            end
        end
        if mod(cn+2,2)==1
            D_clean(sim,cn,round((cn+1+2)/2))=D_clean(sim,cn,round((cn+1+2)/2))*2;
            C_clean(sim,cn,round((cn+1+2)/2))=C_clean(sim,cn,round((cn+1+2)/2))*2;
        end
    end
end

save(strcat('2022results\','result_Beta_BSKIEonly_',datestr(date,'yyyy-mm-dd'),'_',Filename,'_',num2str(stopat*100),'_',num2str(ThermalProgram(2,1))));
beep



