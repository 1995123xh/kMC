%% Making a kinetic MonteCarlo for simulating propane generation from
%%% kerogen molecules
%%%
clear all
Filename='kIA';

disp('loading...')
load(Filename)

nmol=1;         %number of molecules in the graph. 1000-2000 atoms have highest efficiency
sz=100;         %number of graphs in the cluster
enrichfactor=50; %the factor of deuterium concentration
enrichfactor_13C=5;
crkratio=1;    %fraction of bonds going broken
stopat=20;      %percentage of completion
stopat=stopat/100;
sims=10;        %number of MC simulations
rep=12;       %number of reps in the parfor loop
updatestepsize=50;  %step size for refreshing the reaction rates due to temeprature change and increasing capping D. If you don't wish to temperature(e.g. for computing efficiency), use a very large number here)

ThermalProgram=[0 0.2 0.5 0.7 1; [500 500 500 500 500]];  % A thermal program where first row is fraction to completion and second row is temeprature
inittemp=ThermalProgram(2,1);


BoM_initial=BoM;    %transfering the original molecule somewhere else...
atoms_initial=atoms;

graphstable=cell(sims,rep);
gcount=0;

BoM=zeros(nmol*length(BoM_initial),nmol*length(BoM_initial));
for i=1:nmol
    BoM(((i-1)*length(BoM_initial)+1):(i)*length(BoM_initial),((i-1)*length(BoM_initial)+1):(i)*length(BoM_initial))=BoM_initial;
    if i>1
        atoms=strcat(atoms, atoms_initial);
    end
end

BoG=graph(BoM, 'upper','omitselfloops');    %Construct a graph of bond orders
BoG0=BoG;           %the original molecule


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
deg_mega=zeros(sz, length(deg));
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
R0_CAP=10*155.76e-6*(1-20/1000)*enrichfactor;
R0_13C=0.01118*enrichfactor_13C;
F0=R0/(R0+1);
F0_13C=R0_13C/(R0_13C+1);
% T=273.15+180; %kelvin
R=8.314;
aH=3;        %fractionation factor for capping hydrogen

disp('Initializing...')

%% Mark the aromaticity
arobonds=find(abs(BoG.Edges.Weight-1.5)<0.001);
for b=arobonds
    [beginA, endA] = findedge(BoG, b);
    atoms(beginA)='A';
    atoms(endA)='A';
end

%% creat singlebonds and singlebondsmatrix
NodesList=BoG.Edges.EndNodes;
issinglebonds=BoG.Edges.Weight==1;
sgbNodesList=NodesList(issinglebonds,:);
BondTypes=atoms(sgbNodesList);
singlebonds=find(BoG.Edges.Weight==1); %Find all single bonds in the graph
singlebonds=[singlebonds,singlebonds];
AsN=[];
for i=1:length(singlebonds)
    if isequal(BondTypes(i,:),'CC')
        singlebonds(i,2)=1;
    elseif isequal(BondTypes(i,:),'CA')
        singlebonds(i,2)=2;
        AsN=[AsN, sgbNodesList(i,1)];
    elseif isequal(BondTypes(i,:), 'AC')
        singlebonds(i,2)=2;
        AsN=[AsN, sgbNodesList(i,2)];
    elseif isequal(BondTypes(i,:),'CN') || isequal(BondTypes(i,:),'NC')
        singlebonds(i,2)=3;
    elseif isequal(BondTypes(i,:),'CO') || isequal(BondTypes(i,:),'OC')
        singlebonds(i,2)=4;
    elseif isequal(BondTypes(i,:),'CS') || isequal(BondTypes(i,:),'SC')
        singlebonds(i,2)=5;
    else
        singlebonds(i,2)=13; %not included bonds;
    end
end
singlebonds=[singlebonds, sgbNodesList];  %four columns: 1 sgb index 2 bond type 3&4 sgb atoms
atoms(AsN)='R';
BondTypes=atoms(NodesList(issinglebonds,:));
for i=1:size(singlebonds,1)
    if isequal(BondTypes(i,:),'CR') || isequal(BondTypes(i,:),'RC')
        singlebonds(i,2)=6;
    end
end
atoms(AsN)='C';
singlebonds0=singlebonds;
sgbmatrix0=BoM;
for i = 1:size(singlebonds,1)
    [beginA, endA] = findedge(BoG, singlebonds(i,1));
    sgbmatrix0(beginA,endA) = i;
    sgbmatrix0(endA,beginA) = i;
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

hnumber=4-deg_mega0;
dabundance=hnumber.*F0.*(1-F0).^(hnumber-1); % fixed a bug of no exponent
% parpool('local')
tic
% BoG_mega0=parallel.pool.Constant(BoG_mega0);
disp('Starts Monte-Carlo...');
timeM=zeros(sims, rep);
methaneN=[];
for sim=1:sims
    cnoutisotemp_D=cell(rep,1);   %pre-allocate these counters
    cnoutisotemp_C=cell(rep,1);
%      mypool=gcp;
%      addAttachedFiles(mypool,{'getRateMatrix.m','kinPara.m'})
%     ticBytes(gcp);


    parfor iteration=1:rep
         rng('shuffle');
        %data transfer
        BoG_mega=cell(1,sz);
        BoG_mega(:)={BoG};
%         BoM_mega=full(BoM_mega0);
        BoM_mega=cell(1,sz);
        BoM_mega(:)={BoM};
        deg_mega=deg_mega0;
        atomorder=atoms;
        degorder=deg;
%        BoMorder=BoM;
        ThermalProgramPar=ThermalProgram;
%         cutedgelist=[];
         singlebonds=singlebonds0;
         sgbmatrix=sgbmatrix0;
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
        %calculating rate
        [Dgraph,Dindex]=find(Dmarker==1);  %the graph and the index of the deuterium
        [Cgraph,Cindex]=find(Cmarker==1);  %the graph and the index of the deuterium
        [predic, KIE, KIE_13C]=kinPara(inittemp);
        rate_mega=getRateMatrix(BoG_mega, BoM_mega, singlebonds, sgbmatrix, deg_mega, Dmarker, Cmarker, inittemp);
        totalbonds=round(crkratio*sum(sum(rate_mega~=0)));
        
        for step=1:round(totalbonds*stopat) % core loop for breaking bonds
            T=interp1(ThermalProgramPar(1,:), ThermalProgramPar(2,:), step/totalbonds);
            if mod(step, updatestepsize)==0
 %               [predic, KIE, KIE_13C]=kinPara(T);  %update kinetics based on temperature
                new_rate_mega=getRateMatrix(BoG_mega, BoM_mega, singlebonds, sgbmatrix, deg_mega, Dmarker, Cmarker, T);
                rate_mega(rate_mega~=0)=new_rate_mega(rate_mega~=0);
            end
            sumrate=sum(sum(rate_mega));
            cumrate=cumsum(rate_mega(:)); % Generate a cumulative sum of rates
            TimeRand=rand();
            %timeM(sim, iteration)=timeM(sim,iteration)+1/sumrate*log(1/TimeRand);
            
            dice=rand()*sumrate;
            
            i = find(cumrate > dice,1);
            index=ceil(i/sz); %the graph; this is because the 1-D ordering of matrix goes by columns
            g=i-(index-1)*sz;
            
            %cut the bond
            [A1, A2] = findedge(BoG,singlebonds(index,1));
            %             if Dmarker(g, A1) || Dmarker(g, A2) || Cmarker(g, A1) || Cmarker(g, A2)
            %                 disp(step)
            %             end
            %             if or(ismember(A1, prop), ismember(A2, prop))
            %                 break
            %             end
            BoM_mega{g}(A1, A2)=0;
            BoM_mega{g}(A2, A1)=0;
            deg_mega(g,A1)=deg_mega(g,A1)-1;
            deg_mega(g,A2)=deg_mega(g,A2)-1;

            rate_mega(g, index)=0;
            BoG_mega{g}=rmedge(BoG_mega{g}, A1,A2);
            %add capping hydrogen/deterium
            
            for A_13C=[A1,A2] %add capping hydrogen
                dice=rand();
                if dice<R0_CAP/aH/(R0_CAP/aH+1)
                    Dmarker(g, A_13C)=1;           
                end
            end
        end
        
        %Count hydrocarbons
        %         tables=cell(sz,1);
        for g=1:sz
            %             tables{g}=BoG_mega{g}.Edges;
            splitg=conncomp(BoG_mega{g},'OutputForm','cell'); %decompose eachgraph into sepearte parts
            for i=1:length(splitg)
                if length(splitg{i})<9          %ignore larger fragments
                    splitg_i=subgraph(BoG_mega{g},splitg{i});   %construct a seperate graph
                    if length(splitg{i})==2  %Find the ethanes
                        pnodes=splitg{i};
                        if atomorder(pnodes)=='CC'
                            if sgbmatrix(pnodes(1), pnodes(2))~=0 %single C-C bonds
                                Dlabel=[(Dmarker(g, pnodes(1))), (Dmarker(g, pnodes(2)))];
                                if or(isequal(Dlabel, [1, 0]), isequal(Dlabel,[0, 1]))
                                    ethane_D1=ethane_D1+1;
                                elseif isequal(Dlabel, [1, 1])
                                    ethane_D2=ethane_D2+1;
                                else
                                    ethane_Dn=ethane_Dn+1;
                                end
                                Clabel=[Cmarker(g, pnodes(1)), Cmarker(g, pnodes(2))];
                                if or(isequal(Clabel, [1, 0]), isequal(Clabel,[0, 1]))
                                    ethane_13C1=ethane_13C1+1;
                                elseif isequal(Clabel, [1, 1])
                                    ethane_13C2=ethane_13C2+1;
                                else
                                    ethane_13Cn=ethane_13Cn+1;
                                end
                            end
                        end
                    elseif length(splitg{i})==1  %Find the methanes
                        pnodes=splitg{i};
                        if atomorder(pnodes)=='C'
                            methaneN=[methaneN;pnodes];
                            Dlabel=Dmarker(g, pnodes);
                            Clabel=Cmarker(g, pnodes);
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
                        cn=length(splitg{i})-2;   %
                        pnodes=splitg{i};                           %
                        if all(atomorder(pnodes)=='C')                   %no hetero atoms
                            Dlabel=[Dmarker(g, pnodes)];
                            Clabel=[Cmarker(g, pnodes)];
                            splitg_i.Nodes.DLabel=double(Dlabel');
                            splitg_i.Nodes.CLabel=double(Clabel');
                            if isisomorphic(cnoutorder(cn).BoG,splitg_i, 'EdgeVariables','Weight')
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

save(strcat('result_',datestr(date,'yyyy-mm-dd'),'_',Filename,'_',num2str(stopat*100)))

%% Calculate kinetic parameters based on RMG and KIEs based on Tang 2005







