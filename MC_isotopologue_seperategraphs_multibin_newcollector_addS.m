%%% Making a kinetic MonteCarlo for simulating propane generation from
%%% kerogen molecules
%%%
clear all
Filename='kIA.mat'

nmol=1;         %number of molecules in the graph. 1000-2000 atoms have highest efficiency
sz=500;         %number of graphs in the cluster
enrichfactor=50; %the factor of deuterium concentration
enrichfactor_13C=5;
enrichfactor_34S=1;
crkratio=0.2;    %fraction of bonds going broken
sims=50;        %number of MC simulations
rep=168;       %number of reps in the parfor loop

disp('loading...')
load(Filename)
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


BoG_mega=cell(1,sz);   %cell array of graphs
BoG_mega(:)={BoG};
BoG_mega0=BoG_mega;

BoM_mega=cell(1,sz);
BoM_mega(:)={BoM};
BoM0=BoM;
BoM_mega0=BoM_mega;

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

MeSH_D=0;
MeSH_13C=0;
MeSH_34S=0;
MeSH_all=0;

H2S_34S=0;
H2S_n=0;

DMS_34S=0;
DMS_n=0;

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
R0=155.76e-6*enrichfactor;
R0_13C=0.01118*enrichfactor_13C;
R0_34S=0.0441626*enrichfactor_34S;   %VCDT
F0=R0/(R0+1);
F0_13C=R0_13C/(R0_13C+1);
F0_34S=R0_34S/(R0_34S+1);
T=273.15+180; %kelvin
R=8.314;
RH2=(1-836/1000)*R0;             %hydrogen isotopes of the hydrogen gas
%% Calculating kinetic parameters based on and KIEs based on Tang 2005
disp('Initializing...')
predic.CC=6.64671E23*T^-1.99*exp(-375755/R/T);       %rates in s^-1 RMG-Py rate rules
predic.CA=1.58849e20*(T)^-0.76*exp(-442154/8.314/T); %
predic.CN=0; %4.68781e8*T^0.52*exp(-349516.61/R/T);
predic.CO=2.63224e18*T^-0.64*exp(-351864.4/R/T);
predic.CS=4.3545e13*T^0.56*exp(-305100.73/R/T);

DDH=zeros(4,4);
DDH(1,1)=384; % 1 deg carbon
DDH(1,2)=327; % 2 deg carbon
DDH(1,3)=316; % 3 deg carbon

DDH(2,1)=359; % 1 deg C-A
DDH(2,2)=359+(DDH(1,2)-DDH(1,1));
DDH(2,3)=359+(DDH(1,3)-DDH(1,1));

DDH(3,1)=427; % 1 deg C-A
DDH(3,2)=427+(DDH(1,2)-DDH(1,1));
DDH(3,3)=427+(DDH(1,3)-DDH(1,1));

DDH(4,1)=377; % 1 deg C-A
DDH(4,2)=377+(DDH(1,2)-DDH(1,1));
DDH(4,3)=377+(DDH(1,3)-DDH(1,1));
DDH(:,4)=DDH(:,3);
% DDH(4)=359+(DDH(2)-DDH(1)); % 2 deg C-A           %estimated by C-C scaling
% DDH(5)=359+(DDH(3)-DDH(1)); % 3 deg C-A
% DDH(6)=427+(DDH(2)-DDH(1)); % 2 deg C-O
% DDH(7)=427+(DDH(3)-DDH(1)); % 3 deg C-O
% DDH(8)=377+(DDH(2)-DDH(1)); % 2 deg C-S
% DDH(9)=377+(DDH(3)-DDH(1)); % 2 deg C-S
A_H=1.07;
KIE=1./(A_H*exp(-DDH*4.184/8.314/T));
aH=3;
% KIE=[1,1,1,1,1,1,1,1,1];
% aH=1;

% For carbon isotopes from table 2 in Tang et al., 2001
DDC=zeros(4,4);  %four for the 4 types of bonds, and 3 for carbon degrees
A_13C=zeros(4,4);

DDC(1,1)=48.853; % DEA of 1 deg carbon, in cal/mol
A_13C(1,1)=1.03;
DDC(1,2)=38.763; % 2 deg carbon
A_13C(1,2)=1.02;
DDC(1,3)=35.261; % 3 deg carbon
A_13C(1,3)=1.023;

DDC(2,1)=51.562; % 1 deg C-A
A_13C(2,1)=1.034;
DDC(2,2)=28.314; % 2 deg C-A
A_13C(2,2)=1.018;
DDC(2,3)=(28.314+35.261-38.763); % 3 deg C-A
A_13C(2,3)=1.018;

DDC(3,1)=51.539; % 1 deg C-O
A_13C(3,1)=1.034;
DDC(3,2)=38.595; % 2 deg C-O
A_13C(3,2)=1.025;
DDC(3,3)=(38.595+35.261-38.763); % 3 deg C-O
A_13C(3,3)=1.025;

DDC(4,1)=34.962; % 1 deg C-S
A_13C(4,1)=1.026; % 1 deg C-S
DDC(4,2)=36.325; % 2 deg C-S
A_13C(4,2)=1.022; % 2 deg C-S
DDC(4,3)=(36.325+35.261-38.763); % 3 deg C-S
A_13C(4,3)=1.022; % 2 deg C-S

DDC(:,4)=DDC(:,3);
A_13C(:,4)=A_13C(:,3);
KIE_13C=1./(A_13C.*exp(-DDC.*4.184./8.314./T));

% for sulfur isotopes from Tang et al., 2000
DDS=zeros(1,2);
A_34S=zeros(1,2);
DDS(1,1)= 11.14;  % 1deg S-C
DDS(1,2)=9.349;   % 2deg S-C
A_34S(1,1)=1.007;
A_34S(1,2)=1.008;

KIE_34S=1./(A_34S.*exp(-DDS.*4.184./8.314./T));

%% Mark the aromaticity
arobonds=find(abs(BoG.Edges.Weight-1.5)<0.001);
for b=arobonds
    [beginA, endA] = findedge(BoG, b);
    atoms(beginA)='A';
    atoms(endA)='A';
end

%% Determine rate of each cleavage
singlebonds=find(BoG.Edges.Weight==1); %Find all single bonds in the graph

for i = 1:length(singlebonds)          %%Determine the type of bond that could be broken
    [beginA, endA] = findedge(BoG,singlebonds(i));
    if atoms(beginA) == 'C' && atoms(endA) == 'C'
        rate(i)=predic.CC;
    elseif or(strcat(atoms(beginA),atoms(endA))=='CA' , strcat(atoms(beginA),atoms(endA)) == 'AC')
        rate(i)=predic.CA;
    elseif or(strcat(atoms(beginA),atoms(endA))=='CN' , strcat(atoms(beginA),atoms(endA)) == 'NC')
        rate(i)=predic.CN;
    elseif or(strcat(atoms(beginA),atoms(endA))=='CO' , strcat(atoms(beginA),atoms(endA)) == 'OC')
        rate(i)=predic.CO;
    elseif or(strcat(atoms(beginA),atoms(endA))=='CS' , strcat(atoms(beginA),atoms(endA)) == 'SC')
        rate(i)=predic.CS;
    end
end
rate0=rate;

rate_mega=zeros(sz,length(rate));
for i=1:sz
    rate_mega(i,:)=rate;
end

rate_mega0=rate_mega;
sumrate0=sum(sum(rate_mega0));
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

tic
disp('Starts Monte-Carlo...')
for sim=1:sims
    for iteration=1:rep
        BoG_mega=BoG_mega0;
        BoM_mega=BoM_mega0;
        deg_mega=deg_mega0;
        rate_mega=rate_mega0;
        atomorder=atoms;
        degorder=deg;
        BoMorder=BoM;
        KIEpar=KIE;        %hydrogen KIE
        KIEpar_13C=KIE_13C;
        KIEpar_34S=KIE_34S;
        cutedgelist=[];
        b=0;
        cnoutiso_Dpar=zeros(length(cnout),8+1);
        cnoutiso_Cpar=zeros(length(cnout),8+1);
        cnoutorder=cnout;
        %generating random deuterium distribution
        cabundance=F0_13C; %
        randomarray=rand(size(hnumber));
        Dmarker=randomarray<dabundance;
        randomarray=rand(size(hnumber));
        Dmarker(:,atoms~='C')=0;
        Cmarker=randomarray<F0_13C;
        Cmarker(:,atoms~='C')=0;                   % 2019-06-15 fixed the error of wrong indexing here
        randomarray=rand(size(hnumber));
        Smarker=randomarray<F0_34S;
        Smarker(:,atoms~='S')=0;                   % 2019-06-15 fixed the error of wrong indexing here
        %calculating rate
        [Dgraph,Dindex]=find(Dmarker==1);  %the graph and the index of the deuterium
        [Cgraph,Cindex]=find(Cmarker==1);  %the graph and the index of the deuterium
        [Sgraph,Sindex]=find(Smarker==1);  %the graph and the index of the deuterium
        if Dindex
            for i=1:length(Dgraph)
                g=Dgraph(i);
                indx=Dindex(i);
                Dedge=outedge(BoG_mega{g}, Dindex(i));
                for b=Dedge'
                    if BoG_mega{g}.Edges.Weight(b)==1
                        singlebondindx=find(singlebonds==b);
                        [beginA, endA] = findedge(BoG,b);
                        switch strcat(atomorder(beginA),atomorder(endA))
                            case 'CC'
                                rate_mega(g,singlebondindx)=rate_mega(g,singlebondindx)/KIEpar(1, degorder(indx));      %kinetic isotope effect based on degrees
                            case {'CA', 'AC'}
                                rate_mega(g,singlebondindx)=rate_mega(g,singlebondindx)/KIEpar(2, degorder(indx));
                            case {'CO', 'OC'}
                                rate_mega(g,singlebondindx)=rate_mega(g,singlebondindx)/KIEpar(3, degorder(indx));
                            case {'CS', 'SC'}
                                rate_mega(g,singlebondindx)=rate_mega(g,singlebondindx)/KIEpar(4, degorder(indx));
                        end
                    end
                end
            end
        end
        
        %13C KIE
        if Cindex
            for i=1:length(Cgraph)
                g=Cgraph(i);
                indx=Cindex(i);
                Cedge=outedge(BoG_mega{g}, Cindex(i));
                for b=Cedge'
                    if BoG_mega{g}.Edges.Weight(b)==1
                        singlebondindx=find(singlebonds==b);
                        [beginA, endA] = findedge(BoG,b);
                        switch strcat(atomorder(beginA),atomorder(endA))
                            case 'CC'
                                rate_mega(g,singlebondindx)=rate_mega(g,singlebondindx)/KIEpar_13C(1, degorder(indx));      %kinetic isotope effect based on degrees
                            case {'CA', 'AC'}
                                rate_mega(g,singlebondindx)=rate_mega(g,singlebondindx)/KIEpar_13C(1, degorder(indx));
                            case {'CO', 'OC'}
                                rate_mega(g,singlebondindx)=rate_mega(g,singlebondindx)/KIEpar_13C(1, degorder(indx));
                            case {'CS', 'SC'}
                                rate_mega(g,singlebondindx)=rate_mega(g,singlebondindx)/KIEpar_13C(1, degorder(indx));
                        end
                    end
                end
            end
        end
        if Sindex
            for i=1:length(Sgraph)
                g=Sgraph(i);
                indx=Sindex(i);
                Sedge=outedge(BoG_mega{g}, Sindex(i));
                if length(Sedge)<3
                    for b=Sedge'
                        if BoG_mega{g}.Edges.Weight(b)==1
                            singlebondindx=find(singlebonds==b);
                            [beginA, endA] = findedge(BoG,b);
                            switch strcat(atomorder(beginA),atomorder(endA))
                                case {'CS', 'SC'}
                                    rate_mega(g,singlebondindx)=rate_mega(g,singlebondindx)/KIEpar_34S(1, degorder(indx));
                            end
                        end
                    end
                end
            end
        end
        
        for step=1:round(crkratio*numel(rate_mega))
            sumrate=sum(sum(rate_mega));
            cumrate=cumsum(rate_mega(:)); % Generate a cumulative sum of rates
            dice=rand()*sumrate;
            i = find(cumrate > dice,1);
            index=ceil(i/sz); %the graph
            g=i-(index-1)*sz;
            
            %cut the bond
            [A1, A2] = findedge(BoG0,singlebonds(index));
            %             if or(ismember(A1, prop), ismember(A2, prop))
            %                 break
            %             end
            BoM_mega{g}(A1, A2)=0;
            BoM_mega{g}(A2, A1)=0;
            deg_mega(g,A1)=deg_mega(g,A1)-1;
            deg_mega(g,A2)=deg_mega(g,A2)-1;
            %     cleave(A1)=cleave(A1)+1;   %mark the cleavage
            %     cleave(A2)=cleave(A2)+1;
            rate_mega(g, index)=0;
            BoG_mega{g}=rmedge(BoG_mega{g}, A1,A2);
            
            %add capping hydrogen/deterium
            
            for A_13C=[A1,A2]
                dice=rand();
                if dice<R0/aH/(R0/aH+1)
                    Dmarker(g, A_13C)=not(Dmarker(g,A_13C));           %add capping hydrogen
                    if Dmarker(g, A_13C)==1
                        Dedge=outedge(BoG_mega{g}, A_13C);
                        for b=Dedge'
                            if BoG_mega{g}.Edges.Weight(b)==1
                                singlebondindx=find(singlebonds==b);
                                [beginA, endA] = findedge(BoG,b);
                                switch strcat(atomorder(beginA),atomorder(endA))
                                    case 'CC'
                                        rate_mega(g, singlebondindx)=rate_mega(g, singlebondindx)/KIEpar(1, deg_mega(g,A_13C));      %kinetic isotope effect based on degrees
                                    case {'CA', 'AC'}
                                        rate_mega(g, singlebondindx)=rate_mega(g, singlebondindx)/KIEpar(1, deg_mega(g,A_13C));
                                    case {'CO', 'OC'}
                                        rate_mega(g, singlebondindx)=rate_mega(g, singlebondindx)/KIEpar(1, deg_mega(g,A_13C));
                                    case {'CS', 'SC'}
                                        rate_mega(g, singlebondindx)=rate_mega(g, singlebondindx)/KIEpar(1, deg_mega(g,A_13C));
                                end
                            end
                        end
                    end
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
                            if BoM(pnodes(1), pnodes(2))<1.5 %single C-C bonds
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
                        elseif or(atomorder(pnodes)=='CS', atomorder(pnodes)=='SC')
                            if BoM(pnodes(1), pnodes(2))<1.5
                                if sum([Smarker(g, pnodes(1)), Smarker(g, pnodes(2))])==1
                                    MeSH_34S=MeSH_34S+1;
                                end
                                if sum([Cmarker(g, pnodes(1)), Cmarker(g, pnodes(2))])==1           %13C
                                    MeSH_13C=MeSH_13C+1;
                                end
                                if sum([Dmarker(g, pnodes(1)), Dmarker(g, pnodes(2))])==1           %13C
                                    MeSH_meD=MeSH_meD+1;
                                end
                                MeSH_all=MeSH_all+1;
                            end
                        end
                    elseif length(splitg{i})==1  %Find the methanes
                        pnodes=splitg{i};
                        if atomorder(pnodes)=='C'
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
                        elseif atomorder(pnodes)=='S'
                            Slabel=Smarker(g, pnodes);
                            if Slabel
                                H2S_34S=H2S_34S+1;
                            else
                                H2S_n=H2S_n+1;
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
                            if isisomorphic(cnoutorder(cn).BoG,splitg_i)
                                for isotopomer=1:size(cnoutorder(cn).label,1)
                                    if isisomorphic(cnoutorder(cn).labelgraph{isotopomer}, splitg_i, 'NodeVariables', 'DLabel')
                                        cnoutiso_Dpar(cn,isotopomer)=cnoutiso_Dpar(cn,isotopomer)+1;
                                    end
                                    if isisomorphic(cnoutorder(cn).labelgraph{isotopomer}, splitg_i, 'NodeVariables', 'CLabel')
                                        cnoutiso_Cpar(cn,isotopomer)=cnoutiso_Cpar(cn,isotopomer)+1;
                                    end
                                end
                            end
                        elseif sum(atomorder(pnodes)=='S')==1
                            if length(splitg{i})==3 && sum(atomorder(pnodes)=='C')==2
                                for p=pnodes
                                    if neighbors(BoG_mega{g},p)==2 && atomorder(p)=='S'   %find DMS
                                        if Slabel(p)
                                            DMS_34S=DMS_34S+1;
                                        else
                                            DMS_n=DMS_n+1;
                                        end
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
    
    cnoutiso_log_D=cumsum(cnoutiso_D);
    cnoutiso_log_C=cumsum(cnoutiso_C);
    disp(sim)
    toc
end
%% Parsing outcome
ethane_13Cclumplog=ethane_13C2_log.*ethane_13Cn_log./(ethane_13C1_log).^2*4;

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

save(strcat('result_',datestr(date,'yyyy-mm-dd'),'_',Filename))






