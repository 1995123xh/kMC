function [rate_mega, singlebonds]=getRateMatrix_bondM(ConnM, BondM, Dmarker, Cmarker, atoms, Temp)

[Dgraph,Dindex]=find(Dmarker==1);
[Cgraph,Cindex]=find(Cmarker==1);
[predic, KIE, KIE_13C] =kinPara(Temp);

deg=sum(BondM,2);
deg_mega=deg';
deg_mega(deg_mega>4)=4;
% BoM_mega={};
% BoM_mega{1}=adjacency(BoG_mega{1},'weighted');    %note: when this is generated

sgbidxM=zeros(size(ConnM));
[beginN, numc]=find(BondM==1);
endN=beginN;
for i=1:length(endN)
    endN(i)=ConnM(beginN(i), numc(i));
    sgbidxM(beginN(i), numc(i))=i;
end

% NodesList=BoG_mega{1}.Edges.EndNodes;
% BondWeights=BoG_mega{1}.Edges.Weight;
BondTypes=atoms([beginN,endN]);
BondDegs=deg_mega([beginN]).^2+deg_mega([endN]).^2;  % square sum of carbon degrees of a bond, used for calculating CC rates
BondDegs(BondDegs==2)=1;
BondDegs(BondDegs==5)=2;
BondDegs(BondDegs==10)=3;
BondDegs(BondDegs==17)=4;
BondDegs(BondDegs==8)=5;
BondDegs(BondDegs==13)=6;
BondDegs(BondDegs==20)=7;
BondDegs(BondDegs==18)=8;
BondDegs(BondDegs==25)=9;
BondDegs(BondDegs==32)=10;

AsN=endN(BondTypes(:,1)=='A');  %aromatic
atoms(AsN)='R';   %aromatic sidechain

intBond=zeros(size(BondTypes));
intBond(BondTypes=='C')=1;
intBond(BondTypes=='A')=2;
intBond(BondTypes=='N')=3;
intBond(BondTypes=='O')=4;
intBond(BondTypes=='S')=5;
% singlebonds=find(BoG_mega{1}.Edges.Weight==1); %Find all single bonds in the graph
singlebonds=prod(intBond,2);
singlebonds(singlebonds<1)=13;
singlebonds(singlebonds>5)=13;

rate=zeros(length(singlebonds),1);

rate(singlebonds==1)=predic.CC(BondDegs(singlebonds==1));
rate(singlebonds==2)=predic.CA;
rate(singlebonds==3)=predic.CN;
rate(singlebonds==4)=predic.CO;
rate(singlebonds==5)=predic.CS;

BondTypes=atoms([beginN,endN]);  %now including the aromatic side chains
intBond(BondTypes=='C')=1;
intBond(BondTypes=='R')=6;
singlebondsAs=prod(intBond,2);
rate(singlebondsAs==6)=predic.CAs;


rate_mega=rate';

if Dindex
    for i=1:length(Dgraph)
        g=Dgraph(i);
        indx=Dindex(i);
        Dnodes=ConnM(indx, BondM(indx,:)==1);
        Dedge=sgbidxM(indx, BondM(indx,:)==1);
        for b=1:length(Dnodes)
                singlebondindx = Dedge(b);
                bondtype=singlebonds(singlebondindx);
                if bondtype~=13
                    if bondtype<3
                        rate_mega(g,singlebondindx)=rate_mega(g,singlebondindx)/KIE(bondtype, deg_mega(g,indx));
                    elseif bondtype>3
                        rate_mega(g,singlebondindx)=rate_mega(g,singlebondindx)/KIE(bondtype-1, deg_mega(g,indx));
                    end
                end
                        secDnodes = ConnM(Dnodes(b), BondM(Dnodes(b),:)==1); %add secondary KIEs
                        secDedge = sgbidxM(Dnodes(b), BondM(Dnodes(b),:)==1);
                        secDnodes = secDnodes(secDnodes~=indx);
                        secDedge = secDedge(secDnodes~=indx);
                        for bs=1:length(secDedge)
                                singlebondindx = secDedge(bs);
                                rate_mega(g,singlebondindx) = rate_mega(g,singlebondindx)/KIE(5, deg_mega(g, indx));
                        end
        end
    end
end
%13C KIE
if Cindex
    for i=1:length(Cgraph)
        g=Cgraph(i);
        indx=Cindex(i);
        Cnodes=ConnM(indx, BondM(indx,:)==1);
        Cedge=sgbidxM(indx, BondM(indx,:)==1);
        for b=1:length(Cnodes)
                singlebondindx = Cedge(b);
                %                 singlebondindx = sgbmatrix(indx, Cnodes(b));
                bondtype=singlebonds(singlebondindx);
                if bondtype~=13
                    if bondtype<3
                        rate_mega(g,singlebondindx)=rate_mega(g,singlebondindx)/KIE_13C(bondtype, deg_mega(g,indx));
                    elseif bondtype>3
                        rate_mega(g,singlebondindx)=rate_mega(g,singlebondindx)/KIE_13C(bondtype-1, deg_mega(g,indx));
                    end
                end
        end
    end
end

singlebonds=[singlebonds, beginN, endN];
end


