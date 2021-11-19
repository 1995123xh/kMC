function rate_mega=getRateMatrix(BoG_mega, BoM_mega, singlebonds, sgbmatrix, deg_mega, Dmarker, Cmarker, Temp)

[Dgraph,Dindex]=find(Dmarker==1);
[Cgraph,Cindex]=find(Cmarker==1);
[predic, KIE, KIE_13C] =kinPara(Temp);

sgbatoms=singlebonds(:,3:4); %beginN and endN of singlebonds
BondDegs=zeros(size(deg_mega,1),size(singlebonds,1));
for i=1:size(deg_mega,1)
    BondDegs(i,:)=deg_mega(i,sgbatoms(:,1)).^2+deg_mega(i,sgbatoms(:,2)).^2;
end
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
BondDegs(BondDegs==0)=11;   %not actually a bond

predic.CC(11)=0;

rate_mega=zeros(size(BoM_mega,2),size(singlebonds,1));
for i=1:size(rate_mega,1)
    rate=zeros(size(singlebonds,1),1);
    BondDegsSlice=BondDegs(i,:);
    rate(singlebonds(:,2)==1)=predic.CC(BondDegsSlice(singlebonds(:,2)==1));
    rate(singlebonds(:,2)==2)=predic.CA;
    rate(singlebonds(:,2)==3)=predic.CN;
    rate(singlebonds(:,2)==4)=predic.CO;
    rate(singlebonds(:,2)==5)=predic.CS;
    rate(singlebonds(:,2)==6)=predic.CAs;
    rate_mega(i,:)=rate;
end

if Dindex
    for i=1:length(Dgraph)
        g=Dgraph(i);
        indx=Dindex(i);
        [Dedge, Dnodes]=outedges(BoG_mega{g}, Dindex(i));

        for b=1:length(Dnodes)
            if BoM_mega{g}(indx, Dnodes(b))==1
                %                 singlebondedge = findedge(BoG, indx, Dnodes(b));      %needs to relocate where the edge is because some bonds have been broken;
                %                 singlebondindx = find(singlebonds==singlebondedge);
                singlebondindx = sgbmatrix(indx, Dnodes(b));
                bondtype=singlebonds(singlebondindx,2);
                if bondtype~=13
                    if bondtype<3
                        rate_mega(g,singlebondindx)=rate_mega(g,singlebondindx)/KIE(bondtype, deg_mega(g,indx));
                    elseif bondtype>3 && bondtype<6
                        rate_mega(g,singlebondindx)=rate_mega(g,singlebondindx)/KIE(bondtype-1, deg_mega(g,indx));
                    elseif bondtype==6
                        rate_mega(g,singlebondindx)=rate_mega(g,singlebondindx)/KIE(1, deg_mega(g,indx)); %C-R bond; nodata
                    end
                end
                %                 switch strcat(atoms(indx),atoms(Dnodes(b)))
                %                     case 'CC'
                %                         rate_mega(g,singlebondindx)=rate_mega(g,singlebondindx)/KIE(1, deg_mega(g,indx));      %kinetic isotope effect based on degrees
                %                     case {'CA', 'AC'}
                %                         rate_mega(g,singlebondindx)=rate_mega(g,singlebondindx)/KIE(2, deg_mega(g,indx));
                %                     case {'CO', 'OC'}
                %                         rate_mega(g,singlebondindx)=rate_mega(g,singlebondindx)/KIE(3, deg_mega(g,indx));
                %                     case {'CS', 'SC'}
                %                         rate_mega(g,singlebondindx)=rate_mega(g,singlebondindx)/KIE(4, deg_mega(g,indx));
                %                 end
            end
            [secDedge,secDnodes] = outedges(BoG_mega{g}, Dnodes(b));   %add secondary KIEs
            secDnodes = secDnodes(secDedge~=Dedge(b));
            for bs=1:length(secDnodes)
                if BoM_mega{g}(Dnodes(b),secDnodes(bs))==1
                    %                     singlebondedge = findedge(BoG, Dnodes(b), secDnodes(bs));
                    %                     singlebondindx = find(singlebonds==singlebondedge);
                    singlebondindx = sgbmatrix(Dnodes(b), secDnodes(bs));
                    rate_mega(g,singlebondindx) = rate_mega(g,singlebondindx)/KIE(5, deg_mega(g, indx));
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
        [Cedge, Cnodes]=outedges(BoG_mega{g}, Cindex(i));
        for b=1:length(Cnodes)
            if BoM_mega{g}(indx, Cnodes(b))==1
                %                 singlebondedge = findedge(BoG, indx, Cnodes(b));
                %                 singlebondindx = find(singlebonds==singlebondedge);
                singlebondindx = sgbmatrix(indx, Cnodes(b));
                bondtype=singlebonds(singlebondindx,2);
                if bondtype~=13
                    if bondtype<3
                        rate_mega(g,singlebondindx)=rate_mega(g,singlebondindx)/KIE_13C(bondtype, deg_mega(g,indx));
                    elseif bondtype>3 && bondtype<6
                        rate_mega(g,singlebondindx)=rate_mega(g,singlebondindx)/KIE_13C(bondtype-1, deg_mega(g,indx));
                    elseif bondtype==6
                        rate_mega(g,singlebondindx)=rate_mega(g,singlebondindx)/KIE_13C(6, deg_mega(g,indx));
                    end
                end
                %                 switch strcat(atoms(indx),atoms(Cnodes(b)))
                %                     case 'CC'
                %                         rate_mega(g,singlebondindx)=rate_mega(g,singlebondindx)/KIE_13C(1, deg_mega(g,indx));      %kinetic isotope effect based on degrees
                %                     case {'CA', 'AC'}
                %                         rate_mega(g,singlebondindx)=rate_mega(g,singlebondindx)/KIE_13C(1, deg_mega(g,indx));
                %                     case {'CO', 'OC'}
                %                         rate_mega(g,singlebondindx)=rate_mega(g,singlebondindx)/KIE_13C(1, deg_mega(g,indx));
                %                     case {'CS', 'SC'}
                %                         rate_mega(g,singlebondindx)=rate_mega(g,singlebondindx)/KIE_13C(1, deg_mega(g,indx));
                %                 end
            end
            [secCedge,secCnodes] = outedges(BoG_mega{g}, Cnodes(b));
            secCnodes = secCnodes(secCedge~=Cedge(b));
            for bs=1:length(secCnodes)
                if BoM_mega{g}(Cnodes(b),secCnodes(bs))==1
                    %                     singlebondedge = findedge(BoG, Cnodes(b), secCnodes(bs));
                    %                     singlebondindx = find(singlebonds==singlebondedge);
                    singlebondindx = sgbmatrix(Cnodes(b), secCnodes(bs));
                    rate_mega(g,singlebondindx) = rate_mega(g,singlebondindx)/KIE_13C(5, deg_mega(g, indx));
                end
            end
        end
    end
end
end


