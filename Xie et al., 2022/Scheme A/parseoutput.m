%make a parsing script to deal with data output of the newcollector
%load(file)
close all

D_clean_discrete=cnoutiso_D;
C_clean_discrete=cnoutiso_C;

for sim=1:sims   %removing the symmetry factors
    for cn=1:6
        for pos=1:cn+2
            if or(pos==1,pos==cn+2)
                D_clean_discrete(sim,cn,pos)=D_clean_discrete(sim,cn,pos)/3;
            else
                D_clean_discrete(sim,cn,pos)=D_clean_discrete(sim,cn,pos)/2;
            end
        end
        if mod(cn+2,2)==1
            D_clean_discrete(sim,cn,round((cn+1+2)/2))=D_clean_discrete(sim,cn,round((cn+1+2)/2))*2;
            C_clean_discrete(sim,cn,round((cn+1+2)/2))=C_clean_discrete(sim,cn,round((cn+1+2)/2))*2;
        end
    end
end

prop_sp_C_discrete=C_clean_discrete(:,1,2)./C_clean_discrete(:,1,1)*1000-1000;
prop_sp_D_discrete=D_clean_discrete(:,1,2)./D_clean_discrete(:,1,1)*1000-1000;
Cerr=std(prop_sp_C_discrete)/sqrt(length(prop_sp_C_discrete))
Derr=std(prop_sp_D_discrete)/sqrt(length(prop_sp_D_discrete))


prop_sp_C_log=C_clean(:,1,2)./C_clean(:,1,1)*1000-1000;
prop_sp_D_log=D_clean(:,1,2)./D_clean(:,1,1)*1000-1000;

butane_sp_C_discrete=C_clean_discrete(:,2,2)./C_clean_discrete(:,2,1)*1000-1000;
butane_sp_D_discrete=D_clean_discrete(:,2,2)./D_clean_discrete(:,2,1)*1000-1000;
butane_sp_C_log=C_clean(:,2,2)./C_clean(:,2,1)*1000-1000;
butane_sp_D_log=D_clean(:,2,2)./D_clean(:,2,1)*1000-1000;
Cerr_butane=std(butane_sp_C_discrete)/sqrt(length(butane_sp_C_discrete));
Derr_butane=std(butane_sp_D_discrete)/sqrt(length(butane_sp_D_discrete));

R13C(1)=methane_13C/(methane_D+methane_n-methane_13CD);
RD(1)=methane_D/(methane_13C+methane_n-methane_13CD)/4;
R13C(2)=ethane_13C1/ethane_13Cn/2;
RD(2)=ethane_D1/ethane_Dn/6;

nalk=zeros(8,1);
nalk(1)=methane_n+methane_13C+methane_D-methane_13CD;
nalk(2)=ethane_Dn+ethane_D1;
for i=0:5
    nalk(i+3)=sum(cnoutiso_log_C(sims,i+1,1:ceil((i+3)/2)))+cnoutiso_log_C(sims,i+1,i+4);
end
nalk=nalk';
wetness=sum(nalk(2:5))/sum(nalk(1:5));
dryness=nalk(1)/(nalk(2)+nalk(3))


figure()   % the number of isotopic species generated
subplot(1,2,1)
hold on
bar(categorical({'C1','C2','C3','C4','C5','C6','C7','C8'}),[methane_13C; ethane_13C1; cnoutiso_log_C(sims,:,1)']);
title('13C-species')
subplot(1,2,2)
hold on
b=bar(categorical({'C1','C2','C3','C4','C5','C6','C7','C8'}),[methane_D; ethane_D1; cnoutiso_log_D(sims,:,1)']);
title('D-species')

figure()
plot([1:sims], 1000*ethane_13Cclumplog-1000);
yax=ylim;
hold on
title('\Delta 13C-13C ethane')
text(sims*0.75,yax(2)*0.75,strcat('n=',num2str(ethane_13C2_log(sims))))

figure()
subplot(1,2,1)
plot([1:sims],prop_sp_C_log)
title('propane \epsilon13Cc-t')
subplot(1,2,2)
plot([1:sims],prop_sp_D_log)
title('propane \epsilonDc-t')

%% for plotting the butane cracking test
figure()
subplot(1,2,1)
sc=scatter([1:sims]*rep,butane_sp_C_log, 'b','LineWidth',1);
hold on
p11=plot([-10, sims+100]*280,[14.4, 14.4],'r','LineWidth',3);
p11.Color(4)=0.75;
plot([-10, sims+100]*280,[14.4+Cerr_butane, 14.4+Cerr_butane],'k--','LineWidth',1.5);
plot([-10, sims+100]*280,[14.4-Cerr_butane, 14.4-Cerr_butane],'k--','LineWidth',1.5);
xlim([-3000 150000]);
title('butane \epsilon^{13}Cc-t')
ylim([10 25]);
xlabel('#simulations')
ylabel('\epsilon')
box on 
subplot(1,2,2)
title('butane \epsilonDc-t')
scatter([1:sims]*rep,butane_sp_D_log,'b', 'LineWidth',1)
hold on
title('butane \epsilonDc-t')
p=plot([-10,sims+100]*280,[93.0,93.0],'r','LineWidth',3);
p.Color(4)=0.75;
plot([-10,sims+100]*280,[93.0+Derr_butane,93.0+Derr_butane],'k--','LineWidth',3/2);
plot([-10,sims+100]*280,[93.0-Derr_butane,93.0-Derr_butane],'k--','LineWidth',3/2);
xlim([-3000 150000]);
xlabel('#simulations')
ylabel('\epsilon')
box on

%%
figure()
for i=0:5
    nC=i+3;
    nD=6+2*(i+1);
    R13C(i+3)=sum(cnoutiso_log_C(sims,i+1,1:ceil((i+3)/2)))/cnoutiso_log_C(sims,i+1,i+4)/nC;
    RD(i+3)=sum(cnoutiso_log_D(sims,i+1,1:ceil((i+3)/2)))/cnoutiso_log_D(sims,i+1,i+4)/nD;
    subplot(6,2,2*i+1)
    %s1=plot([1:i+3],squeeze(D_clean(sims,i+1,1:i+3)/2/cnoutiso_log_D(sims,i+1,i+4)/(R0/(1-100/1000))*1000-1000), 'o-','LineWidth', .75);
    s1=errorbar([1:i+3], squeeze(D_clean(sims,i+1,1:i+3)/2/cnoutiso_log_D(sims,i+1,i+4)/(R0/(1-100/1000))*1000-1000),1./sqrt(reshape(cnoutiso_log_D(sims,i+1,1:i+3),1,i+3))*1000, 'o-','LineWidth', .75);
    s1.MarkerFaceColor=[0.75 0.75 0.75 ];
    s1.CapSize=0;
    s1.Color=[.05 0 0];
    xlim([0 8])
    if i==5
        xlabel('Position #')
    end
    if i==0
        title('Intramolecular \deltaD')
    end
    subplot(6,2,2*i+2)
    %s2=plot([1:i+3],squeeze(C_clean(sims,i+1,1:i+3)/2/cnoutiso_log_C(sims,i+1,i+4)/(R0_13C/(1-25/1000))*1000-1000), 'o-','LineWidth', .75);
    s2=errorbar([1:i+3], squeeze(C_clean(sims,i+1,1:i+3)/2/cnoutiso_log_C(sims,i+1,i+4)/(R0_13C/(1-25/1000))*1000-1000),1./sqrt(reshape(cnoutiso_log_C(sims,i+1,1:i+3),1,i+3))*1000, 'o-','LineWidth', .75);
    s2.MarkerFaceColor=[.75 .75 .75];
    s2.CapSize=0;
    s2.Color=[.05 0 0];
    xlim([0 8])    
    if i==5
        xlabel('Position #')
    end
    if i==0
        title('Intramolecular \delta^{13}C')
    end
    hold on
end
delta13C=1000*R13C/R0_13C*(1-25/1000)-1000;
deltaD=1000*RD/R0*(1-100/1000)-1000;

figure()
subplot(1,2,1);
bar(delta13C);
title('\delta13C')
subplot(1,2,2);
bar(deltaD);
title('\deltaD')

C13_chung=[1,1/2,1/3,1/4,1/5,1/6,1/7,1/8];
D_chung=[3/4,1/2,3/8,3/10,3/12,3/14,3/16,3/18];
D_chung=[1/2,1/3,1/4,1/5,1/6,1/7,1/8,1/9];

figure()
subplot(1,2,1)
scatter(C13_chung,delta13C)
title('^{13}C Chung plot')
subplot(1,2,2)
scatter(D_chung,deltaD)
title('\deltaD Chung plot')

prop_sp_C_log(end)
prop_sp_D_log(end)
