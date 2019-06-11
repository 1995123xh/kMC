%make a parsing script to deal with data output of the newcollector
function result=parseoutput(file)
load(file)
prop_sp_C_log=C_clean(:,1,2)./C_clean(:,1,1)*1000-1000;
prop_sp_D_log=D_clean(:,1,2)./D_clean(:,1,1)*1000-1000;
R13C(1)=methane_13C/(methane_D+methane_n-methane_13CD);
RD(1)=methane_D/(methane_13C+methane_n-methane_13CD)/4;
R13C(2)=ethane_13C1/ethane_13Cn/2;
RD(2)=ethane_D1/ethane_Dn/6;

figure()   % the number of isotopic species generated
subplot(1,2,1)
hold on
bar(categorical({'C3','C4','C5','C6','C7','C8'}),cnoutiso_log_C(sims,:,1)')
title('13C-species')
subplot(1,2,2)
hold on
b=bar(categorical({'C3','C4','C5','C6','C7','C8'}),cnoutiso_log_D(sims,:,1)')
title('D-species')

figure()
plot([1:sims], 1000*ethane_13Cclumplog-1000)
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

figure()
for i=0:5
    nC=i+3;
    nD=6+2*(i+1);
    R13C(i+3)=sum(cnoutiso_log_C(sims,i+1,1:ceil((i+3)/2)))/cnoutiso_log_C(sims,i+1,i+4)/nC;
    RD(i+3)=sum(cnoutiso_log_D(sims,i+1,1:ceil((i+3)/2)))/cnoutiso_log_D(sims,i+1,i+4)/nD;
    subplot(6,2,2*i+1)
    scatter([1:i+3],squeeze(D_clean(sims,i+1,1:i+3)/mean(D_clean(sims,i+1,1:i+3))*1000-1000))
    subplot(6,2,2*i+2)
    scatter([1:i+3],squeeze(C_clean(sims,i+1,1:i+3)/mean(C_clean(sims,i+1,1:i+3))*1000-1000))
    hold on
end

delta13C=1000*R13C/R0_13C-1000;
deltaD=1000*RD/R0-1000;

figure()
bar(delta13C)
title('\delta13C')
figure()
bar(deltaD)
title('\deltaD')


