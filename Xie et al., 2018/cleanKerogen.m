%clean up kerogen molecules to keep only the large chunks
filename='EFK_50A_140'
load(filename)
BoG=graph(BoM, 'upper','omitselfloops'); 
splitg=conncomp(BoG,'OutputForm','cell');
freelist=[];  %small molecules that needs to be removed
for i=1:length(splitg)
    if length(splitg{i})<10
        freelist=[freelist splitg{i}];
    end
end
BoM(freelist,:)=[];
BoM(:,freelist)=[];
BoG=graph(BoM, 'upper','omitselfloops'); 
atoms(freelist)=[];
ringinfo(freelist)=[];
save(strcat(filename,'_clean.mat'),'BoM','ringinfo','atoms');
