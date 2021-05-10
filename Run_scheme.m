 clear all
 clc
 
%  mex cec14_func.cpp -DWINDOWS

D=10;                                                                      % Dimension of the problem
Xmin=-100; Xmax=100;                                                       % Range of the search space
G_min=[-1400,-1300,-1200,-1100,-1000,-900,-800,-700,-600,-500,-400,-300,-200,-100,100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400];

Subpopulations=30;                                                         % Number of the subpopulations
iter_max=1e4*D;                                                            % MaxFES
runs=25;                                                                   % Repeats
funcnum=28;                                                                % Selection of the test function
fhd=str2func('cec13_func');                                                % CEC 2013 benchmark function suit
Algstr={'MPSS_MQHOA','cm_fr'};
Algnum=1;                                                                  % Selected algorithm
namestr=['.\data\',Algstr{Algnum},['_cec13_result_D',num2str(D),'_P',num2str(Subpopulations),'_',datestr(now, 'yyyymmddHHMMSS')]];

txtstr=[namestr,'.txt'];
fp=fopen(txtstr,'w+');
for i=1:funcnum    
    func_num=i;  
    fprintf('%s \n',datestr(now, 'yyyy/mm/dd|HH:MM:SS'));
    fprintf(fp,'%s \r\n',datestr(now, 'yyyy/mm/dd|HH:MM:SS'));
    for j=1:runs 
        tic;
            switch(Algnum) 
                case 1
                    [gbest,gbestval,FES]=MPSS_MQHOA(fhd,D,Subpopulations,iter_max,Xmin,Xmax,func_num);
                case 2
                    [gbest,gbestval,FES]=mqa_cm_fr(fhd,D,Subpopulations,iter_max,Xmin,Xmax,func_num);
            end
            
        tot_time(j) = toc;                                                 % Time cost
        xbest(i,:)=gbest;                                                  % Position of the global optima
        fbest(i,j)=gbestval;                                               % Fitness value of the glocal optima
        Err(i,j)=fbest(i,j)-G_min(i);                                      % Errors
        FESM(i,j)=FES;                                                     % Function evaluations
        fprintf('Fuc= f%d,  No. %d run, DIM=%d, FE=%d, time=%d, error=%e\n',i,j,D,FES,toc,fbest(i,j)-G_min(i));
        fprintf(fp,'Fuc= f%d,  No. %d run, DIM=%d, Global minimum=%e. FE=%d, time=%d, error=%e\r\n',i,j,D,fbest(i,j),FES,toc,fbest(i,j)-G_min(i));
    end 
   % Print to screen
   fprintf('---------------------------------------------------------------\n');        
   fprintf('Repeat=%d, Mean FE=%1.2e,Meantime=%1.2e\n',runs,mean(FES),mean(tot_time));        
%    fprintf('MeanValue=%1.2e, BestValue=%1.2e, Std=%1.2e, \n',mean(fbest(i,:)),min(fbest(i,:)),std(fbest(i,:)));        
%    fprintf('MeanErr=%1.2e, BestErr=%1.2e, StdErr=%1.2e, \n',mean(Err(i,:)),min(Err(i,:)),std(Err(i,:)));
   fprintf('%s \n',datestr(now, 'yyyy/mm/dd|HH:MM:SS'));
   % Print to 'txt' file 
   fprintf(fp,'---------------------------------------------------------------\r\n');
   fprintf(fp,'Repeat=%d, Mean FE=%1.2e,Meantime=%1.2e\r\n',runs,mean(FES),mean(tot_time));        
   fprintf(fp,'MeanValue=%1.2e, BestValue=%1.2e, Std=%1.2e, \r\n',mean(fbest(i,:)),min(fbest(i,:)),std(fbest(i,:)));        
   fprintf(fp,'MeanErr=%1.2e, BestErr=%1.2e, StdErr=%1.2e, \r\n',mean(Err(i,:)),min(Err(i,:)),std(Err(i,:)));
   fprintf(fp,'%s \r\n',datestr(now, 'yyyy/mm/dd|HH:MM:SS'));  
   
end
% Data Processing


% Mean error
for i=1:funcnum
meanErr(i)=mean(Err(i,:));
end 
% Std error
for i=1:funcnum
stdErr(i)=std(Err(i,:));
end 
% Save to 'xlsx' file
xlsxstr=[namestr,'.xlsx'];
xlswrite(xlsxstr,meanErr,'meanErr');
xlswrite(xlsxstr,stdErr,'stdErr');
staterr=[meanErr',stdErr'];
xlswrite(xlsxstr,staterr,'staterr');
fclose(fp);
