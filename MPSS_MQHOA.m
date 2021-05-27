%------------------------------------------------------------------
% Multiple-population-based superposition state MQHOA
% MPSS-MQHOA algorithm for unimodal global optimization
% Matlab Code by Gang Xin, Peng Wang, Yuwei Jiao (10.05.2021).
% All rights reserved by Parallel Computing Lab.
%
%------------------------------------------------------------------
 function [gbest,gbestval,fitcount]= MPSS_MQHOA(fhd,Dimension,subpop_input,Max_Gen,VRmin,VRmax,varargin)
%clear all;
feature jit off;
format long;
%global DIM;
%--Parameters definition------------------------------%
%repeat          Repeat (default repeat=1)
%DIM             Dimension
%sigmaMin        Minimum search scale
%subpop          Number of subpopulations
%minDomain       Lower bound
%minDomain       Upper bound
%------------------------------------------------------%
repeat=1;
DIM=Dimension;
subpop=subpop_input;
minDomain=VRmin; maxDomain=VRmax;
maxFE=Max_Gen;
swarmNum=300;    % entire population size
%value_min=1/swarmNum;
ranking_coef=3;  % Adjust the size of the optimal dynamic population, the larger, the higher the peak of the population [1,2,3]
%-----------------------------------------------------%

%------Output-----------------------------------------%
%gbestV          Best fitness value
%gfe             Function evaluations            
%tot_time        Time cost
%-----------------------------------------------------%
gbestV=zeros(1,repeat);  
gfe=zeros(1,repeat); 
tot_time = zeros(1,repeat);
%-----------------------------------------------------%

for rep=1:repeat
    tic;
%--Some temporary variables---------------------------%
%funcV             Save the k sampled results (Vector£©
%temp_samplePos    Save the temperary sampled results for subpopulation£¨Matrix£©
%optimalSolution   Save the positions for all optimal solutions of the subpopulation(Matrix£©
%sigma             Search range
%stdPre            Standard deviation of sampling points in each dimension£¨Vector£©
%evolutiontime     Evaluation of the object functions
%-------------------------------------------------------%
    funcV=zeros(1,subpop);
    sigma=maxDomain-minDomain;
    stdPre=zeros(DIM,1);
    optimalSolution=unifrnd(minDomain,maxDomain,DIM,subpop);
    stdPre=std(optimalSolution,1,2) ;
    evolutiontime=0;
    for k=1:subpop 
        funcV(k)=feval(fhd,optimalSolution(:,k),varargin{:});
        evolutiontime=evolutiontime+1;  
    end
%-----------------------------------------------%
        while (evolutiontime<maxFE) % Stop criteria 
           while (evolutiontime<maxFE) 
                change_flag=0; % Stable flag          
%---------------Search scale matrix------------% 
            cov=diag(sigma.^2*ones(DIM,1)); 
%------------------Fitness Ranking-------------%
           [~,fitness_temp]=sort(funcV,'descend'); 
           [~,fitness_index]=sort(fitness_temp);
            sum_ranking=0; 
            for k=1:subpop
                sum_ranking=sum_ranking+(fitness_index(k))^(ranking_coef); 
            end           
                          
            for k=1:subpop 
%Adjust the size of the subpopulation dynamically-%               
               fitness=((fitness_index(k))^(ranking_coef))/(sum_ranking); 
               dy_swarmNum=(fitness)*swarmNum;
               dy_swarmNum=round(dy_swarmNum);
               if dy_swarmNum==0
                   dy_swarmNum=1;
               end
%-Sampling for the ith subpopulation--------------%
                    if isnan(dy_swarmNum)
                        dy_swarmNum=round(swarmNum/subpop); 
                    end
                    temp_samplePos=mvnrnd(optimalSolution(:,k),cov,dy_swarmNum); 
                    temp_samplePos=temp_samplePos';  
%-Infeasible solution processing-----------------%                    
                    for d=1:DIM
                        for L=1:dy_swarmNum
                            if temp_samplePos(d,L)>maxDomain
                                temp_samplePos(d,L)=minDomain+rand.*(maxDomain-minDomain);                 
                            end
                            if temp_samplePos(d,L)<minDomain
                                temp_samplePos(d,L)=minDomain+rand.*(maxDomain-minDomain);  
                            end
                        end
                    end
%--Comparing and replacing the worse soluton-----------%
                    for L=1:dy_swarmNum 
                        sampleValue=feval(fhd,temp_samplePos(:,L),varargin{:});
                        evolutiontime=evolutiontime+1; 
                        if sampleValue<funcV(k)
                            funcV(k)=sampleValue;
                            optimalSolution(:,k)=temp_samplePos(:,L); 
                            change_flag=1;
                        end
                        
                    end  
            end % End the evolution of all subpopulation
            if(change_flag==0)
                  break; 
            end                           
%--Replace the worst soution with mean position--------%
           meanPos=mean(optimalSolution,2);
           [~,index_max]=max(funcV);
           optimalSolution(:,index_max)=meanPos;        
           funcV(index_max)=feval(fhd,meanPos,varargin{:});
           evolutiontime=evolutiontime+1;           
%------------------------------------------------------%
           stdPre=std(optimalSolution,1,2); 
           if max(stdPre)<sigma % ground-state criteria
               break;
           end
           end % End of the evolution of current search scale
           sigma=sigma/1.2; % Annealing process
        end % End of the optimization

    
%--Save the results-----------------------------------%
    tot_time(rep) = toc; 
    [global_min,index]=min(funcV);
    gbestV(rep)=min(funcV);
    gfe(rep)=evolutiontime;
%----------------------------------------------------%
 end % Repeat only once in subprogram MPSS_MQHOA.m
    fitcount=evolutiontime;
    gbestval=global_min;
    gbest=(optimalSolution(:,index))'; 
