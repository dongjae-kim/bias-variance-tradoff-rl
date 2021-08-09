function [out]=SIMUL_Arbitration_v3_GA3(job_opt)

% SpeedyGA is a vectorized implementation of a Simple Genetic Algorithm in Matlab
% Version 1.3
% Copyright (C) 2007, 2008, 2009  Keki Burjorjee
% Created and tested under Matlab 7 (R14).

%  Licensed under the Apache License, Version 2.0 (the "License"); you may
%  not use this file except in compliance with the License. You may obtain
%  a copy of the License at
%
%  http://www.apache.org/licenses/LICENSE-2.0
%
%  Unless required by applicable law or agreed to in writing, software
%  distributed under the License is distributed on an "AS IS" BASIS,
%  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%  See the License for the specific language governing permissions and
%  limitations under the License.

%  Acknowledgement of the author (Keki Burjorjee) is requested, but not required,
%  in any publication that presents results obtained by using this script

%  Without Sigma Scaling, Stochastic Universal Sampling, and the generation of mask
%  repositories, SpeedyGA faithfully implements the specification of a simple genetic
%  algorithm given on pages 10,11 of M. Mitchell's book An Introduction to
%  Genetic Algorithms, MIT Press, 1996). Selection is fitness
%  proportionate.

%clear all
close all
warning('off')




%% Parameters
Opt.SIMUL_PC=2; % 1:local   2:gpu
if(Opt.SIMUL_PC==1)    
Opt.str_path_div='\';   
addpath('/home/swlee/simul/Model_RL');
else     
Opt.str_path_div='/'; 
end


% 0.experiment
Opt.exp_name0=job_opt.name;
Opt.TARGET_FOLDER_NAME='result_simul';
Opt.ind_included=job_opt.list; %job_opt=[1:3,5:12,14:20];
Opt.EXP_NAME=['ANY_start_' sprintf('%02d_Sbj%02d',Opt.exp_name0,Opt.ind_included)];
% 1.parameters 
% [PE_tolerance_m1, Time_Step,...
% A_12, B_12, A_21, B_21,a...
% param_sarsa/fwd.alpha, myArbitrator.p, .tau(softmax),...
% myArbitrator.ind_active_model, PE_tolerance_m2]

% [NOTE] the codingbition for param_fwd.alpha is disabled. Instead,
% 'param_fwd.alpha'='param_sarsa.alpha'
Opt.IsCodingLog=[0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]; %Opt.IsCodingLog=[0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0];
Opt.bits=[4, 12, 6, 6, 6, 6, 6, 4, 6, 1, 4];
minval=[0.1, 1e-2, 0.2, .1e1, 1.5, 1.3e1, 0.1, 1, 0.1, 1, 3];
maxval=[0.5, 1e2, 2.5, 1.5e1, 5.0, 3.0e1, 1.0, 1e1, 6.0, 2, 10];

%% automatically determined
if(Opt.SIMUL_PC==1)
    Opt.path_ext=['F:' Opt.str_path_div '0-Program' Opt.str_path_div 'MATLABwork' Opt.str_path_div 'work' Opt.str_path_div 'ModelRL'];
else
    Opt.path_ext=['/home/swlee/simul/ModelRL'];
end
Opt.minmax=[minval; maxval];
%if(exist([Opt.path_ext Opt.str_path_div Opt.TARGET_FOLDER_NAME Opt.str_path_div Opt.EXP_NAME],'dir')==7)
%    rmdir([Opt.path_ext Opt.str_path_div Opt.TARGET_FOLDER_NAME Opt.str_path_div Opt.EXP_NAME],'s');
%end
mkdir([Opt.path_ext Opt.str_path_div Opt.TARGET_FOLDER_NAME],Opt.EXP_NAME);


%% Parameters for GA
len1=sum(Opt.bits); % The length of the genomes
popSize1=100;   % The size of the population (must be an even number)
maxGens=500;              % The maximum number of generations allowed in a run
probCrossover1=0.1;           % The probability of crossing over.
probMutation1=0.05;        % The mutation probability (per bit)
sigmaScalingFlag=1;        % Sigma Scaling is described on pg 168 of M. Mitchell's
% GA book. It often improves GA performance.
sigmaScalingCoeff=0.5;       % Higher values => less fitness pressure

SUSFlag=1;                 % 1 => Use Stochastic Universal Sampling (pg 168 of
%      M. Mitchell's GA book)
% 0 => Do not use Stochastic Universal Sampling
%      Stochastic Universal Sampling almost always
%      improves performance

crossoverType=2;           % 0 => no crossover
% 1 => 1pt crossover
% 2 => uniform crossover

visualizationFlag=0;       % 0 => don't visualize bit frequencies
% 1 => visualize bit frequencies

verboseFlag=0;             % 1 => display details of each generation
% 0 => run quietly

useMaskRepositoriesFlag=1; % 1 => draw uniform crossover and mutation masks from
%      a pregenerated repository of randomly generated bits.
%      Significantly improves the speed of the code with
%      no apparent changes in the behavior of
%      the SGA
% 0 => generate uniform crossover and mutation
%      masks on the fly. Slower.





% crossover masks to use if crossoverType==0.
mutationOnlycrossmasks1=false(popSize1,len1);

% pre-generate two “repositories?of random binary digits from which the
% the masks used in mutation and uniform crossover will be picked.
% maskReposFactor determines the size of these repositories.

maskReposFactor1=5;
uniformCrossmaskRepos1=rand(popSize1/2,(len1+1)*maskReposFactor1)<0.5;
mutmaskRepos1=rand(popSize1,(len1+1)*maskReposFactor1)<probMutation1;

% preallocate vectors for recording the average and maximum fitness in each
% generation
avgFitnessHist=zeros(1,maxGens+1);
maxFitnessHist=zeros(1,maxGens+1);


eliteIndiv1=[];
eliteFitness=-realmax;


% the population is a popSize by len matrix of randomly generated boolean
% values
pop1=rand(popSize1,len1)<.5;

Hist_fitnessVals=[];

for gen=0:maxGens
    
    Opt.generation_id=gen+1;
    disp(sprintf('[%s] %d-st generation:: fitness value computation... #####################################',Opt.EXP_NAME,gen+1));
    
    % evaluate the fitness of the population. The vector of fitness values
    % returned must be of dimensions "1 x popSize".
    fitnessVals=test_fit_v3(pop1,Opt);
    
    Hist_fitnessVals=[Hist_fitnessVals; fitnessVals];
    
    [maxFitnessHist(1,gen+1),maxIndex]=max(fitnessVals);
    avgFitnessHist(1,gen+1)=mean(fitnessVals);
    if eliteFitness<maxFitnessHist(gen+1)
        eliteFitness=maxFitnessHist(gen+1);
        eliteIndiv1=pop1(maxIndex,:);
    end
    
    % display the generation number, the average Fitness of the population,
    % and the maximum fitness of any individual in the population
    %     if verboseFlag
    %         display(['gen=' num2str(gen,'%.3d') '   avgFitness=' ...
    %             num2str(avgFitnessHist(1,gen+1),'%3.3f') '   maxFitness=' ...
    %             num2str(maxFitnessHist(1,gen+1),'%3.3f')]);
    %     end
    % Conditionally perform bit-frequency visualization
    %     if visualizationFlag
    %         figure(1)
    %         set (gcf, 'color', 'w');
    %         hold off
    %         bitFreqs1=sum(pop1)/popSize1;     bitFreqs2=sum(pop2)/popSize2;
    %         plot(1:len,bitFreqs1, '.'); plot(1:len,bitFreqs2, '.');
    %         axis([0 len 0 1]);
    %         title(['Generation = ' num2str(gen) ', Average Fitness = ' sprintf('%0.3f', avgFitnessHist(1,gen+1))]);
    %         ylabel('Frequency of the Bit 1');
    %         xlabel('Locus');
    %         drawnow;
    %     end
    
    disp(sprintf('### %d-st generation:: selection/mutation/crossover/etc... #####################################',gen+1));
    
    % Conditionally perform sigma scaling
    if sigmaScalingFlag
        sigma=std(fitnessVals);
        if sigma~=0;
            fitnessVals=1+(fitnessVals-mean(fitnessVals))/...
                (sigmaScalingCoeff*sigma);
            fitnessVals(fitnessVals<=0)=0;
        else
            fitnessVals=ones(1,popSize1);
        end
    end
    
    
    % Normalize the fitness values and then create an array with the
    % cumulative normalized fitness values (the last value in this array
    % will be 1)
    cumNormFitnessVals=cumsum(fitnessVals/sum(fitnessVals));
    
    % Use fitness proportional selection with Stochastic Universal or Roulette
    % Wheel Sampling to determine the indices of the parents
    % of all crossover operations
    if SUSFlag
        markers1=rand(1,1)+[1:popSize1]/popSize1; 
        markers1(markers1>1)=markers1(markers1>1)-1;
    else
        markers1=rand(1,popSize1);
    end
    [temp1 parentIndices1]=histc(markers1,[0 cumNormFitnessVals]);    
    parentIndices1=parentIndices1(randperm(popSize1));    
    
    % deterimine the first parents of each mating pair
    firstParents1=pop1(parentIndices1(1:popSize1/2),:);    
    
    % determine the second parents of each mating pair
    secondParents1=pop1(parentIndices1(popSize1/2+1:end),:);    
    
    % create crossover masks
    if crossoverType==0
        masks1=mutationOnlycrossmasks1;
        masks2=mutationOnlycrossmasks2;
    elseif crossoverType==1
        masks1=false(popSize1/2, len1);
        masks2=false(popSize2/2, len2);
        temp1=ceil(rand(popSize1/2,1)*(len1-1));
        temp2=ceil(rand(popSize2/2,1)*(len2-1));
        for i=1:popSize1/2
            masks1(i,1:temp1(i))=true;
        end
        for i=1:popSize2/2
            masks2(i,1:temp2(i))=true;
        end
    else
        if useMaskRepositoriesFlag
            temp1=floor(rand*len1*(maskReposFactor1-1));
            masks1=uniformCrossmaskRepos1(:,temp1+1:temp1+len1);            
        else
            masks1=rand(popSize1/2, len1)<.5;            
        end
    end
    
    % determine which parent pairs to leave uncrossed
    reprodIndices1=rand(popSize1/2,1)<1-probCrossover1;    
    masks1(reprodIndices1,:)=false;
        
    % implement crossover
    firstKids1=firstParents1;    
    firstKids1(masks1)=secondParents1(masks1);
    secondKids1=secondParents1;    
    secondKids1(masks1)=firstParents1(masks1);    
    pop1=[firstKids1; secondKids1];
        
    % implement mutation
    if useMaskRepositoriesFlag
        temp1=floor(rand*len1*(maskReposFactor1-1));    
        masks1=mutmaskRepos1(:,temp1+1:temp1+len1);        
    else
        masks1=rand(popSize1, len1)<probMutation1;
    end
    pop1=xor(pop1,masks1);    
    
    % keep the best
    pop1(1,:)=eliteIndiv1;
end
% if verboseFlag
%     figure(2)
%     %set(gcf,'Color','w');
%     hold off
%     plot([0:maxGens],avgFitnessHist,'k-');
%     hold on
%     plot([0:maxGens],maxFitnessHist,'c-');
%     title('Maximum and Average Fitness')
%     xlabel('Generation')
%     ylabel('Fitness')
% end


out=1;
end
