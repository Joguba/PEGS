% Jonathan's Photoelastic Disk Solver
%
% Takes input from preprocessing script joDiskPrep.m as of 2016/09/27
% inspired from peDiskSolve by James Puckett (Phd-Thesis 2012)
% http://nile.physics.ncsu.edu
%
% If you use this sovler please cite the follwoing paper
% K.E. Daniels, J. E. Kollmer & J. G. Puckett, "Photoelastic force measurements in granular materials", Rev. Sci. Inst. (201X)
% DOI: XXXXXX
%
% The generation of the synthetic force images is parallelized 
% in a way that each row of the ouput image can be its own worker/thread
% it is usually limited by your number of processing cores though. 
%
% Running this script may take a long time (up to one minute per particle
% on a 2016 desktop computer) so it makes sense to send this off to a high
% performance computer. 
%
% Depending on your experimental data you can play around with the fit
% options, which might speed up processing considerabely.
%
% Last edit on 2016/09/28 by Jonathan Kollmer (jekollme@ncsu.edu)

function PeGSDiskNewton3HPC(fileName)


status = 0;
while status == 0  
    [status,~] = license('checkout','Optimization_Toolbox');
    % early break (saves 60 sec. if checkout works on 1st try)
    if status == 1
        continue
    end
    % otherwise wait
    pause(60)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%User defined values:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%which files are we processing ?
%directory = '/Users/jekollmer/Desktop/TestEphraim/';
%directory = '/eno/jekollme/DATA/arneReanalyzed/data01/';
%directory = 'DATA/test/';
%directory = '/eno/jekollme/DATA/20160711/pegs/step09/';
%directory = '/eno/jekollme/DATA/20170616/03-10x12/solved/';
%directory = '/eno/jekollme/DATA/Gardner/20170616/hpc3/03/';
%directory = '/eno/jekollme/DATA/20170616/04-15x12/solved/';

directory = './';
%files = dir([directory, '*Tracked.mat']); 
files = dir('*solved.mat'); 
%files = dir([directory, '*newtoCleaned.mat']); 

%directory = '/eno/jekollme/DATA/20170206/2/Cycle07/';
%files = dir([directory, 'DSC*_adjusted.mat']); 

%how much of the particle diameter is used to fit the synthetic image 
%(1 = use everything). Change this parameter only if the fit doesn't work 
%correctly because of imaging artifacts at the particle boundary.
maskradius = 0.96;% 
scaling = 1; %scale the image by this factor before doing the fit

%do we want to see each particle on screen while it is fitted ?
verbose = false; 
errorThreshold = 1000; %if the fit error is lager than this, completely abandon the forces and replace by opposing ones.

swap1=0;
swap2=0;
sub1=0;
sub2=0;
swap1vec=[];
swap2vec=[];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%There should be no need for user input below this line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
nFrames = length(files); %how many files are we processing ?
for frame = 1:nFrames %loop over these frames 

    fileName = [directory,files(frame).name]; %which file/frame are we processing now ?
    maskradius = maskradius / 2; %I did an unwise choice in naming this
    load(fileName); %load the particle data file to process
    pres = PeGS_SortByID(pres) ;
    particle=pres; 
    N = length(particle); %number of particles in this frame
    
    display(['processing file ',fileName, 'containing ' ,num2str(N), 'particles']); %status indicator
    for annealingCycle=1:2
    for n=1:N
        %keyboard
        display(['frame: ',num2str(frame),', checking force(s) on particle ',num2str(n)]); %status indicator
        if (particle(n).z > 1 )
                    %bookkeeping
                    fsigma = particle(n).fsigma;
                    %fsigma = 100;
                    rm = particle(n).rm;

%                     %This is the Camera Image
%   %template = im2double(particle(n).forceImage);
                     template = particle(n).forceImage;
                     template = imresize(template,scaling);
%                     template = template-0.1;
%                     %template = (template -0.2); %fine tunes the image, this should happen in preprocessing!
%                     template = template.*(template > 0); %fine tunes the image, this should happen in preprocessing!
%                     template = template*3; %fine tunes the image, this should happen in preprocessing!

                    %size of the force image
                    px = size(template,1); 
                    
                    %plot the experimental image that is going to be fitted
                    %onto
                    if verbose
                        figure(1)
                        subplot(1,2,1)
                        title('camera image');
                        imshow(template);
                    end
                    
                    %Create initial guesses for each contact force, based
                    %on the gradient squared in the contact region (also
                    %better guess lower than too high)
                    z = particle(n).z;
                    forces = zeros(z,1);
                    
                    %Fitting results
                    forces = particle(n).forces;
                    alphas = particle(n).alphas;
                    beta = particle(n).betas;

                    %apply force balance to the initial guesses
                    %[alphas,forces] = forceBalance(forces,alphas,beta);
                    

                    %create a circular mask
                    cx=px/2;cy=px/2;ix=px;iy=px;r=maskradius*px;
                    [x,y]=meshgrid(-(cx-1):(ix-cx),-(cy-1):(iy-cy));
                    c_mask=((x.^2+y.^2)<=r^2);   

                    %This is the function we want to fit, i.e. a synthetic
                    %version of the force image with free parameters f (here called par(1:z)) and
                    %alpha (here called par(z+1:z*z) since we have to stuff
                    %them into one vector.
                    func = @(par) joForceImgNoFB(z, par(1:z),par(z+1:z+z), beta(1:z), fsigma, rm, px, verbose); %+par(2*z+1); %this is the function I want to fit (i.e. synthetic stres image), the fitting paramters are in vector par
                    %This is the error function we are actually fitting,
                    %that is, the distance between our fit function and the
                    %real particle image in terms of the sum of squares of the pixelwise differnce.
                    %Also a mask is applied to crop
                    %out only the circular part of the particle. 
                    err = @(par) sum(sum( ( (template-func(par)).^2) )); %BUG: for some reason I sometimes get imaginary results, this should not happen
                    
                    %Set up initial guesses
                    p0(1:z) = forces;
                    p0(z+1:2*z) = alphas;

                    %resudual
                    fitError = err(p0);
                    score(1) = fitError;
                    bestFit = fitError;
                    
                    %generate an image with the original parameters
                    imgFit = joForceImgNoFB(z, forces, alphas, beta, fsigma, rm, px*(1/scaling), verbose);
                    

                    z = length(forces);
                    for m = 1:z
                        
                        forces = particle(n).forces; %take the original forces again
                        alphas = particle(n).alphas;

                        if(verbose)
                            figure(2)
                            subplot(z,8,(m-1)*8+1)
                            imshow(template)
                            title('Camera')
                            xlabel('0')
                            subplot(z,8,(m-1)*8+2)
                            imshow(imgFit)
                            title('PEGS fit')
                            xlabel(['residual: ',num2str(score(1))],'FontSize',20);
                        end
                        otherParticle = particle([particle.id]==particle(n).neighbours(m));
                        reciprocalContact = ([otherParticle.neighbours] == particle(n).id);
                        if ~isempty(reciprocalContact)
                        if ~isempty(otherParticle.forces)
                        
                            %%% SWAP CONTACT FORCE
                            forces(m) = otherParticle.forces(reciprocalContact);
                            alphas(m) = otherParticle.alphas(reciprocalContact);

                            p0(1:z) = forces;
                            p0(z+1:2*z) = alphas;
                            score(2) = err(p0);                  
                            imgFit2 = joForceImgNoFB(z, forces, alphas, beta, fsigma, rm, px*(1/scaling), verbose);
                            if(verbose)
                                figure(2)
                                subplot(z,8,(m-1)*8+3)
                                imshow(imgFit2)
                                title('F_{AB}->F_{BA}')
                                xlabel(['residual: ',num2str(score(2))],'FontSize',20);
                                drawnow
                            end
                            if score(2)<bestFit %Did this give a better fit than the original solution ?
                                pres(n).forces(m) = forces(m);
                                pres(n).alphas(m) = alphas(m);
                                bestFit = score(2);
                                display(['swapping contact force',num2str(m), ' in annealing cycle ', num2str(annealingCycle)]);
                                swap1=swap1+1;
                            end
                            
%                             %%% USE MEAN CONTACT FORCE
%                             forces(m) = (forces(m) + otherParticle.forces(reciprocalContact) ) /2;
%                             alphas(m) = (alphas(m) + otherParticle.alphas(reciprocalContact) ) /2;
% 
%                             p0(1:z) = forces;
%                             p0(z+1:2*z) = alphas;
%                             score(3) = err(p0);                  
%                             imgFit3 = joForceImgNoFB(z, forces, alphas, beta, fsigma, rm, px*(1/scaling), verbose);
%                             if(verbose)
%                                 figure(2)
%                                 subplot(z,8,(m-1)*8+4)
%                                 imshow(imgFit3)
%                                 title('(F_{AB} + F_{BA})/2')
%                                 xlabel(['residual: ',num2str(score(3))],'FontSize',20);
%                                 drawnow
%                             end
%                             
%                             if score(3)<bestFit %Did this give a better fit than the original solution or the one before ?
%                                 pres(n).forces(m) = forces(m); %TODO: Think about weather to carry this over to the next particle or not (this right now does not, would need to sawp pres to particle here
%                                 pres(n).alphas(m) = alphas(m);
%                                 bestFit = score(3);
%                                 display(['substituting mean contact force ',num2str(m), ' in annealing cycle ', num2str(annealingCycle)]);
%                                 sub1=sub1+1;
%                             end
                            
                            %%% NOW EVERYTHING AGAIN BUT USE ALREADY UPDATED CONTACT FORCE
                            
                            otherParticle = pres([particle.id]==particle(n).neighbours(m));
                            forces = pres(n).forces;
                            alphas = pres(n).alphas;

                            %%% SWAP CONTACT FORCE
                            forces(m) = otherParticle.forces(reciprocalContact);
                            alphas(m) = otherParticle.alphas(reciprocalContact);

                            p0(1:z) = forces;
                            p0(z+1:2*z) = alphas;
                            score(4) = err(p0);                  
                            imgFit2 = joForceImgNoFB(z, forces, alphas, beta, fsigma, rm, px*(1/scaling), verbose);
                            if(verbose)
                                figure(2)
                                subplot(z,8,(m-1)*8+5)
                                imshow(imgFit2)
                                title('Updated F_{AB}->F_{BA}')
                                xlabel(['residual: ',num2str(score(4))],'FontSize',20);
                                drawnow
                            end
                            if score(4)<bestFit %Did this give a better fit than the original solution ?
                                pres(n).forces(m) = forces(m);
                                pres(n).alphas(m) = alphas(m);
                                bestFit = score(4);
                                display(['swapping contact force from updated contact',num2str(m), ' in annealing cycle ', num2str(annealingCycle)]);
                                swap2=swap2+1;
                            end
                            
%                              %%% USE MEAN CONTACT FORCE
%                             forces(m) = (forces(m) + otherParticle.forces(reciprocalContact) ) /2;
%                             alphas(m) = (alphas(m) + otherParticle.alphas(reciprocalContact) ) /2;
% 
%                             p0(1:z) = forces;
%                             p0(z+1:2*z) = alphas;
%                             score(5) = err(p0);                  
%                             imgFit3 = joForceImgNoFB(z, forces, alphas, beta, fsigma, rm, px*(1/scaling), verbose);
%                             if(verbose)
%                                 figure(2)
%                                 subplot(z,8,(m-1)*8+6)
%                                 imshow(imgFit3)
%                                 title('Updated (F_{AB} + F_{BA})/2')
%                                 xlabel(['residual: ',num2str(score(5))],'FontSize',20);
%                                 drawnow
%                             end
%                             
%                             if score(5)<bestFit %Did this give a better fit than the original solution ?
%                                 pres(n).forces(m) = forces(m);
%                                 pres(n).alphas(m) = alphas(m);
%                                 bestFit = score(5);
%                                 display(['substituting mean contact force from upadted contact ',num2str(m), ' in annealing cycle ', num2str(annealingCycle)]);
%                                 sub2=sub2+1;
%                             end
                            
                        end
                            
                        else
                           display('Warning: opposing contact not found');
                        end
                    
                            p0(1:z) = pres(n).forces;
                            p0(z+1:2*z) = pres(n).alphas;
                            finalScore = err(p0); 

                            imgFit4 = joForceImgNoFB(z, pres(n).forces, pres(n).alphas, beta, fsigma, rm, px*(1/scaling), verbose);
                            pres(n).synthImg = imgFit4;
                            pres(n).fitError = finalScore;
                            if (finalScore<1)
                                %keyboard
                                display('this fit seems too good to be true');
                            end
                            if(verbose)
                                figure(2)
                                subplot(z,8,(m-1)*8+8)
                                imshow(imgFit4)
                                title('Best')
                                xlabel(['residual: ',num2str(pres(n).fitError)],'FontSize',20);
                                drawnow
                                pause(0.1)
                            end
                            
                    
                    end
                    
                    
                    
                    %generate an image with the new parameters
                    %figure(3)
                    %imgFit = joForceImgNoFB(z, forces, alphas, beta, fsigma, rm, px*(1/scaling), verbose);
                    
                    %since the image gnerator also forces force balance we
                    %have to explicitly do it once more to the values we
                    %are gonna save 
                    %[alphas,forces] = forceBalance(forces,alphas,beta);
                    
%                     %store the new information into a new particle vector 
%                     pres(n).fitError = fitError;
%                     pres(n).forces = forces;
%                     pres(n).alphas = alphas;
%                     pres(n).synthImg = imgFit;
        end
    end
    %save the result
    particle=pres;
        swap1vec(annealingCycle)=swap1;
        swap1=0;
        swap2vec(annealingCycle)=swap2;
        swap2=0;
    end
    save([fileName(1:end-10),'newtonized2.mat'],'pres');
    
    swap1vec
    swap2vec
    
    
    
%end

end

%%% 
                           %%%%% DO WE TRUST THE CURRENT FIT OR RATHER
%                             %%%%% SWAP COMPLETELY ?
%                              if(bestFit > errorThreshold && ~swapped)
%                                    
%                                     display('fit error too high, completely swapping forces');
%                                     swapped = true;
%                                     
%                                     forces(m) = otherParticle.forces(reciprocalContact);
%                                     alphas(m) = otherParticle.alphas(reciprocalContact);
% 
%                                     p0(1:z) = forces;
%                                     p0(z+1:2*z) = alphas;
%                                     score(6) = err(p0);                  
%                                     imgFit6 = joForceImgNoFB(z, forces, alphas, beta, fsigma, rm, px*(1/scaling), verbose);
%                                     if(verbose)
%                                         figure(2)
%                                         subplot(z,8,(m-1)*8+7)
%                                         imshow(imgFit6)
%                                         title('FullContactSwap')
%                                         xlabel(['residual: ',num2str(score(5))],'FontSize',20);
%                                         drawnow
%                                     end
% 
%                                     pres(n).forces(m) = forces(m);
%                                     pres(n).alphas(m) = alphas(m);
%                                     bestFit = score(6);
% 
%                              end