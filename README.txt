  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   
 / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \
      `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'
| |											      											  				  | | 
				+----------------------------------------------+                      		      
| |				|  +---------------------------------------+   |		      		  		  	  		  	  | |	
				|  |  MATLAB R2025a                        |   |		      
| |				|  | ------------------------------------- |   |		      		  		   	  		  	  | |
				|  | +---+-------------------------------+ |   |		       		      
| |				|  | |~~~|    Computation of periodic    | |   |		      		  		   	  		 	  | | 
				|  | |~~~|    solutions for the 3 x 3    | |   |		       
| |				|  | |~~~|  compressible Euler equations | |   |	    ________________________________	  | |  
				|  | |~~~|                               | |   |	   /                               /\
| |				|  | |~~~| [part of the dissertation of  | |   |      / part of the dissertation of  _/ /\    | |  
				|  | |~~~|   Andry Brinsko, University   | |   |	 /   Andry Brinsko, University    \/      
| |				|  | |~~~|   of Massachusetts Amherst]   | |   |    /    of Massachusetts Amherst    / \	  | |  
				|  | +---+-------------------------------+ |   |   /________________________________/ /       
| |				|  | |   | © Andry Brinsko, 2023-2025    | |   |   \________________________________\/	      | |  
				|  | +---+-------------------------------+ |   |	\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \	       
| |				|  +---------------------------------------+   |		      		  		   	  		  	  | |
				+----------------------------------------------+		    _________________________  
| |       	    |\  +------------------------------------------+ \		   /\     © Andry Brinsko    \  	  | |
				 \\  \_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\__\ \		  /_/\       2023-2025        \
| |	 		 	  \\  \___\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\____\ \	    \/\________________________\      | |
	  		  	   \\  \___\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\ ___  \ \	     \/________________________\      
| |	   		   		\\  \_____\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\ \_\ _\ \		 / / / / / / / / / / / / / /      | |
	    		     \\  \_\_\_\_\_\_\_________________\___\_\_\_\_ \ \		      
| |	     		      \\  +-----------------------------------------+  \		      	  		   	      	  | |
  	                   \\                                               \	      
| |               	    \\                                               \	      	  		  	  		  	  | |
					     \+----------------------------------------------+|	      
| |		 		 		  ------------------------------------------------	      	  		   	  		  	  | |

| |											      											  		  		  | |
      .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.  .-.-.   
 \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / /  
  `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'

For any questions or communications otherwise, send an emai



 _.--.__.-'""`-.__.--.__.-'""`-.__.--.__.-'""`-.__.--.__.-'""`-._.--.__.-'""`-.__.--.__.-'""`-.__.--.__.-'""`-._
 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
_/""7__/""7__/""7__/""7__/""7__/""7__/""7__/""7__/""7__/""7__/""7__/""7__/""7__/""7__/""7__/""7_/""7__/""7__/""7_
 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 "`--'""`-.__.-'""`--'""`-.__.-'""`--'""`-.__.-'""`--'""`-.__.-'"`--'""`-.__.-'""`--'""`-.__.-'""`--'""`-.__.-'"

          /\
         /  \
       /_ %%=O=%     _____________        
        %  - -%     |        '\\\\\\       
   _____c%   > __   |        ' ____|_      
  (_|. .  % ` % .'  |   +    '||::::::   To run a computation of a periodic solution of the Euler equations
   ||. ___)%%%%_.'  |        '||_____|    yourself, there is one main top-level script that allows for the user
   ||.(  \ ~ / ,)'  \'_______|_____|      to specify some of the main parameters of the problem at runtime:
   || /|  \'/  |\   ___/____|___\___     
  _,,,;!___*_____\_|    _    '  <<<:|     		%% \Euler time evolution\TimeEvoMain.m %% 
 /     /|          |_________'___o_o| 
/_____/ /   __________________________
|:____|/   | Running a computation of |
           |  a solution on your own  |
            ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
	⚠ NOTE ⚠ Currently, it is up to you to manually set α, the scaling coefficient that
	determines how big of a perturbation from background pressure we are starting with in

		p^{(0)}(x,0) = \bar{p} + (α \bar{p}) φ_k
	
	When setting α, keep in mind that in the author's expierence, the code begins have a 
	higher risk of breaking once α passes 0.4 to 0.5. This may seem rather slow for the code to
	start falling apart, but such a perturbation from atmospheric background pressure already 
	corresponds to painfully loud sounds, on the level of a chainsaw by a listener's ear.

  .     .     .     .     .     .     .     .     .     .     .     .     .     .     .     .     .     
   `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` 



                  								|
                 							   /‾\
                							  |  :|`-._
                							  |  :|`-._`-._
               								 /   ::\   `-._`-._
              								/     ::\      `-(_)
             							   |_________|      / /
                						      \  ,/        / /
                 							   ‾‾         / /
                           								 / /
                          								/ /
         							   ________________/,_&___________
        							  /8P'      ________________ Y888/
       								 /P'      /   A guide to   / Y8/
      								/'  /\   /   the runtime  /   /
     							   /  . \ \ / initialization /   /
    						      /  //  \ \       process  /   /       |\      _,,,---,,_
   								 /  //    \ \______________/   /        /,`.-'`'    -.  ;-;;,_
  								/ ///      \_\        __      /        |,4-  ) )-,_..;\ (  `'-'
 							   /8 `'                 /_/    ./        '---''(_/--'  `-'\_)
							  /88b._______________________.8/


   `'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'
      `     `     `     `     `     `     `     `     `     `     `     `     `     `     `     `  

 ____         _      ____     
/\  _\      /' \    /\__ \    	   
\ \ \/     /\_, \   \/_/\ \       Upon running %% TimeEvoMain.m %% the user will first be asked to
 \ \ \     \/_/\ \     \ \ \       specify the number of dimensons to compute a solution in
  \ \ \_      \ \ \     \_\ \      (in regards to the spatial domain)
   \ \___\     \ \_\    /\___\
    \/___/      \/_/    \/___/

   Upon running %% TimeEvoMain.m %% the user will first be asked to specify the number of dimensons 
   to compute a solution in (in regards to the spatial domain)


	⚠ NOTE ⚠ As of 11/11/2025, only the 2D disc spatial domain is fully integrated into this
		  overhaul of the code that seeks to handle all cases from one main script, and
		  so instead of the previously mentioned script, until this is updated, run

		     %% \Euler time evolution\TwoD_Disk_TimeEvoMain.m %%

 .     .     .     .     .     .     .     .     .     .     .     .     .     .     .     .     .     
   `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` 
              
 ____          ___       ____      
/\  _\       /'___`\    /\__ \    
\ \ \/      /\_\ /\ \   \/_/\ \      The user is prompted for the number of distinct entropy levels
 \ \ \      \/_/// /__     \ \ \       to use for the piecewise constant entropy profile s(x)
  \ \ \_       // /_\ \     \_\ \ 
   \ \___\    /\______/     /\___\
    \/___/    \/_____/      \/___/

	Options: [nonisentropic] Choosing 2+ distinct entropy levels is valid for all spatial
                  dimensions. In this case, the user will be prompted to choose a background
                  temperature for the region corresponding to the first entropy level, which is 
                  thus set to s = 0. 
                  The user is then prompted to specify a temperature for each other distinct region
                  of the piecewise constant s(x). The corresponding value of s in that region is
                  then calculated with data tables for sea level air, relative to the first region.

		 [isentropic] Choosing that the entropy profile have just one distinct entropy
                  level, s = s_0, is only a valid option for the available 2 and 3 dimensional
                  spatial domains. These include the disc, the annulus, and the sphere.


	⚠ NOTE ⚠ Again, as of 11/11/2025, only the 2D disc spatial domain is fully integrated into
                   this overhaul of the code that seeks to handle all cases from one main script


  `'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'
      `     `     `     `     `     `     `     `     `     `     `     `     `     `     `     `  
                                
 ____         __       ____     
/\  _\      /'__`\    /\__ \       [For the disc/annulus or the sphere] The user is then asked to
\ \ \/     /\_\L\ \   \/_/\ \       indicate whether or not to use a uniform discretization of
 \ \ \     \/_/_\_<_     \ \ \      the radial interval [0,R]
  \ \ \_     /\ \L\ \     \_\ \ 
   \ \___\   \ \____/     /\___
    \/___/    \/___/      \/___/

	Options: [uniform] Input 1 to use an equally spaced discretization of [0,R]. Using this
		  option allows for the use of the usual Simpson's 1/3 rule, which we see as
		  striking a good balance between high accuracy and fast computation. However,
		  the choice is still offered as this does come with the downside of points
		  being more sparse as the radius increases.

		 [nonuniform] Input 0 to use a nonuniform discretization of [0,R]. If this is
		  chosen, then there is another choice to make between two different nonuniform
		  grid types.

			Sub options: [exponential]

				     [polynomial]

	⚠ NOTE ⚠ As of 11/28/2025, the polynomial nonuniform discretization is hard-coded in,
		option to use exponential discretization in writing, soon to be added to code.

  .     .     .     .     .     .     .     .     .     .     .     .     .     .     .     .     .     
   `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` 

 ____      __ __        ____     
/\  _\    /\ \\ \      /\__ \    
\ \ \/    \ \ \\ \     \/_/\ \      The user is prompted to specify the size of the eigenfunction
 \ \ \     \ \ \\ \_      \ \ \      bases to use (i.e. to specify the size of the finite subset of
  \ \ \_    \ \__ ,__\     \_\ \     the infinite dimensional eigenfunction bases to use)
   \ \___\   \/_/\_\_/     /\___\
    \/___/      \/_/       \/___/

	⚠ IMPORTANT ⚠ This is a very consequential choice, being one of the two sets of
		     parameters that have the greatest effect on the runtime length. The other
		     is, unsurprisingly, the number of points to use for the discretization(s)
		     of each variable's interval. (So in the case of the disc, this is to refer
		     to the N point discretization of [0,2pi] and the M point discretization 
		     of [0,R]. Currently, these values, M and N, are set within the code itself
		     and are not set by user input at runtime.
		     If you would like to change these, they are located at the top of the main
		     script for domain/entropy case (ex. disc-isentropic-main)

   
       /                      The author primarily ran this code on an Intel 13900K and
      /                       higher end calculations used approximately thetaBasisSize = 28,
     /____________________    rBasisSize = 20 on the disc, where the eigenfunctions are indexed
     |________  __________    by two variables, making the eigenfunction basis consist of
     /_____  /||   |          		
    |".___."| ||   |          		thetaBasisSize * rBasisSize + 1
    |_______|/ |   |          
     || .___."||  /           elements, with the +1 coming from the (0,0) eigenfunction.
     ||_______|| /            (keep in mind, for plotting and Hs norm analysis afterwards,
     |_________|/               only half the modes are used)

   With this basis size and a discretization of the disc using rPoints = 500, thetaPoints = 300,
   the computation can take up to 1 to 2 hours on similar hardware. However, this does depend on
   the next choice the user is prompted to make.

   (GPU version of the code being debugged/worked on at the moment, to be included later)

  `'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'
      `     `     `     `     `     `     `     `     `     `     `     `     `     `     `     `  
                                
 ____      ______      ____     
/\  _\    /\  ___\    /\__ \    
\ \ \/    \ \ \__/    \/_/\ \      Next to be specified is the particular element of our (finite)
 \ \ \     \ \___``\     \ \ \      eigenfunction basis to use as our fixed k-mode in the solution.
  \ \ \_    \/\ \L\ \     \_\ \ 
   \ \___\   \ \____/     /\___\
    \/___/    \/___/      \/___/

   These solutions are found by first setting our initial data to be the sum of some background
   pressure pBar and a chosen eigenfunction dubbed the k-mode:
	p^(0) (x,0) = pBar + alpha*phi_k
   Once the residual is computed, corrections are made only in the form of j-modes with j ~= k.
   (so j = 0 is acceptable, which would change the background pressure, we just want k fixed
    as we are trying to compute a nontrivial solution)



          ______________
         /             /|
        /             / |
       /____________ /  |
      |  _________  |   |____________________
      | |         | |   |/        /|,       /|
      | |     ..  | |   /        / /9      / |
      | |  .      | |  /_______ / /9      /  |
      | |_________| | |  ____ +| /9      /   |
      |________++___|/|________|/9      /    |
         ________________     ,9`      /   / |
        /  -/      /-   /|  ,9        /   /| |
       /______________ //|,9         /   / | |
      |       ______  ||,9          /   /  | |
      |  -+  |_9366_| ||/          /   /|  | |
      |_______________|/__________/   / |  | |
      /////----------/|           |  /__|  | |___
      |o     o  \o|  \|           |  |  |  | |
      |o    \|_  ||  o|______     |  |__|  | |_____
      |o \_  |   ||  o|      |    |  |  |  | /
      |o /   |\  /|  o|      |    |  |  |__|/
      |o             o|      |    |  |
      |o-------------o|      |    |  |
      |o   /\/\      o|      |    |  |
      |o  / o o|     o|      |    |  |
      |o / \_+_/     o|      |    |  |
      |o |\     \    o|      |    |  |
      |o | |+ +-|    o|      |    |  |
      |o-------------o|      |    |  |
      |o     /|      o|      |    | /   
       \/|/|/ |/\/|/\/       |____|/





 _.--.__.-'""`-.__.--.__.-'""`-.__.--.__.-'""`-.__.--.__.-'""`-._.--.__.-'""`-.__.--.__.-'""`-._
 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
_/""7__/""7__/""7__/""7__/""7__/""7__/""7__/""7__/""7__/""7__/""7__/""7__/""7__/""7__/""7__/""7_
 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 "`--'""`-.__.-'""`--'""`-.__.-'""`--'""`-.__.-'""`--'""`-.__.-'"`--'""`-.__.-'""`--'""`-.__.-'"



    _____________________________
   | \___===__________________()_\
   | |                            |   
   | |   _______________________  |   The problem's parameters are already set so as to best
   | |  |       iTunes         |  |    model air at sea level, at a temperature specified by
   | |  | > Generating and     |  |    the user within the range [50F,80F]
   | |  | > exporting sound    |  |
   | |  | > waves with         |  |   If we fix the position of a listener in the, then the
   | |  | > pressure data      |  |    changing pressure values over time at that spot should
   | |  | > from computed      |  |    allow us to create the sound heard by the the listener.
   | |  | > solutions          |  |
   | |  |                      |  |   For this reason, there is a dedicated script to isolate 
   | |  |______________________|  |    the varying values of pressure over time at a fixed
   | |                            |    position for each iterated solution.           
   | |             ##             |    
   | |         :........:         |   [which is needed as, if you are looking at the author's
   | |       --::......:---       |    data packs, they include all variables except for the
   | |      ----::....:---==      |    actual matrix of pressure values, which is too large
   | |     ====--      =+++**     |    to export in the default way. The pressure matrix
   | |     %%%%%%      %%%%%%     |    can be successfully exported though, if desired, by
   | |     **++==      :=++**     |    specifying -v7.3 when saving pNLMatrix]
   | |      =--:::....:--::-      |   
   | |       ::::......::--       |   Note that the generated audio can change a lot depending
   | |         :........:         |   on the position chosen to fix the "listener" at, this
   | |             ##             |   is still an area of ongoing research and coding work.
   | |                            |
    \|____________________________|


 ____         _      ____       
/\  _\      /' \    /\__ \      
\ \ \/     /\_, \   \/_/\ \       Assuming you are not looking to run an entire computation
 \ \ \     \/_/\ \     \ \ \       yourself and instead will be using precomputed data, start
  \ \ \_      \ \ \     \_\ \      by loading in a saved data pack within
   \ \___\     \ \_\    /\___\  
    \/___/      \/_/    \/___/  		%% \Euler time evolution\data packs %%


   `'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'
      `     `     `     `     `     `     `     `     `     `     `     `     `     `     `     `  

 ____          ___       ____     
/\  _\       /'___`\    /\__ \    
\ \ \/      /\_\ /\ \   \/_/\ \     
 \ \ \      \/_/// /__     \ \ \     After all the variables have been loaded into your workspace,
  \ \ \_       // /_\ \     \_\ \     run the script indicated below:
   \ \___\    /\______/     /\___\   	
    \/___/    \/_____/      \/___/

     	%% \Euler time evolution\analysis\TwoD_Disk_TimeEvo_CreateFixedPositionProfiles.m %%

  .     .     .     .     .     .     .     .     .     .     .     .     .     .     .     .     .     
   `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` 

                                
 ____         __       ____     
/\  _\      /'__`\    /\__ \    
\ \ \/     /\_\L\ \   \/_/\ \     The user will first be prompted to specify the location to fix
 \ \ \     \/_/_\_<_     \ \ \     the "listener" at, in terms of a theta-index
  \ \ \_     /\ \L\ \     \_\ \    (between 1 and thetaPoints) and an r-index (between 1 and rPoints)
   \ \___\   \ \____/     /\___\
    \/___/    \/___/      \/___/
                                

   `'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'
      `     `     `     `     `     `     `     `     `     `     `     `     `     `     `     `  

 ____      __ __        ____     
/\  _\    /\ \\ \      /\__ \    
\ \ \/    \ \ \\ \     \/_/\ \     The user is next prompted to specify which iteration step to
 \ \ \     \ \ \\ \_      \ \ \     go up to when reconstructing these fixed position pressure
  \ \ \_    \ \__ ,__\     \_\ \    profiles (we start at step 0, and stop at the indicated step j)
    \/___/      \/_/       \/___/

	This script creates the follow arrays:
		(a) pressureProfilesFixedFullPeriod(1:timePoints,1:chosenIterationStep)
		(b) pressureProfilesFixedFullPeriodNoBackground(1:timePoints,1:chosenIterationStep)
   	where the second array is, of course, the first but without the constant background pressure.
   
   Although the actual computations are done to the half period T (for 2T-periodic functions),
   these pressure profiles are extended to the time grid [0, dt, 2dt, ..., T-dt] so that
   the profiles are ready to be analyzed with something like the FFT

  .     .     .     .     .     .     .     .     .     .     .     .     .     .     .     .     .     
   `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` 

 ____      ______      ____      
/\  _\    /\  ___\    /\__ \      The script lastly prompts the user on whether or not to create
\ \ \/    \ \ \__/    \/_/\ \      extended versions of these profiles for sound generation. 
 \ \ \     \ \___``\     \ \ \     
  \ \ \_    \/\ \L\ \     \_\ \   (typically the period of the computed solutions seems to be 
   \ \___\   \ \____/     /\___\   measured in milliseconds, so to play as a sound many, many copies
    \/___/    \/___/      \/___/    of the profile need to be added on)


   If the user indicates that yes, they would like to make an extended copy of the profiles for 
   sound generation, the user will be prompted to input how many periods to make the extended
   versions. Good options seems to be 1024, 2048, or even 4096 periods, which can last anywhere
   from a second to a few seconds.


   `'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'
      `     `     `     `     `     `     `     `     `     `     `     `     `     `     `     `  

 ____        ____      ____      
/\  _\      /'___\    /\__ \     
\ \ \/     /\ \__/    \/_/\ \     To save the soundwave corresponding to a specific profile
 \ \ \     \ \  _``\     \ \ \     (i.e. for a specific iteration step), we usually input
  \ \ \_    \ \ \L\ \     \_\ \    directly in the command window:
   \ \___\   \ \____/     /\___\ 
    \/___/    \/___/      \/___/ 

      	audiowrite('Nov_15_2025_K23Alpha004_Step0Audio2048.wav', ...
		pressureProfileFixedExtendedResampled(:,1)/ ...
		max(abs(pressureProfileFixedExtendedResampled(:,1))),96000)
   
   The inputs into "audiowrite" are:
	(a) 'filename.wav'
	(b) pressureProfile(:,stepJ)/max(abs(pressureProfile(:,stepJ))
	    	this is the actual pressure profile, normalized so as to prevent clipping when
	    	writing the .wav file ("soundsc" handles scaling itself, but not "audiowrite")
	(c) 96000
		When writing an audio file, we want to use one of a number of typical audio 
		sample rates that most computers can handle. When the script is finished running,
		it will have already resampled the extended profiles to have sample rate 96000
		


	⚠ NOTE ⚠ If you are using my data packs, then you'll notice that in addition to the
		  usual data pack for any given date, which have names that look something like
		  Nov_16_2025_K14_Alpha004_LargeBases.mat, there is sometimes a separate .mat
		  file solely for the variable pNLMatrix, like with
		  Nov_16_2025_K14_Alpha004_LargeBases_pNLMatrix.mat

		  Sometimes I go to the trouble of saving this massive variable (pNLMatrix is
		  such that normally, it's file is larger than the normal data pack) for the
		  sake of thorough archives, but it is unnessecary (in this section) to either
		  load the .mat file for pNLMatrix or to run a script to reconstruct pNLMatrix


   `'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'
      `     `     `     `     `     `     `     `     `     `     `     `     `     `     `     `  


