  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.  
 / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ 
      `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   
| |											      | | 
			+----------------------------------------------+                      		      
| |			|  +---------------------------------------+   |		      | |	
			|  |  MATLAB R2025a                        |   |		      
| |			|  | ------------------------------------- |   |		      | |
			|  | +---+-------------------------------+ |   |		       		      
| |			|  | |~~~|    Computation of periodic    | |   |		      | | 
			|  | |~~~|    solutions for the 3 x 3    | |   |		       
| |			|  | |~~~|  compressible Euler equations | |   |		      | |  
			|  | |~~~|                               | |   |		        
| |			|  | |~~~| [part of the dissertation of  | |   |		      | |  
			|  | |~~~|   Andry Brinsko, University   | |   |		       
| |			|  | |~~~|   of Massachusetts Amherst]   | |   |		      | |  
			|  | +---+-------------------------------+ |   |		        
| |			|  | |   | © Andry Brinsko, 2023-2025    | |   |		      | |  
			|  | +---+-------------------------------+ |   |		       
| |			|  +---------------------------------------+   |		      | |
			+----------------------------------------------+		      
| |       	        |\  +------------------------------------------+ \		      | |
			\\  \_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\__\ \		      
| |	 		 \\  \___\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\____\ \		      | |
	  		  \\  \___\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\ ___  \ \		      
| |	   		   \\  \_____\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\ \_\ _\ \		      | |
	    		    \\  \_\_\_\_\_\_\_________________\___\_\_\_\_ \ \		      
| |	     		     \\  +-----------------------------------------+  \		      | |
  	             	      \\                                               \	      
| |               	       \\                                               \	      | |
				\+----------------------------------------------+|	      
| |		 		 ------------------------------------------------	      | |

| |											      | |
      .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.    
 \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / /
  `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-' 



 _.--.__.-'""`-.__.--.__.-'""`-.__.--.__.-'""`-.__.--.__.-'""`-._.--.__.-'""`-.__.--.__.-'""`-._
 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 "`--'""`-.__.-'""`--'""`-.__.-'""`--'""`-.__.-'""`--'""`-.__.-'"`--'""`-.__.-'""`--'""`-.__.-'"


      ._________________.       
      |.---------------.|       	
      ||   Running a   ||       	
      ||  computation  ||          To run a computation of a periodic solution of the Euler
      || of a solution ||          equations yourself, there is one main top-level script that
      ||  on your own  ||          allows for the user to specify some of the main parameters
      ||_______________||      	   of the problem at runtime:	
      /.-.-.-.-.-.-.-.-.\       		
     /.-.-.-.-.-.-.-.-.-.\      	     %% \Euler time evolution\TimeEvoMain.m %%
    /.-.-.-.-.-.-.-.-.-.-.\     	
   /______/__________\___o_\	
   \_______________________/    


1. Upon running %% TimeEvoMain.m %% the user will first be asked to specify the number of dimensons 
   to compute a solution in (in regards to the spatial domain)
	
	Options: 1, 2, and 3 dimensions


	⚠ NOTE ⚠ As of 11/11/2025, only the 2D disc spatial domain is fully integrated into this
		  overhaul of the code that seeks to handle all cases from one main script, and
		  so instead of the previously mentioned script, until this is updated, run

		     %% \Euler time evolution\TwoD_Disk_TimeEvoMain.m %%


   `'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'
      `     `     `     `     `     `     `     `     `     `     `     `     `     `     `     `  

2. The user is prompted for the number of distinct entropy levels to use for the piecewise constant
   entropy profile s(x)

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

  .     .     .     .     .     .     .     .     .     .     .     .     .     .     .     .     .     
   `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` 


3. [For the disc/annulus or the sphere] The user is then asked to indicate whether or not to use a
   uniform discretization of the radial interval [0,R]

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


   `'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'
      `     `     `     `     `     `     `     `     `     `     `     `     `     `     `     `  

4. The user is prompted to specify the size of the eigenfunction bases to use (i.e. to specify
   the size of the finite subset of the infinite dimensional eigenfunction bases to use)

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

  .     .     .     .     .     .     .     .     .     .     .     .     .     .     .     .     .     
   `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` 


5. Next to be specified is the particular element of our (finite) eigenfunction basis to use as
   our fixed k-mode in the solution.

   These solutions are found by first setting our initial data to be the sum of some background
   pressure pBar and a chosen eigenfunction dubbed the k-mode:
	p^(0) (x,0) = pBar + alpha*phi_k
   Once the residual is computed, corrections are made only in the form of j-modes with j ~= k.
   (so j = 0 is acceptable, which would change the background pressure, we just want k fixed
    as we are trying to compute a nontrivial solution)




 _.--.__.-'""`-.__.--.__.-'""`-.__.--.__.-'""`-.__.--.__.-'""`-._.--.__.-'""`-.__.--.__.-'""`-._
 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 "`--'""`-.__.-'""`--'""`-.__.-'""`--'""`-.__.-'""`--'""`-.__.-'"`--'""`-.__.-'""`--'""`-.__.-'"



    _____________________________
   | \___===__________________()_\
   | |                            |   
   | |   _______________________  |   The problem's parameters are already set so as to best
   | |  |                      |  |   model air at sea level, at a temperature specified by
   | |  | > Generating and     |  |   the user within the range [50F,80F]
   | |  | > exporting sound    |  |
   | |  | > waves with         |  |   If we fix the position of a listener in the, then the
   | |  | > pressure data      |  |   changing pressure values over time at that spot should
   | |  | > from computed      |  |   allow us to create the sound heard by the the listener.
   | |  | > solutions          |  |
   | |  |                      |  |   For this reason, there is a dedicated script to isolate
   | |  |                      |  |   the varying values of pressure over time at a fixed
   | |  |                      |  |   position for each iterated solution.
   | |  |______________________|  |
   | |                            |   [which is needed as, if you are looking at the author's
   | |             ##             |    data packs, they include all variables except for the
   | |         :........:         |    actual matrix of pressure values, which is too large
   | |       --::......:---       |    to export in the default way. The pressure matrix
   | |      ----::....:---==      |    can be successfully exported though, if desired, by
   | |     ====--      =+++**     |    specifying -v7.3 when saving pNLMatrix]
   | |     %%%%%%      %%%%%%     |
   | |     **++==      :=++**     |   
   | |      =--:::....:--::-      |   
   | |       ::::......::--       |   Note that the generated audio can change a lot depending
   | |         :........:         |   on the position chosen to fix the "listener" at, this
   | |             ##             |   is still an area of ongoing research and coding work.
   | |                            |
    \|____________________________|


1. Assuming you are not looking to run an entire computation yourself and instead will be using
   precomputed data, start by loading in a saved data pack within

	%\Euler time evolution\data packs


   `'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'
      `     `     `     `     `     `     `     `     `     `     `     `     `     `     `     `  

2. After all the variables have been loaded into your workspace, run the following script

     %% \Euler time evolution\analysis\TwoD_Disk_TimeEvo_CreateFixedPositionProfiles.m %%

  .     .     .     .     .     .     .     .     .     .     .     .     .     .     .     .     .     
   `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` 


3. The user will first be prompted to specify the location to fix the "listener" at, in terms
   of a theta-index (between 1 and thetaPoints) and an r-index (between 1 and rPoints)


   `'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'
      `     `     `     `     `     `     `     `     `     `     `     `     `     `     `     `  

4. The user is next prompted to specify which iteration step to go up to when reconstructing
   these fixed position pressure profiles (we start at step 0, and stop at the indicated step j)

  .     .     .     .     .     .     .     .     .     .     .     .     .     .     .     .     .     
   `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` 


5. This script creates the follow arrays:
	(a) pressureProfilesFixedFullPeriod(1:timePoints,1:chosenIterationStep)
	(b) pressureProfilesFixedFullPeriodNoBackground(1:timePoints,1:chosenIterationStep)
   where the second array is, of course, the first but without the constant background pressure.
   
   Although the actual computations are done to the half period T (for 2T-periodic functions),
   these pressure profiles are extended to the time grid [0, dt, 2dt, ..., T-dt] so that
   the profiles are ready to be analyzed with something like the FFT


   `'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'`'. .'
      `     `     `     `     `     `     `     `     `     `     `     `     `     `     `     `  

6. The script lastly prompts the user on whether or not to create extended versions of these 
   profiles for the purposes of sound generation. 
   (typically the period of the computed solutions seems to be measured in milliseconds, so to
   play as a sound many, many copies of the profile need to be added on)

   If the user indicates that yes, they would like to make an extended copy of the profiles for 
   sound generation, the user will be prompted to input how many periods to make the extended
   versions. Good options seems to be 1024, 2048, or even 4096 periods, which can last anywhere
   from a second to a few seconds.

  .     .     .     .     .     .     .     .     .     .     .     .     .     .     .     .     .     
   `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` `._.` 


7. To save the soundwave corresponding to a specific profile (i.e. for a specific iteration step),
   we usually input directly in the command window:

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
