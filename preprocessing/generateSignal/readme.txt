This code generates simulated events to the MiniAOD level using the model defined in ./models. I currently only have this set up to produce 2018 events.

To set this up for another user I think all that should need to happen is to change all of the hardcoded paths. The easiest way to do this in my opinion would be to copy this whole directory somewhere on your afs space and then use some code editor (I use Visual Studio Code) to search for all instances of "atownse2". If you want to copy my directory mappings you could just replace "atownse2" with your username and set up the same directories.

- Temporary directory : XRootD can't write directly to Hadoop so the files need to be written somewhere else first. I use the scratch365 space for this.
- Gridpack directory : Gridpacks make event generation more scalable, to generate signal first we need to make the gridpacks. I store these on Hadoop.
- Output directory : I use Hadoop for this as well
- Test directory : Optional but I find it to be nice to set up a test directory somewhere with a CMSSW instance so I can make sure the event generation is doing ok.
- Condor directories: For condor output, they are specified in generateSignal.py and by default are placed in the same directory as the script.

I think you will need to modify paths in:
- run_event_generation.sh
- generateSignal.py
- genproductions/bin/MadGraph5_aMCatNLO/gridpack_generation.sh

To run this code you only need to use generateSignal.py, modifying the arguments for your use case. This generate a gridpack (interactively) if it doesn't exist and then submits jobs to condor to generate and process the events. An example would be:

python generateSignal.py --m_moe 500 0.1 -n 10 (-b)

This will generate 10 events with a Bkk mass of 500 and Radion rest mass/energy of 0.1. The "-b" argument says to submit to condor. I usually run a test run of 10 events interactively (i.e. without submitting to condor) just to make sure that everything is running correctly. There are other ways to specify what you want to run.

I have a mass grid defined elsewhere that I use as a shorthand to generate a bunch of mass points at once, you can comment out this section (L7-9,40-41).

You might have some issues with python, I use python3 and I don't think the default version on earth comes with the htcondor package installed. The easiest thing to do might be to just reformat all the f strings so that you can use the earth default which is python2. If you want to use python3 I would suggest either using a distribution from CMSSW_12+ or setting up a mamba/conda environment with htcondor installed (though you also need to copy the configs to get it to work right with the mamba/conda route). If you go the mamba/conda route let me know and I can share how I got it to work.
