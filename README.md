# FP-SI-code-analysis_Georgiou et al. 
The MATLAB code was used for analysis of the Fiber photometry data during social interaction
We made use of code from Martianova et al. (2019), _J Vis Exp_ and Zhang et al. (2010), _Analyst_
Details of the analyses and use of core are provided in the Fiber Photometry methods section in the manuscript
The get_zdff function is not directly called by the code but included for completeness - all its component code is reproduced in the master code. That code is broken into sections: signal processing , then sorting the keystrokes (with d / f referring to the 2 types of interactions) and finding the interaction epochs, then running analysis first not as timelocked data, and second as timelocked data. Although multiple baseline methods are outputted for the timelocked analysis we only used the data from the method indicated (3B) for subsequent analysis
