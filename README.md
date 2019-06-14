# Collective-oscillations

This is the distributed code for the manuscript:  Shou-Wen Wang & Lei-Han Tang, Emergence of collective oscillations in adaptive cells, https://arxiv.org/abs/1809.06997


## Structure:
1) **Final model for glycolsyis**: mathematica-based code for glycolysis related models;				

   a) **glycolysis_dupreez2 _unperturb_19_6_7**:  the original dupreez2 model, no perturbation is done here
   
   b) **glycolysis_dupreez2 _perturbed_19_6_7**:  the modified dupreeze2 model.  Response of the system to stepwise or periodic change of ACE       is computed and a phase diagram of the whole system is generated. In particular, it includes the following two models
      
         i) glycolysisACE: perturbation of the original model by changing the ACE level
      
         ii) glycolysisACEremoveGLYCO: perturbation of the model without the GLYO reaction
      
   c) **glycoACEpert_Phase_diagram**: data generated by glycolysis_dupreez2 _perturbed_19_6_7 for making the phase diagram
   
   d) **minimum_Model_glycolysis_19_6_7**: the minimum model we constructed based on the response features of the full model.  
   
         i) only the model for the intracellular circuit
      
         ii) intracellular circuit + ACE diffusion to allow cell-to-cell communication, assuming fast signal diffusion
      
         iii) intracellular circuit + ACE diffusion to allow cell-to-cell communication, explicitly treat signal inside and outside the cell,            with the diffusion speed controlled by the diffusion constant
  
2) **collective oscillations**: matlab-based code for all other models and calculations in this manuscript

	a) **startup**: a script that add the files to the searching path of matlab

   b) **adaptation**: related to the adaptive circuit illustrated in the main text
   
       i) Main_fre_amp_prediction_comparison_adaptation: numerically compute the oscillation amplitude and frequency from a given adaptive          model and compare it to the prediction based from the analytic approximation of the adaptive spectrum
   
       ii) main_predicting_epsilon_tau_star: compute the relationship between the adaptation error epsilon and the minimum signal timescale          tau that is needed to induce collective oscillations (whatever the cell density is needed)
      
       iii) paper_fig_plot_adaptation_oscillation: explore various functions, compute the numerically exact response spectrum etc.
       
       iV) functions: includes all sorts functions
       
   c) **FitzHugh-Nagumo Model**
   
       i) Main_fre_amp_prediction_comparison: compute the amplitude and frequency of the collective oscillations, and compare it to that            predicted from the numerically computed response spectrum
         
       ii) paper_fig_plot_FNH_oscillation: explore various functions, compute the numerically exact response spectrum etc.
       
   d) **My_Common_Function**: common functions for plot setting, and computation of the response spectrum or correlation spectrum etc. 
     
   
