package de.uni_erlangen.lstm.models.adm1;

/**
 * Specifies the default digester initial and influent state variables
 * Modified from the BSM2 adjusted model for IAWQ AD Model No 1.
 * 
 * Special thanks to  Ulf Jeppsson, Christian Rosen and Darko Vrecko
 * for use of their Matlab code of the ADM1, 
 * developed when (around 2004) they were all working together at the 
 * Department of Industrial Electrical Engineering and Automation (IEA), Lund University, Sweden.
 * 
 * @author liampetti
 *
 */
public class BSM2Defaults {
	
	/*
	 * Dynamic State Variables
	 */
	private double S_su, S_aa, S_fa, S_va, S_bu, S_pro, S_ac, S_h2, S_ch4, S_IC, S_IN, S_I;
	private double X_xc, X_ch, X_pr, X_li, X_su, X_aa, X_fa, X_c4, X_pro, X_ac, X_h2, X_I;
	private double S_cat, S_an, S_hva, S_hbu, S_hpro, S_hac, S_hco3, S_nh3, S_gas_h2, S_gas_ch4, S_gas_co2;
	private double Q_D, T_D, gas_ch4, gas_vol, ph, temp1, temp2;
	
	public double[] DigesterInit() {
		/*
		 * Digester Initialisation - Outputs (42 Variables)
		 */
		// array position, (ADM1 id), description, (units)
		S_su =  0.012; // 0.	(1) monosaccharides (kg COD/m3)[g COD/L]
		S_aa =  0.0053; // 1. 	(2) amino acids (kg COD/m3)[g COD/L]
		S_fa =  0.099; // 2. 	(3) total long chain fatty acids (LCFA) (kg COD/m3)[g COD/L]
		S_va =  0.012; // 3. 	(4) total valerate (kg COD/m3)[g COD/L]
		S_bu =  0.013; // 4. 	(5) total butyrate (kg COD/m3)[g COD/L]
		S_pro = 0.016; // 5. 	(6) total propionate (kg COD/m3)[g COD/L]
		S_ac =  0.2; // 6. 		(7) total acetate (kg COD/m3)[g COD/L]
		S_h2 =  2.30e-7; // 7. 	(8) hydrogen gas (kg COD/m3)[g COD/L]
		S_ch4 = 0.055; // 8. 	(9) methane gas (kg COD/m3)[g COD/L]
		S_IC = 0.15; // 9. 		(10) inorganic carbon (kmole C/m3)
		S_IN = 0.13; // 10. 	(11) inorganic nitrogen (kmole N/m3)
		S_I = 0.033; // 11. 	(12) soluble inerts (kg COD/m3)[g COD/L]
		// Particulates ---->
		X_xc = 0.31; // 12. 	(13) composites (kg COD/m3)[g COD/L]
		X_ch = 0.028; // 13. 	(14) carbohydrates (kg COD/m3)[g COD/L]
		X_pr = 0.1; // 14. 		(15) proteins (kg COD/m3)[g COD/L]
		X_li = 0.029; // 15. 	(16) lipids (kg COD/m3)	[g COD/L]	
		X_su = 0.42; // 16. 	(17) sugar degraders (kg COD/m3)[g COD/L]
		X_aa = 1.18; // 17. 	(18) amino acid degraders (kg COD/m3)[g COD/L]
		X_fa = 0.24; // 18.  	(19) LCFA degraders (kg COD/m3)[g COD/L]
		X_c4 = 0.43; // 19. 	(20) valerate and butyrate degraders (kg COD/m3)[g COD/L]
		X_pro = 0.14; // 20. 	(21) propionate degraders (kg COD/m3)[g COD/L]
		X_ac = 0.76; // 21. 	(22) acetate degraders (kg COD/m3)[g COD/L]
		X_h2 =  0.32; // 22. 	(23) hydrogen degraders (kg COD/m3)[g COD/L]
		X_I =  25.6; // 23. 	(24) particulate inerts (kg COD/m3)[g COD/L]
		// <---- Particulates
		S_cat =  0.04; // 24. cations (metallic ions, strong base) (kmole/m3)
		S_an =  0.02; // 25. anions (metallic ions, strong acid) (kmole/m3)
		S_hva = 0.011;   // 26. is actually Sva- = valerate (kg COD/m3)[g COD/L]
		S_hbu = 0.013;   // 27. is actually Sbu- = butyrate (kg COD/m3)[g COD/L]
		S_hpro =  0.016; // 28. is actually Spro- = propionate (kg COD/m3)[g COD/L]
		S_hac =  0.2;   // 29. is actually Sac- = acetate (kg COD/m3)[g COD/L]
		S_hco3 = 0.14; // 30. bicarbonate (kmole C/m3)
		S_nh3 = 0.0041; // 31. ammonia (kmole N/m3)
		S_gas_h2 = 1.02e-005; // 32. hydrogen concentration in gas phase (kg COD/m3)[g COD/L]
		S_gas_ch4 = 1.63; // 33. methane concentration in gas phase (kg COD/m3)[g COD/L]
		S_gas_co2 = 0.014; // 34. carbon dioxide concentration in gas phase (kmole C/m3)
		Q_D = 0.0; // 35. flow rate (m3/d) - Set by influent
		T_D = 35.0; // 36. digester temperature (deg C)
		gas_ch4 = 0.0; // 37. methane volume (m3/d)
		gas_vol = 0.0; // 38. gas volume (m3/d)
		ph = 0.0; // 39. pH
		temp1 = 0.0; //40. Optional
		temp2 = 0.0; // 41. Optional
		
		return new double[] { S_su, S_aa, S_fa, S_va, S_bu, S_pro, S_ac, S_h2, S_ch4,
				S_IC, S_IN, S_I, X_xc, X_ch, X_pr, X_li, X_su, X_aa, X_fa, X_c4, X_pro, X_ac,
				X_h2, X_I, S_cat, S_an, S_hva, S_hbu, S_hpro, S_hac, S_hco3, S_nh3, S_gas_h2, S_gas_ch4,
				S_gas_co2, Q_D, T_D, gas_ch4, gas_vol, ph, temp1, temp2 };
	}
	
	public double[] Influent() {
		/*
		 * Digester Influent (42 Variables)
		 */
		// array position, (ADM1 id), description, (units)
		S_su =  0.01; // 0.		(1) monosaccharides (kg COD/m3)[g COD/L]
		S_aa =  0.001; // 1. 	(2) amino acids (kg COD/m3)[g COD/L]
		S_fa =  0.001; // 2. 	(3) total long chain fatty acids (LCFA) (kg COD/m3)[g COD/L]
		S_va =  0.001; // 3. 	(4) total valerate (kg COD/m3)[g COD/L]
		S_bu =  0.001; // 4. 	(5) total butyrate (kg COD/m3)[g COD/L]
		S_pro = 0.001; // 5. 	(6) total propionate (kg COD/m3)[g COD/L]
		S_ac =  0.001; // 6. 	(7) total acetate (kg COD/m3)[g COD/L]
		S_h2 =  1.0e-8; // 7. 	(8) hydrogen gas (kg COD/m3)[g COD/L]
		S_ch4 = 1.0e-5; // 8. 	(9) methane gas (kg COD/m3)[g COD/L]
		S_IC = 0.04; // 9. 		(10) inorganic carbon (kmole C/m3)
		S_IN = 0.01; // 10. 	(11) inorganic nitrogen (kmole N/m3)
		S_I = 0.02; // 11. 		(12) soluble inerts (kg COD/m3)[g COD/L]
		// Particulates ---->
		X_xc = 2.0; // 12. 		(13) composites (kg COD/m3)[g COD/L]
		X_ch = 5.0; // 13. 		(14) carbohydrates (kg COD/m3)[g COD/L]
		X_pr = 20.0; // 14. 	(15) proteins (kg COD/m3)[g COD/L]
		X_li = 5.0; // 15. 		(16) lipids (kg COD/m3)	[g COD/L]	
		X_su = 0.0; // 16. 		(17) sugar degraders (kg COD/m3)[g COD/L]
		X_aa = 0.01; // 17. 	(18) amino acid degraders (kg COD/m3)[g COD/L]
		X_fa = 0.01; // 18.  	(19) LCFA degraders (kg COD/m3)[g COD/L]
		X_c4 = 0.01; // 19. 	(20) valerate and butyrate degraders (kg COD/m3)[g COD/L]
		X_pro = 0.01; // 20. 	(21) propionate degraders (kg COD/m3)[g COD/L]
		X_ac = 0.01; // 21. 	(22) acetate degraders (kg COD/m3)[g COD/L]
		X_h2 =  0.01; // 22. 	(23) hydrogen degraders (kg COD/m3)[g COD/L]
		X_I =  25.0;; // 23. 	(24) particulate inerts (kg COD/m3)[g COD/L]
		// <---- Particulates
		S_cat =  0.04; // 24. cations (metallic ions, strong base) (kmole/m3)
		S_an =  0.02; // 25. anions (metallic ions, strong acid) (kmole/m3)
		// Set from init --->
		S_hva = 0.0;   // 26. is actually Sva- = valerate (kg COD/m3)[g COD/L]
		S_hbu = 0.0;   // 27. is actually Sbu- = butyrate (kg COD/m3)[g COD/L]
		S_hpro =  0.0; // 28. is actually Spro- = propionate (kg COD/m3)[g COD/L]
		S_hac =  0.0;   // 29. is actually Sac- = acetate (kg COD/m3)[g COD/L]
		S_hco3 = 0.0; // 30. bicarbonate (kmole C/m3)
		S_nh3 = 0.0; // 31. ammonia (kmole N/m3)
		// <--- Set from init
		S_gas_h2 = 0.0; // 32. hydrogen concentration in gas phase (kg COD/m3)[g COD/L]
		S_gas_ch4 = 0.0; // 33. methane concentration in gas phase (kg COD/m3)[g COD/L]
		S_gas_co2 = 0.0; // 34. carbon dioxide concentration in gas phase (kmole C/m3)		
		Q_D = 170.0; // 35. flow rate (m3/d)
		T_D = 0.0; // 36. temperature (deg C) - SET BY DIGESTER
		gas_ch4 = 0.0; // 37. Optional
		gas_vol = 0.0; // 38. gas volume (m3/d)
		ph = 0.0; // 39. pH
		temp1 = 0.0; // 40. Optional
		temp2 = 0.0; // 41. Optional
		
		return new double[] { S_su, S_aa, S_fa, S_va, S_bu, S_pro, S_ac, S_h2, S_ch4,
				S_IC, S_IN, S_I, X_xc, X_ch, X_pr, X_li, X_su, X_aa, X_fa, X_c4, X_pro, X_ac,
				X_h2, X_I, S_cat, S_an, S_hva, S_hbu, S_hpro, S_hac, S_hco3, S_nh3, S_gas_h2, S_gas_ch4,
				S_gas_co2, Q_D, T_D, gas_ch4, gas_vol, ph, temp1, temp2 };
	}

}
