/*
 * jADM1 -- Java Implementation of Anaerobic Digestion Model No 1
 * ===============================================================
 *
 * Copyright 2015 Liam Pettigrew
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *        http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * ********************************************************************************************
 */

package de.uni_erlangen.lstm.models.adm1;

import java.util.logging.Logger;

import de.uni_erlangen.lstm.file.CSVReader;
import de.uni_erlangen.lstm.file.CSVWriter;

/**
 * Specifies the fixed digester parameters
 * Modified from the BSM2 adjusted model adm1_ODE_bsm2.c for IAWQ AD Model No 1.
 * 
 * Special thanks to  Ulf Jeppsson, Christian Rosen and Darko Vrecko
 * for use of their Matlab code of the ADM1, 
 * developed when (around 2004) they were all working together at the 
 * Department of Industrial Electrical Engineering and Automation (IEA), Lund University, Sweden.
 * 
 * @author liampetti
 *
 */
public class DigesterParameters {
	public final static Logger LOGGER = Logger.getLogger(DigesterParameters.class.getName());
	
	/*
	 * Digestor Parameters
	 */
	private double f_sI_xc, f_xI_xc, f_ch_xc, f_pr_xc, f_li_xc, N_xc, N_I, N_aa, C_xc, C_sI, C_ch;
	private double C_pr, C_li, C_xI, C_su, C_aa, f_fa_li, C_fa, f_h2_su, f_bu_su, f_pro_su, f_ac_su;
	private double N_bac, C_bu, C_pro, C_ac, C_bac, Y_su, f_h2_aa, f_va_aa, f_bu_aa, f_pro_aa, f_ac_aa;
	private double C_va, Y_aa, Y_fa, Y_c4, Y_pro, C_ch4, Y_ac, Y_h2;
	private double k_dis, k_hyd_ch, k_hyd_pr, k_hyd_li, K_S_IN, k_m_su, K_S_su, pH_UL_aa, pH_LL_aa;
	private double k_m_aa, K_S_aa, k_m_fa, K_S_fa, K_Ih2_fa, k_m_c4, K_S_c4, K_Ih2_c4, k_m_pro, K_S_pro;
	private double K_Ih2_pro, k_m_ac, K_S_ac, K_I_nh3, pH_UL_ac, pH_LL_ac, k_m_h2, K_S_h2, pH_UL_h2, pH_LL_h2;
	private double k_dec_Xsu, k_dec_Xaa, k_dec_Xfa, k_dec_Xc4, k_dec_Xpro, k_dec_Xac, k_dec_Xh2;
	private double T_base, T_op, pK_w_base, pK_a_va_base, pK_a_bu_base, pK_a_pro_base, pK_a_ac_base, pK_a_co2_base, pK_a_IN_base;
	private double k_A_Bva, k_A_Bbu, k_A_Bpro, k_A_Bac, k_A_Bco2, k_A_BIN;
	private double k_P, kLa, K_H_h2o_base, K_H_co2_base, K_H_ch4_base, K_H_h2_base;
	private double V_liq, V_gas;
	
	/** 
	 * Default settings according to what you would typically find for sludge digesters
	 */
	public DigesterParameters () {
		/*
		 * Fixed Digester Parameters
		 */
		T_base=298.15;    		// 0.  Base temp: 25 degC = 298.15 K (273.15 + 25)
		T_op=308.15;      		// 1.  Operational temperature of AD and interfaces, 35 degC, should be an input - NOT USED
		pK_w_base=14.0; 		// 2.  Log10 of acid-base equilibrium coefficient
		pK_a_va_base=4.86;		// 3.  Log10 of acid-base equilibrium coefficient valerate
		pK_a_bu_base=4.82;		// 4.  Log10 of acid-base equilibrium coefficient butyrate
		pK_a_pro_base=4.88;		// 5.  Log10 of acid-base equilibrium coefficient propionate
		pK_a_ac_base=4.76;		// 6.  Log10 of acid-base equilibrium coefficient acetate
		pK_a_co2_base=6.35;		// 7.  Log10 of acid-base equilibrium coefficient carbon dioxide
		pK_a_IN_base=9.25;		// 8.  Log10 of acid-base equilibrium coefficient inorganic nitrogen
		K_H_h2_base=7.8e-4;		// 9.  Henry's law coefficient - hydrogen
		K_H_ch4_base=0.0014;	// 10. Henry's law coefficient - methane
		K_H_co2_base=0.035;		// 11. Henry's law coefficient - carbon dioxide
		K_H_h2o_base=0.0313; 	// 12. Henry's law coefficient - water
		pH_UL_aa=5.5;			// 13. Amino acids pH inhibition upper limit
		pH_LL_aa=4.0;			// 14. Amino acids pH inhibition lower limit
		pH_UL_ac=7.0;			// 15. Acetate pH inhibition upper limit
		pH_LL_ac=6.0;			// 16. Acetate pH inhibition lower limit
		pH_UL_h2=6.0;			// 17. Hydrogen pH inhibition upper limit
		pH_LL_h2=5.0;			// 18. Hydrogen pH inhibition lower limit
		K_S_IN=1.0e-4;			// 19. Half saturation value inorganic nitrogen
		K_Ih2_fa=5.0e-6;		// 20. 50% inhibitory concentration long chain fatty acids
		K_Ih2_c4=1.0e-5;		// 21. 50% inhibitory concentration valerate and butyrate
		K_Ih2_pro=3.5e-6;		// 22. 50% inhibitory concentration propionate
		K_I_nh3=0.0018;			// 23. Half saturation value ammonia
		k_dis=0.5;				// 24. First order disintegration rate
		k_hyd_ch=10.0;			// 25. First order hydrolysis parameter carbohydrates
		k_hyd_pr=10.0;			// 26. First order hydrolysis parameter proteins
		k_hyd_li=10.0;			// 27. First order hydrolysis parameter lipids
		k_m_su=30.0;			// 28. Monod maximum specific uptake rate monosaccharides
		K_S_su=0.5;				// 29. Half saturation value monosaccharides
		k_m_aa=50.0;			// 30. Monod maximum specific uptake rate amino acids
		K_S_aa=0.3;				// 31. Half saturation value amino acids
		k_m_fa=6.0;				// 32. Monod maximum specific uptake rate long chain fatty acids
		K_S_fa=0.4;				// 33. Half saturation value long chain fatty acids
		k_m_c4=20.0;			// 34. Monod maximum specific uptake rate valerate and butyrate
		K_S_c4=0.2;				// 35. Half saturation value valerate and butyrate
		k_m_pro=13.0;			// 36. Monod maximum specific uptake rate propionate
		K_S_pro=0.1;			// 37. Half saturation value propionate
		k_m_ac=8.0;				// 38. Monod maximum specific uptake rate acetate
		K_S_ac=0.15;			// 39. Half saturation value acetate
		k_m_h2=35.0;			// 40. Monod maximum specific uptake rate hydrogen
		K_S_h2=7.0e-6;			// 41. Half saturation value hydrogen
		k_dec_Xsu=0.02;			// 42. First order decay rate particulates monosaccharides
		k_dec_Xaa=0.02;			// 43. First order decay rate particulates amino acids
		k_dec_Xfa=0.02;			// 44. First order decay rate particulates long chain fatty acids
		k_dec_Xc4=0.02;			// 45. First order decay rate particulates valerate and butyrate
		k_dec_Xpro=0.02;		// 46. First order decay rate particulates propionate
		k_dec_Xac=0.02;			// 47. First order decay rate particulates acetate
		k_dec_Xh2=0.02;			// 48. First order decay rate insoluble hydrogen
		k_A_Bva=1.0e10;			// 49. Acid base kinetic parameter valerate
		k_A_Bbu=1.0e10;			// 50. Acid base kinetic parameter butyrate
		k_A_Bpro=1.0e10;		// 51. Acid base kinetic parameter propionate
		k_A_Bac=1.0e10;			// 52. Acid base kinetic parameter acetate
		k_A_Bco2=1.0e10;		// 53. Acid base kinetic parameter carbon dioxide
		k_A_BIN=1.0e10;			// 54. Acid base kinetic parameter inorganic nitrogen
		kLa=200.0;				// 55. Gas-liquid transfer coefficient
		C_xc=0.02786;			// 56. Carbon content of composite material
		f_sI_xc=0.1;			// 57. Yield (catabolism only) of soluble inerts on composite material
		C_sI=0.03;				// 58. Carbon content of soluble inerts
		f_ch_xc=0.2;			// 59. Yield (catabolism only) of carbohydrates on composite material
		C_ch=0.0313;			// 60. Carbon content of carbohydrates
		f_pr_xc=0.2;			// 61. Yield (catabolism only) of proteins on composite material
		C_pr=0.03;				// 62. Carbon content of proteins
		f_li_xc=0.3;			// 63. Yield (catabolism only) of lipids on composite material
		C_li=0.022;				// 64. Carbon content of lipids
		f_xI_xc=0.2;			// 65. Yield (catabolism only) of particulate inerts on composite material
		C_xI=0.03;				// 66. Carbon content of particulate inerts
		C_su=0.0313;			// 67. Carbon content of monosaccharides
		C_aa=0.03;				// 68. Carbon content of amino acids
		f_fa_li=0.95;			// 69. Yield (catabolism only) of long chain fatty acids on lipids
		C_fa=0.0217;			// 70. Carbon content of long chain fatty acids
		Y_su=0.1;				// 71. Yield of biomass on monosaccharides
		f_bu_su=0.13;			// 72. Yield (catabolism only) of butyrate on monosaccharides
		C_bu=0.025;				// 73. Carbon content of butyrate
		f_pro_su=0.27;			// 74. Yield (catabolism only) of propionate on monosaccharides
		C_pro=0.0268;			// 75. Carbon content of butyrate
		f_ac_su=0.41;			// 76. Yield (catabolism only) of acetate on monosaccharides
		C_ac=0.0313;			// 77. Carbon content of acetate
		C_bac=0.0313;			// 78. Carbon content of base acetate
		Y_aa=0.08;				// 79. Yield of biomass on amino acids
		f_va_aa=0.23;			// 80. Yield (catabolism only) of valerate on amino acids
		C_va=0.024;				// 81. Carbon content of valerate
		f_bu_aa=0.26;			// 82. Yield (catabolism only) of butyrate on amino acids
		f_pro_aa=0.05;			// 83. Yield (catabolism only) of propionate on amino acids
		f_ac_aa=0.40;			// 84. Yield (catabolism only) of acetate on amino acids
		Y_fa=0.06;				// 85. Yield of biomass on long chain fatty acids
		Y_c4=0.06;				// 86. Yield of biomass on valerate and butyrate
		Y_pro=0.04;				// 87. Yield of biomass on propionate
		Y_ac=0.05;				// 88. Yield of biomass on acetate
		C_ch4=0.0156;			// 89. Carbon content of methane
		Y_h2=0.06;				// 90. Yield of biomass on hydrogen
		f_h2_su=0.19;			// 91. Yield (catabolism only) of hydrogen on monosaccharides
		f_h2_aa=0.06;			// 92. Yield (catabolism only) of hydrogen on amino acids
		N_xc=0.0026857143;		// 93. Nitrogen component of composite material
		N_I=0.0042857143;		// 94. Nitrogen component of inorganic material
		N_aa=0.007;				// 95. Nitrogen component of amino acids
		N_bac=0.0057142857;		// 96. Nitrogen component of base acetate
		k_P=5.0e4;				// 97. Pipe resistance coefficient (m3 d-1 bar-1)
		V_liq=3400.0; 			// 98. Size of AD liquid portion in m3
		V_gas=300.0; 			// 99. Size of AD gas portion in m3
	}
	
	/**
	 * Read the parameters from a given CSV file
	 * 
	 * @param filename
	 */
	public void readParameters(String filename) {
		CSVReader reader = new CSVReader(filename, ";");
		String[] param = reader.getStrings();
		double[] p = new double[param.length];
		for (int i=0;i<p.length;i++) {
			p[i] = Double.parseDouble(param[i]);
		}
		setParameters(p);
	}
	
	/**
	 * Writes the current parameters to a CSV file
	 * 
	 * @param filename
	 */
	public void writeParameters(String filename) {
		double[] p = getParameters();
		CSVWriter writer = new CSVWriter();
		writer.WriteArray(filename, p);
	}
	
	/**
	 * Retrieves the parameters as an array
	 */
	public double[] getParameters() {
		return new double[] {
				T_base, T_op, pK_w_base,
				pK_a_va_base, pK_a_bu_base, pK_a_pro_base, pK_a_ac_base, pK_a_co2_base, pK_a_IN_base, K_H_h2_base,
				K_H_ch4_base, K_H_co2_base, K_H_h2o_base, pH_UL_aa, pH_LL_aa, pH_UL_ac, pH_LL_ac, pH_UL_h2,
				pH_LL_h2, K_S_IN, K_Ih2_fa, K_Ih2_c4, K_Ih2_pro, K_I_nh3, k_dis, k_hyd_ch, k_hyd_pr, k_hyd_li,
				k_m_su, K_S_su, k_m_aa, K_S_aa, k_m_fa, K_S_fa, k_m_c4, K_S_c4, k_m_pro, K_S_pro, k_m_ac, K_S_ac, 
				k_m_h2, K_S_h2, k_dec_Xsu, k_dec_Xaa, k_dec_Xfa, k_dec_Xc4, k_dec_Xpro, k_dec_Xac, k_dec_Xh2,
				k_A_Bva, k_A_Bbu, k_A_Bpro, k_A_Bac, k_A_Bco2, k_A_BIN, kLa, C_xc, f_sI_xc, C_sI, f_ch_xc, C_ch, 
				f_pr_xc, C_pr, f_li_xc, C_li, f_xI_xc, C_xI, C_su, C_aa, f_fa_li, C_fa, Y_su, f_bu_su, C_bu,
				f_pro_su, C_pro, f_ac_su, C_ac, C_bac, Y_aa, f_va_aa, C_va, f_bu_aa, f_pro_aa, f_ac_aa, Y_fa,
				Y_c4, Y_pro, Y_ac, C_ch4, Y_h2, f_h2_su, f_h2_aa, N_xc, N_I, N_aa, N_bac, k_P, V_liq, V_gas,
				};
	}
	
	/**
	 * Sets the parameters from an array
	 */
	public void setParameters(double[] param) {
		T_base=param[0];
		T_op=param[1];
		pK_w_base=param[2];
		pK_a_va_base=param[3];
		pK_a_bu_base=param[4];
		pK_a_pro_base=param[5];
		pK_a_ac_base=param[6];
		pK_a_co2_base=param[7];
		pK_a_IN_base=param[8];
		K_H_h2_base=param[9];
		K_H_ch4_base=param[10];
		K_H_co2_base=param[11];
		K_H_h2o_base=param[12];
		pH_UL_aa=param[13];
		pH_LL_aa=param[14];
		pH_UL_ac=param[15];
		pH_LL_ac=param[16];
		pH_UL_h2=param[17];
		pH_LL_h2=param[18];
		K_S_IN=param[19];
		K_Ih2_fa=param[20];
		K_Ih2_c4=param[21];
		K_Ih2_pro=param[22];
		K_I_nh3=param[23];
		k_dis=param[24];
		k_hyd_ch=param[25];
		k_hyd_pr=param[26];
		k_hyd_li=param[27];
		k_m_su=param[28];
		K_S_su=param[29];
		k_m_aa=param[30];
		K_S_aa=param[31];
		k_m_fa=param[32];
		K_S_fa=param[33];
		k_m_c4=param[34];
		K_S_c4=param[35];
		k_m_pro=param[36];
		K_S_pro=param[37];
		k_m_ac=param[38];
		K_S_ac=param[39];
		k_m_h2=param[40];
		K_S_h2=param[41];
		k_dec_Xsu=param[42];
		k_dec_Xaa=param[43];
		k_dec_Xfa=param[44];
		k_dec_Xc4=param[45];
		k_dec_Xpro=param[46];
		k_dec_Xac=param[47];
		k_dec_Xh2=param[48];
		k_A_Bva=param[49];
		k_A_Bbu=param[50];
		k_A_Bpro=param[51];
		k_A_Bac=param[52];
		k_A_Bco2=param[53];
		k_A_BIN=param[54];
		kLa=param[55];
		C_xc=param[56];
		f_sI_xc=param[57];
		C_sI=param[58];
		f_ch_xc=param[59];
		C_ch=param[60];
		f_pr_xc=param[61];
		C_pr=param[62];
		f_li_xc=param[63];
		C_li=param[64];
		f_xI_xc=param[65];
		C_xI=param[66];
		C_su=param[67];
		C_aa=param[68];
		f_fa_li=param[69];
		C_fa=param[70];
		Y_su=param[71];
		f_bu_su=param[72];
		C_bu=param[73];
		f_pro_su=param[74];
		C_pro=param[75];
		f_ac_su=param[76];
		C_ac=param[77];
		C_bac=param[78];
		Y_aa=param[79];
		f_va_aa=param[80];
		C_va=param[81];
		f_bu_aa=param[82];
		f_pro_aa=param[83];
		f_ac_aa=param[84];
		Y_fa=param[85];
		Y_c4=param[86];
		Y_pro=param[87];
		Y_ac=param[88];
		C_ch4=param[89];
		Y_h2=param[90];
		f_h2_su=param[91];
		f_h2_aa=param[92];
		N_xc=param[93];
		N_I=param[94];
		N_aa=param[95];
		N_bac=param[96];
		k_P=param[97];
		V_liq=param[98];
		V_gas=param[99];
	}
}
