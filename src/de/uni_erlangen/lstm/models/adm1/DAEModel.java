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

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;


/**
 * Modified from the BSM2 adjusted model for IAWQ AD Model No 1.
 * adm1_ODE.c, adm1_DAE1.c, adm1_DAE2.c, pHsolv.c, Sh2solv.c
 * 
 * Special thanks to  Ulf Jeppsson, Christian Rosen and Darko Vrecko
 * for use of their Matlab code of the ADM1, 
 * developed when (around 2004) they were all working together at the 
 * Department of Industrial Electrical Engineering and Automation (IEA), Lund University, Sweden.
 * 
 * @author liampetti
 *
 */
public class DAEModel implements FirstOrderDifferentialEquations {
	public final static Logger LOGGER = Logger.getLogger(DAEModel.class.getName());
	
	private boolean shDAE;
	private boolean sh2DAE;

	private double eps, phi, S_H_ion;
	private double proc1, proc2, proc3, proc4, proc5, proc6, proc7, proc8, proc9, proc10, proc11, proc12, proc13;
	private double proc14, proc15, proc16, proc17, proc18, proc19;
	private double procT8, procT9, procT10;
	private double I_pH_aa, I_pH_ac, I_pH_h2, I_IN_lim, I_h2_fa, I_h2_c4, I_h2_pro, I_nh3;
	private double reac1, reac2, reac3, reac4, reac5, reac6, reac7, reac8, reac9, reac10, reac11, reac12, reac13;
	private double reac14, reac15, reac16, reac17, reac18, reac19, reac20, reac21, reac22, reac23, reac24;
	private double stoich1, stoich2, stoich3, stoich4, stoich5, stoich6, stoich7, stoich8, stoich9, stoich10, stoich11, stoich12, stoich13;
	private double p_gas_h2o, P_gas, p_gas_h2, p_gas_ch4, p_gas_co2, q_gas;
	private double pHLim_aa, pHLim_ac, pHLim_h2, n_aa, n_ac, n_h2;
	private double K_w, K_a_va, K_a_bu, K_a_pro, K_a_ac, K_a_co2, K_a_IN, K_H_co2, K_H_ch4, K_H_h2;

	private double[] inhib;
	private double[] param;
	private double[] u; // influent
	private double[] xtemp;
	private double factor, R, P_atm;
	private double fix_pH;
	
	/** 
	 * Initiates the model using the defined parameters and pre-calculates the stoichiometry parameter values for use in the water phase
	 * 
	 * @param influent The influent
	 * @param initial The initial digester outputs
	 * @param parameters The digester parameters
	 * @param sh Initial S_H_ion value
	 * @param dae Turn on or off the dae system
	 */
	public DAEModel(double[] influent, double[] parameters, double sh, boolean dae, double ph) {		
		shDAE = dae; // Turn on algebraic equations for SH+ (ph and ion states)
		sh2DAE = dae; // Turn on algebraic equations for SH2
		S_H_ion = sh;
		fix_pH = ph;
		inhib = new double[6]; // holds inhibition functions
		//u = influent; // Influent pointer
		u = new double[influent.length];		
		for (int i=0;i<influent.length;i++) {
			u[i] = influent[i];
		}
		xtemp = new double[u.length]; // holds temporary zero values
		param = parameters; // digester parameters
		eps = 0.000001; // Small constant in case of poor choice of initial conditions for proc8,9
		P_atm = 1.013;	// bar
		R = 0.083145;	// universal gas constant dm3*bar/(mol*K) = 8.3145 J/(mol*K)

		// Stoichiometry for use in water phase equations
		// stoich1 = -C_xc+f_sI_xc*C_sI+f_ch_xc*C_ch+f_pr_xc*C_pr+f_li_xc*C_li+f_xI_xc*C_xI
		stoich1 = -param[56]+param[57]*param[58]+param[59]*param[60]+param[61]*param[62]+param[63]*param[64]+param[65]*param[66];
		// stoich2 = -C_ch+C_su
		stoich2 = -param[60]+param[67];
		// stoich3 = -C_pr+C_aa
		stoich3 = -param[62]+param[68];
		// stoich 4 = -C_li+(1.0-f_fa_li)*C_su+f_fa_li*C_fa
		stoich4 = -param[64]+(1.0-param[69])*param[67]+param[69]*param[70];
		// stoich5 = -C_su+(1.0-Y_su)*(f_bu_su*C_bu+f_pro_su*C_pro+f_ac_su*C_ac)+Y_su*C_bac
		stoich5 = -param[67]+(1.0-param[71])*(param[72]*param[73]+param[74]*param[75]+param[76]*param[77])+param[71]*param[78];
		// stoich6 = -C_aa+(1.0-Y_aa)*(f_va_aa*C_va+f_bu_aa*C_bu+f_pro_aa*C_pro+f_ac_aa*C_ac)+Y_aa*C_bac
		stoich6 = -param[68]+(1.0-param[79])*(param[80]*param[81]+param[82]*param[73]+param[83]*param[75]+param[84]*param[77])+param[79]*param[78];
		// stoich7 = -C_fa+(1.0-Y_fa)*0.7*C_ac+Y_fa*C_bac
		stoich7 = -param[70]+(1.0-param[85])*0.7*param[77]+param[85]*param[78];
		// stoich8 = -C_va+(1.0-Y_c4)*0.54*C_pro+(1.0-Y_c4)*0.31*C_ac+Y_c4*C_bac
		stoich8 = -param[81]+(1.0-param[86])*0.54*param[75]+(1.0-param[86])*0.31*param[77]+param[86]*param[78];
		// stoich9 = -C_bu+(1.0-Y_c4)*0.8*C_ac+Y_c4*C_bac
		stoich9 = -param[73]+(1.0-param[86])*0.8*param[77]+param[86]*param[78];
		// stoich10 = -C_pro+(1.0-Y_pro)*0.57*C_ac+Y_pro*C_bac
		stoich10 = -param[75]+(1.0-param[87])*0.57*param[77]+param[87]*param[78];
		// stoich11 = -C_ac+(1.0-Y_ac)*C_ch4+Y_ac*C_bac
		stoich11 = -param[77]+(1.0-param[88])*param[89]+param[88]*param[78];
		// stoich12 = (1.0-Y_h2)*C_ch4+Y_h2*C_bac
		stoich12 = (1.0-param[90])*param[89]+param[90]*param[78];
		// stoich13 = -C_bac+C_xc
		stoich13 = -param[78]+param[56];
		
		// pH Inhibition
		pHLim_aa = Math.pow(10,(-(param[13] + param[14])/2.0));
		pHLim_ac = Math.pow(10,(-(param[15] + param[16])/2.0));
		pHLim_h2 = Math.pow(10,(-(param[17] + param[18])/2.0));
		n_aa = 3.0/(param[13]-param[14]);
		n_ac = 3.0/(param[15]-param[16]);
		n_h2 = 3.0/(param[17]-param[18]);
	}
	
	// Function for retrieving the current variables from the model
	public double[] getDimensions() {		
		return xtemp;
	}
	
	@Override
	public void computeDerivatives(double t, double[] x, double[] dx)
			throws MaxCountExceededException, DimensionMismatchException {	
		for (int i=0;i<x.length;i++) {
			if (x[i]<0) {
				xtemp[i] = 0.0;
			} else {
				xtemp[i] = x[i];
			}
		}
		
		// Adjustments for acid-base equations
		factor = (1.0/(param[0]) - 1.0/(273.15+xtemp[36]))/(100.0*R);
		K_w = Math.pow(10,-param[2])*Math.exp(55900.0*factor); // T adjustment for K_w 
		K_a_co2 = Math.pow(10,-param[7])*Math.exp(7646.0*factor); // T adjustment for K_a_co2 
		K_a_IN = Math.pow(10,-param[8])*Math.exp(51965.0*factor); // T adjustment for K_a_IN 		
		K_H_h2 = param[9]*Math.exp(-4180.0*factor);     // T adjustment for K_H_h2
		K_H_ch4 = param[10]*Math.exp(-14240.0*factor);  // T adjustment for K_H_ch4
		K_H_co2 = param[11]*Math.exp(-19410.0*factor);  // T adjustment for K_H_co2
		p_gas_h2o = param[12]*Math.exp(5290.0*(1.0/(param[0]) - 1.0/(273.15+xtemp[36])));  // T adjustment for water vapour saturation pressure	
			
		K_a_va = Math.pow(10,-param[3]);
		K_a_bu = Math.pow(10,-param[4]);
		K_a_pro = Math.pow(10,-param[5]);
		K_a_ac = Math.pow(10,-param[6]);
		
		if (fix_pH >= 0) {
			// S_H_ion based on set pH
			shDAE = false;
			S_H_ion = Math.pow(10, -fix_pH);
			
			// Run the DAE functions
			runDAE();
		} else {
			// Run the DAE functions
			runDAE();
			
			// SH+ Equation (pH and ion states)
			if (!shDAE) {
				// Scat+(S_IN-Snh3)-hco3-(Sac/64)-(Spro/112)-(Sbu/160)-(Sva/208)-San
				phi = xtemp[24]+(xtemp[10]-xtemp[31])-xtemp[30]-(xtemp[29]/64.0)-(xtemp[28]/112.0)-
						(xtemp[27]/160.0)-(xtemp[26]/208.0)-xtemp[25];
				S_H_ion = (-phi*0.5)+0.5*Math.sqrt(phi*phi+(4.0*K_w)); // SH+
			} 
		}
		
		// Adjustments for gas pressure
		p_gas_h2 = xtemp[32]*R*(273.15+xtemp[36])/16.0;
		p_gas_ch4 = xtemp[33]*R*(273.15+xtemp[36])/64.0;
		p_gas_co2 = xtemp[34]*R*(273.15+xtemp[36]);
		P_gas = p_gas_h2 + p_gas_ch4 + p_gas_co2 + p_gas_h2o;
				
		// pH Inhibition
		I_pH_aa = Math.pow(pHLim_aa,n_aa)/(Math.pow(S_H_ion,n_aa)+Math.pow(pHLim_aa ,n_aa));
		I_pH_ac = Math.pow(pHLim_ac,n_ac)/(Math.pow(S_H_ion,n_ac)+Math.pow(pHLim_ac ,n_ac));
		I_pH_h2 = Math.pow(pHLim_h2,n_h2)/(Math.pow(S_H_ion,n_h2)+Math.pow(pHLim_h2 ,n_h2));
		
		I_IN_lim = 1.0/(1.0+param[19]/xtemp[10]); // 1.0/(1.0+K_S_IN/S_IN)
		I_h2_fa = 1.0/(1.0+xtemp[7]/param[20]); // 1.0/(1.0+S_h2/K_Ih2_fa)
		I_h2_c4 = 1.0/(1.0+xtemp[7]/param[21]); // 1.0/(1.0+S_h2/K_Ih2_c4)
		I_h2_pro = 1.0/(1.0+xtemp[7]/param[22]); // 1.0/(1.0+S_h2/K_Ih2_pro)
		I_nh3 = 1.0/(1.0+xtemp[31]/param[23]); // 1.0/(1.0+S_nh3/K_I_nh3)
		
		// Inhibitors
		inhib[0] = I_pH_aa*I_IN_lim; // Inhibition Equation 5 & 6
		inhib[1] = inhib[0]*I_h2_fa; // Inhibition Equation 7
		inhib[2] = inhib[0]*I_h2_c4; // Inhibition Equation 8 & 9
		inhib[3] = inhib[0]*I_h2_pro; // Inhibition Equation 10
		inhib[4] = I_pH_ac*I_IN_lim*I_nh3; // Inhibition Equation 11
		inhib[5] = I_pH_h2*I_IN_lim; // Inhibition Equation 12	
		
		// Biochemical process rates
		proc1 = param[24]*xtemp[12]; // k_dis*X_xc, Disintegration
		proc2 = param[25]*xtemp[13]; // k_hyd_ch*X_ch, Hydrolysis of carbohydrates
		proc3 = param[26]*xtemp[14]; // k_hyd_pr*X_pr, Hydrolysis of proteins
		proc4 = param[27]*xtemp[15]; // k_hyd_li*X_li, Hydrolysis of lipids
		proc5 = param[28]*xtemp[0]/(param[29]+xtemp[0])*xtemp[16]*inhib[0]; // k_m_su*(S_su/(K_S_su+S_su))*X_su*inhib_5, Uptake of sugars
		proc6 = param[30]*xtemp[1]/(param[31]+xtemp[1])*xtemp[17]*inhib[0]; // k_m_aa*(S_aa/(K_S_aa+S_aa))*X_aa*inhib_6, Uptake of amino acids
		proc7 = param[32]*xtemp[2]/(param[33]+xtemp[2])*xtemp[18]*inhib[1]; // k_m_fa*(S_fa/(K_S_fa+S_fa))*X_aa*inhib_7, Uptake of LCFA
		proc8 = param[34]*xtemp[3]/(param[35]+xtemp[3])*xtemp[19]*xtemp[3]/(xtemp[3]+xtemp[4]+eps)*inhib[2]; // k_m_c4*(S_va/(K_S_c4+S_va))*X_c4*(S_va/(S_bu+S_va+eps))*inhib_8, Uptake of valerate
		proc9 = param[34]*xtemp[4]/(param[35]+xtemp[4])*xtemp[19]*xtemp[4]/(xtemp[3]+xtemp[4]+eps)*inhib[2]; // k_m_c4*(S_bu/(K_S_c4+S_bu))*X_c4*(S_bu/(S_va+S_bu+eps))*inhib_9, Uptake of butyrate
		proc10 = param[36]*xtemp[5]/(param[37]+xtemp[5])*xtemp[20]*inhib[3]; // k_m_pro*(S_pro/(K_S_pro+S_pro))*X_pro*inhib_10, Uptake of propionate
		proc11 = param[38]*xtemp[6]/(param[39]+xtemp[6])*xtemp[21]*inhib[4]; // k_m_ac*(S_ac/(K_S_ac+S_ac))*X_ac*inhib_11, Uptake of acetate
		proc12 = param[40]*xtemp[7]/(param[41]+xtemp[7])*xtemp[22]*inhib[5]; // k_m_h2*(S_h2/(K_S_h2+S_h2))*X_h2*inhib_12, Uptake of hydrogen
		proc13 = param[42]*xtemp[16]; // k_dec_Xsu*X_su, Decay of X_su
		proc14 = param[43]*xtemp[17]; // k_dec_Xaa*X_aa, Decay of X_aa
		proc15 = param[44]*xtemp[18]; // k_dec_Xfa*X_fa, Decay of X_fa
		proc16 = param[45]*xtemp[19]; // k_dec_Xc4*X_c4, Decay of X_c4
		proc17 = param[46]*xtemp[20]; // k_dec_Xpro*X_pro, Decay of X_pro
		proc18 = param[47]*xtemp[21]; // k_dec_Xac*X_ac, Decay of X_ac
		proc19 = param[48]*xtemp[22]; // k_dec_Xh2*X_h2, Decay of X_h2
		
		// Gas transfer rates
		procT8 = param[55]*(xtemp[7]-16.0*K_H_h2*p_gas_h2); // kLa*(S_h2-16.0*K_H_h2*p_gas_h2)
		procT9 = param[55]*(xtemp[8]-64.0*K_H_ch4*p_gas_ch4); // kLa*(S_ch4-64.0*K_H_ch4*p_gas_ch4)
		procT10 = param[55]*((xtemp[9]-xtemp[30])-K_H_co2*p_gas_co2); // kLa*((S_IC-S_hco3)-K_H_co2*p_gas_co2)
		
		// Reactions
		// reac1 = proc2+(1.0-f_fa_li)*proc4-proc5;
		reac1 = proc2+(1.0-param[69])*proc4-proc5;
		reac2 = proc3-proc6;
		// reac3 = f_fa_li*proc4-proc7;
		reac3 = param[69]*proc4-proc7;
		// reac4 = (1.0-Y_aa)*f_va_aa*proc6-proc8;
		reac4 = (1.0-param[79])*param[80]*proc6-proc8;
		//reac5 = (1.0-Y_su)*f_bu_su*proc5+(1.0-Y_aa)*f_bu_aa*proc6-proc9;
		reac5 = (1.0-param[71])*param[72]*proc5+(1.0-param[79])*param[82]*proc6-proc9;
		// reac6 = (1.0-Y_su)*f_pro_su*proc5+(1.0-Y_aa)*f_pro_aa*proc6+(1.0-Y_c4)*0.54*proc8-proc10;
		reac6 = (1.0-param[71])*param[74]*proc5+(1.0-param[79])*param[83]*proc6+(1.0-param[86])*0.54*proc8-proc10;
		// reac7 = (1.0-Y_su)*f_ac_su*proc5+(1.0-Y_aa)*f_ac_aa*proc6+(1.0-Y_fa)*0.7*proc7+(1.0-Y_c4)*0.31*proc8+(1.0-Y_c4)*0.8*proc9+(1.0-Y_pro)*0.57*proc10-proc11;
		reac7 = (1.0-param[71])*param[76]*proc5+(1.0-param[79])*param[84]*proc6+(1.0-param[85])*0.7*proc7+(1.0-param[86])*0.31*proc8+(1.0-param[86])*0.8*proc9+(1.0-param[87])*0.57*proc10-proc11;
		// reac8 = (1.0-Y_su)*f_h2_su*proc5+(1.0-Y_aa)*f_h2_aa*proc6+(1.0-Y_fa)*0.3*proc7+(1.0-Y_c4)*0.15*proc8+(1.0-Y_c4)*0.2*proc9+(1.0-Y_pro)*0.43*proc10-proc12-procT8;
		reac8 = (1.0-param[71])*param[91]*proc5+(1.0-param[79])*param[92]*proc6+(1.0-param[85])*0.3*proc7+(1.0-param[86])*0.15*proc8+(1.0-param[86])*0.2*proc9+(1.0-param[87])*0.43*proc10-proc12-procT8;
		// reac9 = (1.0-Y_ac)*proc11+(1.0-Y_h2)*proc12-procT9;
		reac9 = (1.0-param[88])*proc11+(1.0-param[90])*proc12-procT9;		
		reac10 = -stoich1*proc1-stoich2*proc2-stoich3*proc3-stoich4*proc4-stoich5*proc5-stoich6*proc6-stoich7*proc7-stoich8*proc8-stoich9*proc9-stoich10*proc10-stoich11*proc11-stoich12*proc12-stoich13*proc13-stoich13*proc14-stoich13*proc15-stoich13*proc16-stoich13*proc17-stoich13*proc18-stoich13*proc19-procT10;
		// reac11 = (N_xc-f_xI_xc*N_I-f_sI_xc*N_I-f_pr_xc*N_aa)*proc1-Y_su*N_bac*proc5+(N_aa-Y_aa*N_bac)*proc6-Y_fa*N_bac*proc7-Y_c4*N_bac*proc8-Y_c4*N_bac*proc9-Y_pro*N_bac*proc10-Y_ac*N_bac*proc11-Y_h2*N_bac*proc12+(N_bac-N_xc)*(proc13+proc14+proc15+proc16+proc17+proc18+proc19);
		reac11 = -(param[93]-param[65]*param[94]-param[57]*param[94]-param[61]*param[95])*proc1-param[71]*param[96]*proc5+(param[95]-param[79]*param[96])*proc6-param[85]*param[96]*proc7-param[86]*param[96]*proc8-param[86]*param[96]*proc9-param[87]*param[96]*proc10-param[88]*param[96]*proc11-param[90]*param[96]*proc12+(param[96]-param[93])*(proc13+proc14+proc15+proc16+proc17+proc18+proc19);
		// reac12 = f_sI_xc*proc1;
		reac12 = param[57]*proc1;
		reac13 = -proc1+proc13+proc14+proc15+proc16+proc17+proc18+proc19;
		// reac14 = f_ch_xc*proc1-proc2;
		reac14 = param[59]*proc1-proc2;
		// reac15 = f_pr_xc*proc1-proc3;
		reac15 = param[61]*proc1-proc3;
		// reac16 = f_li_xc*proc1-proc4;
		reac16 = param[63]*proc1-proc4;
		// reac17 = Y_su*proc5-proc13;
		reac17 = param[71]*proc5-proc13;
		// reac18 = Y_aa*proc6-proc14;
		reac18 = param[79]*proc6-proc14;
		// reac19 = Y_fa*proc7-proc15;
		reac19 = param[85]*proc7-proc15;
		// reac20 = Y_c4*proc8+Y_c4*proc9-proc16;
		reac20 = param[86]*proc8+param[86]*proc9-proc16;
		// reac21 = Y_pro*proc10-proc17;
		reac21 = param[87]*proc10-proc17;
		// reac22 = Y_ac*proc11-proc18;
		reac22 = param[88]*proc11-proc18;
		// reac23 = Y_h2*proc12-proc19;
		reac23 = param[90]*proc12-proc19;
		// reac24 = f_xI_xc*proc1;
		reac24 = param[65]*proc1;
		
		q_gas = param[97]*(P_gas-P_atm);
		if (q_gas < 0)
		   q_gas = 0.0;
			   
		// DE's -> Soluble matter
		// dSsu/dt = Qad/Vad,liq(Ssu,i-Ssu)+reac1
		dx[0] = (x[35]/param[98])*(u[0]-x[0])+reac1; // Ssu
		dx[1] = (x[35]/param[98])*(u[1]-x[1])+reac2; // Saa
		dx[2] = (x[35]/param[98])*(u[2]-x[2])+reac3; // Sfa
		dx[3] = (x[35]/param[98])*(u[3]-x[3])+reac4; // Sva
		dx[4] = (x[35]/param[98])*(u[4]-x[4])+reac5; // Sbu
		dx[5] = (x[35]/param[98])*(u[5]-x[5])+reac6; // Spro
		dx[6] = (x[35]/param[98])*(u[6]-x[6])+reac7; // Sac

		if (!sh2DAE) {	
			dx[7] = (x[35]/param[98])*(u[7]-x[7])+reac8; // Sh2
		} 
				
		dx[8] = (x[35]/param[98])*(u[8]-x[8])+reac9; // Sch4
		dx[9] = (x[35]/param[98])*(u[9]-x[9])+reac10;    // SIC
		dx[10] = (x[35]/param[98])*(u[10]-x[10])+reac11; // SIN
		dx[11] = (x[35]/param[98])*(u[11]-x[11])+reac12; // SI
		
		// DE's -> Particulate matter
		dx[12] = (x[35]/param[98])*(u[12]-x[12])+reac13; // Xc
		dx[13] = (x[35]/param[98])*(u[13]-x[13])+reac14; // Xch
		dx[14] = (x[35]/param[98])*(u[14]-x[14])+reac15; // Xpr
		dx[15] = (x[35]/param[98])*(u[15]-x[15])+reac16; // Xli
		dx[16] = (x[35]/param[98])*(u[16]-x[16])+reac17; // Xsu
		dx[17] = (x[35]/param[98])*(u[17]-x[17])+reac18; // Xaa
		dx[18] = (x[35]/param[98])*(u[18]-x[18])+reac19; // Xfa
		dx[19] = (x[35]/param[98])*(u[19]-x[19])+reac20; // Xc4
		dx[20] = (x[35]/param[98])*(u[20]-x[20])+reac21; // Xpro
		dx[21] = (x[35]/param[98])*(u[21]-x[21])+reac22; // Xac
		dx[22] = (x[35]/param[98])*(u[22]-x[22])+reac23; // Xh2
		dx[23] = (x[35]/param[98])*(u[23]-x[23])+reac24; // XI

		dx[24] = (x[35]/param[98])*(u[24]-x[24]); // Scat+
		dx[25] = (x[35]/param[98])*(u[25]-x[25]); // San-
		
		// Acid-base process rates for ODE
		//k_A_Bva*(S_hva*(K_A_va+S_H_ion)-K_a_va*S_va)
		if (!shDAE) {
			dx[26] = -(param[49]*(xtemp[26]*(K_a_va+S_H_ion)-K_a_va*xtemp[3]));  // Sva-
			dx[27] = -(param[50]*(xtemp[27]*(K_a_bu+S_H_ion)-K_a_bu*xtemp[4]));  // Sbu-
			dx[28] = -(param[51]*(xtemp[28]*(K_a_pro+S_H_ion)-K_a_pro*xtemp[5]));  // Spro-
			dx[29] = -(param[52]*(xtemp[29]*(K_a_ac+S_H_ion)-K_a_ac*xtemp[6]));  // Sac-
			dx[30] = -(param[53]*(xtemp[30]*(K_a_co2+S_H_ion)-K_a_co2*xtemp[9])); // SHCO3-
			dx[31] = -(param[54]*(xtemp[31]*(K_a_IN+S_H_ion)-K_a_IN*xtemp[10])); // SNH3	
		} else {
			xtemp[26] = K_a_va*x[3]/(K_a_va+S_H_ion);
			xtemp[27] = K_a_bu*x[4]/(K_a_bu+S_H_ion);
			xtemp[28] = K_a_pro*x[5]/(K_a_pro+S_H_ion);
			xtemp[29] = K_a_ac*x[6]/(K_a_ac+S_H_ion);
			xtemp[30] = K_a_co2*x[9]/(K_a_co2+S_H_ion);
			xtemp[31] = K_a_IN*x[10]/(K_a_IN+S_H_ion);
		}

		dx[32] = -xtemp[32]*q_gas/param[99]+procT8*param[98]/param[99]; // Sgas,h2
		dx[33] = -xtemp[33]*q_gas/param[99]+procT9*param[98]/param[99]; // Sgas,ch4
		dx[34] = -xtemp[34]*q_gas/param[99]+procT10*param[98]/param[99]; // Sgas,co2

		dx[35] = 0; // Flow
		dx[36] = 0; // Temp
		
		// Correction by factor of 64.0 due to COD basis of Sgas,ch4  // Methane gas (m3/d)
		//xtemp[37] = (q_gas*xtemp[33])*R*(273.15+xtemp[36])/64.0; // Calculate methane flow from concentration in gas phase
		xtemp[37] = q_gas*(p_gas_ch4/P_gas); // Calculate methane flow from partial pressures
			
		xtemp[38] = q_gas;// Gas production (m3/d)
		
		
		xtemp[39] = -Math.log10(S_H_ion); // pH
		
		xtemp[40] = 0; 
		xtemp[41] = 0; 
	}
	
	public void runDAE() {			
		double prevS_H_ion = S_H_ion;
		
		// Stiffness below
		double shDelta = 1.0;
		double shGradEqu = 1.0;
		double sh2Delta = 1.0;
		double sh2GradEqu = 1.0;

		// Newton-Raphson Variables
		double TOL = 1e-12;
		double maxSteps = 1000;
		int i = 1;
		int j = 1;
		
		// SH+ Equation (pH and ion states)
		if (shDAE) {
			while ( (shDelta > TOL || shDelta < -TOL) && (i <= maxSteps) ) {				
				shDelta = xtemp[24]+xtemp[10]
						-(K_a_IN*xtemp[10]/(K_a_IN+S_H_ion))+S_H_ion
						-(K_a_co2*xtemp[9]/(K_a_co2+S_H_ion))
						-(K_a_ac*xtemp[6]/(K_a_ac+S_H_ion))/64.0
						-(K_a_pro*xtemp[5]/(K_a_pro+S_H_ion))/112.0						
						-(K_a_bu*xtemp[4]/(K_a_bu+S_H_ion))/160.0						
						-(K_a_va*xtemp[3]/(K_a_va+S_H_ion))/208.0						
						-(K_w/S_H_ion)-xtemp[25];
				
				shGradEqu = 1+K_a_IN*xtemp[10]/((K_a_IN+S_H_ion)*(K_a_IN+S_H_ion))
			            +K_a_co2*xtemp[9]/((K_a_co2+S_H_ion)*(K_a_co2+S_H_ion))          
			            +1/64.0*K_a_ac*xtemp[6]/((K_a_ac+S_H_ion)*(K_a_ac+S_H_ion))
			            +1/112.0*K_a_pro*xtemp[5]/((K_a_pro+S_H_ion)*(K_a_pro+S_H_ion))
			            +1/160.0*K_a_bu*xtemp[4]/((K_a_bu+S_H_ion)*(K_a_bu+S_H_ion))
			            +1/208.0*K_a_va*xtemp[3]/((K_a_va+S_H_ion)*(K_a_va+S_H_ion))
			            +K_w/(S_H_ion*S_H_ion);
				
				S_H_ion = S_H_ion - shDelta/shGradEqu;
				
				if (S_H_ion <= 0) {
		            S_H_ion = TOL;
		        }
				i++;
			}			
		}

		// SH2 Equation
		if (sh2DAE) {
			while ( (sh2Delta > TOL || sh2Delta < -TOL) && (j <= maxSteps) ) {
				// Calculate ahead within loop	
				I_pH_aa = Math.pow(pHLim_aa,n_aa)/(Math.pow(prevS_H_ion,n_aa)+Math.pow(pHLim_aa ,n_aa));
				I_pH_h2 = Math.pow(pHLim_h2,n_h2)/(Math.pow(prevS_H_ion,n_h2)+Math.pow(pHLim_h2 ,n_h2));
				
				I_IN_lim = 1.0/(1.0+param[19]/xtemp[10]); // 1.0/(1.0+K_S_IN/S_IN)
				I_h2_fa = 1.0/(1.0+xtemp[7]/param[20]); // 1.0/(1.0+S_h2/K_Ih2_fa)
				I_h2_c4 = 1.0/(1.0+xtemp[7]/param[21]); // 1.0/(1.0+S_h2/K_Ih2_c4)
				I_h2_pro = 1.0/(1.0+xtemp[7]/param[22]); // 1.0/(1.0+S_h2/K_Ih2_pro)
				
				// Inhibitors
				inhib[0] = I_pH_aa*I_IN_lim; // Inhibition Equation 5 & 6
				inhib[1] = inhib[0]*I_h2_fa; // Inhibition Equation 7
				inhib[2] = inhib[0]*I_h2_c4; // Inhibition Equation 8 & 9
				inhib[3] = inhib[0]*I_h2_pro; // Inhibition Equation 10
				inhib[5] = I_pH_h2*I_IN_lim; // Inhibition Equation 12	
				
				proc5 = param[28]*xtemp[0]/(param[29]+xtemp[0])*xtemp[16]*inhib[0]; // k_m_su*(S_su/(K_S_su+S_su))*X_su*inhib_5, Uptake of sugars
				proc6 = param[30]*xtemp[1]/(param[31]+xtemp[1])*xtemp[17]*inhib[0]; // k_m_aa*(S_aa/(K_S_aa+S_aa))*X_aa*inhib_6, Uptake of amino acids
				proc7 = param[32]*xtemp[2]/(param[33]+xtemp[2])*xtemp[18]*inhib[1]; // k_m_fa*(S_fa/(K_S_fa+S_fa))*X_aa*inhib_7, Uptake of LCFA
				proc8 = param[34]*xtemp[3]/(param[35]+xtemp[3])*xtemp[19]*xtemp[3]/(xtemp[3]+xtemp[4]+eps)*inhib[2]; // k_m_c4*(S_va/(K_S_c4+S_va))*X_c4*(S_va/(S_bu+S_va+eps))*inhib_8, Uptake of valerate
				proc9 = param[34]*xtemp[4]/(param[35]+xtemp[4])*xtemp[19]*xtemp[4]/(xtemp[3]+xtemp[4]+eps)*inhib[2]; // k_m_c4*(S_bu/(K_S_c4+S_bu))*X_c4*(S_bu/(S_va+S_bu+eps))*inhib_9, Uptake of butyrate
				proc10 = param[36]*xtemp[5]/(param[37]+xtemp[5])*xtemp[20]*inhib[3]; // k_m_pro*(S_pro/(K_S_pro+S_pro))*X_pro*inhib_10, Uptake of propionate
				
				proc12 = param[40]*xtemp[7]/(param[41]+xtemp[7])*xtemp[22]*inhib[5]; // k_m_h2*(S_h2/(K_S_h2+S_h2))*X_h2*inhib_12, Uptake of hydrogen
					
				p_gas_h2 = xtemp[32]*R*(273.15+xtemp[36])/16.0;
				procT8 = param[55]*(xtemp[7]-16.0*K_H_h2*p_gas_h2); // kLa*(S_h2-16.0*K_H_h2*p_gas_h2)
				
				reac8 = (1.0-param[71])*param[91]*proc5+(1.0-param[79])*param[92]*proc6+(1.0-param[85])*0.3*proc7+(1.0-param[86])*0.15*proc8+(1.0-param[86])*0.2*proc9+(1.0-param[87])*0.43*proc10-proc12-procT8;
				
				sh2Delta = (xtemp[35]/param[98])*(u[7]-xtemp[7])+reac8;
			               //-1/V_liq**u[26]
				sh2GradEqu = -1/param[98]*xtemp[35]
					  //-3.0/10.0*(1-Y_fa)*k_m_fa**u[2]/(K_S_fa+*u[2])**u[18]*I_pH_aa/(1+K_S_IN/(*u[10]))/((1+x[0]/K_Ih2_fa)*(1+x[0]/K_Ih2_fa))/K_Ih2_fa
						-3.0/10.0*(1-param[85])*param[32]*xtemp[2]/(param[33]+xtemp[2])*xtemp[18]*I_pH_aa/(1+param[19]/(xtemp[10]))/((1+xtemp[7]/param[20])*(1+xtemp[7]/param[20]))/param[20]
					  //-3.0/20.0*(1-Y_c4)*k_m_c4**u[3]**u[3]/(K_S_c4+*u[3])**u[19]/(*u[4]+*u[3]+eps)*I_pH_aa/(1+K_S_IN/(*u[10]))/((1+x[0]/K_Ih2_c4)*(1+x[0]/K_Ih2_c4))/K_Ih2_c4    
						-3.0/20.0*(1-param[86])*param[34]*xtemp[3]*xtemp[3]/(param[35]+xtemp[3])*xtemp[19]/(xtemp[4]+xtemp[3]+eps)*I_pH_aa/(1+param[19]/(xtemp[10]))/((1+xtemp[7]/param[21])*(1+xtemp[7]/param[21]))/param[21]
			          //-1.0/5.0*(1-Y_c4)*k_m_c4**u[4]**u[4]/(K_S_c4+*u[4])**u[19]/(*u[4]+*u[3]+eps)*I_pH_aa/(1+K_S_IN/(*u[10]))/((1+x[0]/K_Ih2_c4)*(1+x[0]/K_Ih2_c4))/K_Ih2_c4  
						-1.0/5.0*(1-param[86])*param[34]*xtemp[4]*xtemp[4]/(param[35]+xtemp[4])*xtemp[19]/(xtemp[4]+xtemp[3]+eps)*I_pH_aa/(1+param[19]/(xtemp[10]))/((1+xtemp[7]/param[21])*(1+xtemp[7]/param[21]))/param[21]
			          //-43.0/100.0*(1-Y_pro)*k_m_pro**u[5]/(K_S_pro+*u[5])**u[20]*I_pH_aa/(1+K_S_IN/(*u[10]))/((1+x[0]/K_Ih2_pro)*(1+x[0]/K_Ih2_pro))/K_Ih2_pro
						-43.0/100.0*(1-param[87])*param[36]*xtemp[5]/(param[37]+xtemp[5])*xtemp[20]*I_pH_aa/(1+param[19]/(xtemp[10]))/((1+xtemp[7]/param[22])*(1+xtemp[7]/param[22]))/param[22]
			          //-k_m_h2/(K_S_h2+x[0])**u[22]*I_pH_h2/(1+K_S_IN/(*u[10]))+k_m_h2*x[0]/((K_S_h2+x[0])*(K_S_h2+x[0]))**u[22]*I_pH_h2/(1+K_S_IN/(*u[10]))
						-param[40]/(param[41]+xtemp[7])*xtemp[22]*I_pH_h2/(1+param[19]/(xtemp[10]))+param[40]*xtemp[7]/((param[41]+xtemp[7])*(param[41]+xtemp[7]))*xtemp[22]*I_pH_h2/(1+param[19]/(xtemp[10]))
			          //-kLa;
						-param[55];
				
				xtemp[7] = xtemp[7]-sh2Delta/sh2GradEqu;
				
				if (xtemp[7] <= 0) {
		            xtemp[7] = TOL;
		        }
				
				j++;
			}
		}
	}

	@Override
	public int getDimension() {
		return 42;
	}
}
