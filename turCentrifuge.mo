package turbineCentrifuge
  import Cons = Modelica.Constants;
  import SI = Modelica.SIunits;
  import Math = Modelica.Math;

  model turbine
    extends turbinePartielle;
    parameter SI.Length diamEntree "B10";
    parameter SI.Length diamSortie "B12";
    parameter SI.Area epaisseurEntree "B11";
    parameter SI.Angle angleEntree(displayUnit = "°") "B13";
    parameter SI.Angle angleSortie(displayUnit = "°") "B14";
    // etat du fluide
    Medium.ThermodynamicState etat1 = Medium.setState_phX(port_a.p, inStream(port_a.h_outflow), inStream(port_a.Xi_outflow)) "Equation etat entree turbine";
    Medium.ThermodynamicState etat2 = Medium.setState_pTX(port_b.p, T3, inStream(port_a.Xi_outflow)) "Equation etat sortie turbine";
    // variable (Etage 1 : Inducer, Etage 2:Entree rotor, Etage 3:sortie rotor)
    SI.AngularVelocity N(displayUnit = "rpm");
    //SI.Temperature T1 = Medium.temperature(etat1);
    //SI.Temperature T2;
    SI.Temperature T3;
    SI.Temperature T02, T03;
    SI.Pressure P01, P03;
    SI.Enthalpy h01, h03;
    SI.Velocity C1;
    SI.Velocity C2, U2, W2;
    SI.Velocity C3, U3, W3;
    SI.Enthalpy W;
    SI.Power Pmeca;
    //SI.VolumeFlowRate Q;
    Medium.SpecificHeatCapacity Cp_e = Medium.specificHeatCapacityCp(etat1);
    Medium.SpecificHeatCapacity Cv_e = Medium.specificHeatCapacityCv(etat1);
    Medium.SpecificHeatCapacity Cp_s = Medium.specificHeatCapacityCp(etat2);
    Medium.Density rho_e = Medium.density(etat1);
    Medium.Density rho_s = Medium.density(etat2);
    // Variable du compresseur
  equation
    C2 = mt_flow / rho_e / (2 * Cons.pi * diamEntree / 2 * epaisseurEntree) / Math.cos(SI.Conversions.from_deg(angleEntree));
    C1 = C2 * Math.cos(SI.Conversions.from_deg(angleEntree));
    U2 = C2 * Math.sin(SI.Conversions.from_deg(angleEntree));
    N = SI.Conversions.to_rpm(U2 / (diamEntree / 2));
    U3 = SI.Conversions.from_rpm(N) * diamSortie / 2;
    W2 = U2 / Math.tan(SI.Conversions.from_deg(angleEntree));
    C3 = U3 / Math.tan(SI.Conversions.from_deg(angleSortie));
    W3 = U3 / Math.sin(SI.Conversions.from_deg(angleSortie));
    W = 0.5 * (U2 ^ 2 - U3 ^ 2 - (W2 ^ 2 - W3 ^ 2) + C2 ^ 2 - C3 ^ 2);
    h01 = port_a.h_outflow + 0.5 * C1 ^ 2;
    T02 = h01 / Cp_e;
    h03 = h01 - W;
    T03 = h03 / Cp_s;
    P01 = port_a.p + 0.5 * rho_e * C1 ^ 2;
    P03 = P01 * (T03 / T02) ^ (Cp_e / Cv_e / (Cp_e / Cv_e - 1));
    port_b.p = P03 - 0.5 * rho_s * C3 ^ 2;
    T3 = T03 - 0.5 * C3 ^ 2 / Cp_s;
    Pmeca = mt_flow * W;
    port_b.h_outflow = Medium.specificEnthalpy(etat2);
    //Warning
    // assert(C2 > Medium.velocityOfSound(etat1), "Overload, vitesse fluide > Vitesse du son");
  end turbine;

  partial model turbinePartielle
    extends Modelica.Fluid.Interfaces.PartialTwoPort(port_a(p(start = p_a_start), m_flow(start = m_flow_start, min = 0)), port_b(p(start = p_b_start), m_flow(start = -m_flow_start, max = 0)));
    // Initialization
    parameter Medium.AbsolutePressure p_a_start = system.p_start "Guess value for inlet pressure";
    parameter Medium.AbsolutePressure p_b_start = p_a_start "Guess value for outlet pressure";
    parameter Medium.MassFlowRate m_flow_start = system.m_flow_start "Guess value of m_flow = port_a.m_flow";
    // Variables
    Medium.MassFlowRate mt_flow(start = m_flow_start) "Mass flow rate";
    // Assumptions
    parameter Boolean checkValve = true "to prevent reverse flow";
  equation
    // Conservation masse
    mt_flow = port_a.m_flow;
    port_a.m_flow + port_b.m_flow = 0;
    // Transport de substances
    port_a.Xi_outflow = inStream(port_b.Xi_outflow);
    port_b.Xi_outflow = inStream(port_a.Xi_outflow);
    port_a.C_outflow = inStream(port_b.C_outflow);
    port_b.C_outflow = inStream(port_a.C_outflow);
    // Equilibre Thermodynamique
    port_a.h_outflow = inStream(port_a.h_outflow);
  end turbinePartielle;
end turbineCentrifuge;