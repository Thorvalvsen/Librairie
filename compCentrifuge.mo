package compresseurCentrifuge
  import Cons = Modelica.Constants;
  import SI = Modelica.SIunits;
  import Math = Modelica.Math;

  model compCourbe
    extends turboMachinePartielle;
    parameter SI.AngularVelocity N;
    parameter SI.Area section1;
    parameter SI.Length diam;
    parameter SI.Angle Beta;
    // etat du fluide
    Medium.ThermodynamicState etat1 = Medium.setState_phX(port_a.p, inStream(port_a.h_outflow), inStream(port_a.Xi_outflow)) "Equation etat entree compresseur";
    Medium.ThermodynamicState etat2 = Medium.setState_pTX(port_b.p, T2, inStream(port_a.Xi_outflow)) "Equation etat sortie compresseur";
    // variable
    SI.Temperature T1;
    SI.Temperature T2;
    SI.Temperature T01, T02;
    SI.Pressure P01, P02;
    //SI.VolumeFlowRate Q;
    Medium.SpecificHeatCapacity Cp = Medium.specificHeatCapacityCp(etat1);
    Medium.SpecificHeatCapacity Cv = Medium.specificHeatCapacityCv(etat1);
    Medium.SpecificHeatCapacity Cp2 = Medium.specificHeatCapacityCp(etat2);
    //Medium.MolarMass wm = Medium.molarMass(etat1);
    Medium.Density rho1 = Medium.density(etat1);
    Medium.Density rho2 = Medium.density(etat2);
    SI.Velocity C1(start = 289) "Vitesse frontal Cx en entree de compresseur";
    SI.Velocity C2 "Vitesse totale en sortie de compresseur";
    SI.Velocity cRad, cTot "Vitesse radiale et totale en sortie de compresseur";
    Real Coeff1, Coeff2 "Coefficient de comportement du compresseur";
    //   Real Coeff1s, Coeff2s "Coefficient simplifie";
    SI.Velocity U = SI.Conversions.from_rpm(N) * diam / 2;
    Real gamma = Cp / Cv;
    SI.Length Head;
    Real N_adim, M_adim(start = 30);
    //   Real N_adims;
  equation
    //Preparation
    T1 = Medium.temperature(etat1);
    //Q = mt_flow / rho1;
    //Energie totale
    // C1 = vitesseEntree(section1, Q);
    C1 = mt_flow / rho1 / section1;
    T01 = tempTotale(T1, C1, Cp);
    P01 = SI.Conversions.to_bar(pressTotale(port_a.p, rho1, C1));
    //Equation adimensionnee compresseur
    N_adim = SI.Conversions.from_rpm(N) / (T01 ^ 0, 5);
    //   N_adims = SI.Conversions.from_rpm(N) / T1 ^ 0.5;
    (Coeff1, Coeff2) = courbeComp(N_adim);
    //   (Coeff1s, Coeff2s) = courbeComp(N_adims);
    //   mt_flow = homotopy(M_adim * P01 / T01 ^ 0.5, (port_b.p - Coeff2s * port_a.p) / (Coeff1s * T1 ^ 0.5));
    M_adim = mt_flow * T01 ^ 0.5 / P01;
    P02 = P01 * (Coeff1 * M_adim + Coeff2);
    //Equation de la compression
    T02 = T01 * (P02 / P01) ^ ((gamma - 1) / gamma) "Manque rendement polytropique";
    Head = Cp * (T02 - T01);
    (cRad, cTot) = vitesseSortie(Head, U, Beta);
    if cRad < 0 then
      C2 = C1;
    else
      C2 = cTot;
    end if;
    //Retour etat du fluide
    port_b.p = SI.Conversions.from_bar(P02) - 0.5 * rho2 * C2 ^ 2;
    T2 = T02 - 0.5 * C2 ^ 2 / Cp;
    port_b.h_outflow = Medium.specificEnthalpy(etat2);
  end compCourbe;

  model compCaract "Compresseur centrifuge sans courbe"
    extends turboMachinePartielle;
    parameter SI.AngularVelocity N(displayUnit = "rpm");
    parameter SI.Area section1;
    parameter SI.Length diam;
    parameter SI.Angle Beta(displayUnit = "°");
    parameter SI.Power Pmeca;
    parameter SI.Efficiency Poly = 0.8;
    // etat du fluide
    Medium.ThermodynamicState etat1 = Medium.setState_phX(port_a.p, inStream(port_a.h_outflow), inStream(port_a.Xi_outflow)) "Equation etat entree compresseur";
    Medium.ThermodynamicState etat2 = Medium.setState_pTX(port_b.p, T2, inStream(port_a.Xi_outflow)) "Equation etat sortie compresseur";
    // variable
    SI.Temperature T1 = Medium.temperature(etat1);
    SI.Temperature T2;
    SI.Temperature T01, T02;
    SI.Pressure P01, P02;
    //SI.VolumeFlowRate Q;
    Medium.SpecificHeatCapacity Cp = Medium.specificHeatCapacityCp(etat1);
    Medium.SpecificHeatCapacity Cv = Medium.specificHeatCapacityCv(etat1);
    Medium.SpecificHeatCapacity Cp2 = Medium.specificHeatCapacityCp(etat2);
    Medium.MolarMass wm = Medium.molarMass(etat1);
    Medium.Density rho1 = Medium.density(etat1);
    Medium.Density rho2 = Medium.density(etat2);
    // Variable du compresseur
    SI.Velocity C1, Co, Cr, C2;
    SI.Velocity U;
  equation
    U = SI.Conversions.from_rpm(N) * diam / 2;
    C1 = mt_flow / rho1 / section1;
    T01 = tempTotale(T1, C1, Cp);
    P01 = pressTotale(port_a.p, rho1, C1);
    T02 = T01 + Pmeca / mt_flow / Cp;
    P02 = P01 * (T02 / T01) ^ (Poly * (Cp / Cv) / (Cp / Cv - 1));
    Co = Pmeca / mt_flow / U;
    // if U > Co then
    Cr = (U - Co) / Math.tan(SI.Conversions.from_deg(Beta));
    // else
    //  Cr = 0;
    // end if;
    C2 = (Co ^ 2 + Cr ^ 2) ^ 0.5;
    T2 = T02 - 0.5 / Cp * C2 ^ 2;
    port_b.p = P02 - 0.5 * rho2 * C2 ^ 2;
    port_b.h_outflow = Medium.specificEnthalpy(etat2);
    //Warning
    assert(U < Co, "Regime instable vitesse trop faible ou puissance trop forte");
    assert(C2 > Medium.velocityOfSound(etat1), "Overload, vitesse fluide > Vitesse du son");
  end compCaract;

  partial model turboMachinePartielle
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
  end turboMachinePartielle;

  function vitesseEntree
    input SI.Area section "Section d'entree du compresseur";
    input SI.VolumeFlowRate Q;
    output SI.Velocity C1;
  algorithm
    C1 := Q / section;
  end vitesseEntree;

  function vitesseSortie
    input SI.Height Head;
    input SI.Velocity U;
    input SI.Angle beta;
    output SI.Velocity cRadiale;
    output SI.Velocity C2;
  protected
    SI.Angle betaRad = SI.Conversions.from_deg(beta);
    SI.Velocity cAng;
  algorithm
    cAng := Head / U;
    cRadiale := (U - cAng) / Math.tan(betaRad);
    C2 := cAng ^ 2 + cRadiale ^ 2;
  end vitesseSortie;

  function courbeComp
    input Real Nadim;
    output Real ak;
    output Real bk;
  protected
    Real aSurge;
    Real bSurge;
    Real aOverload;
    Real bOverload;
    Real coeffSurge1[2] = {0.011, 0.765};
    Real coeffSurge2[2] = {0.3494, -3.8489};
    Real[2] coeffOver1 = {0.0039, 0.9298};
    Real[2] coeffOver2 = {0.4976, 4.1357};
  algorithm
    aSurge := coeffSurge1[1] * Nadim + coeffSurge1[2];
    bSurge := coeffSurge2[1] * Nadim + coeffSurge2[2];
    aOverload := coeffOver1[1] * Nadim + coeffOver1[2];
    bOverload := coeffOver2[1] * Nadim + coeffOver2[2];
    ak := (aOverload - aSurge) / (bOverload - bSurge);
    bk := aSurge - ak * bSurge;
  end courbeComp;

  function tempTotale
    input SI.Temperature T;
    input SI.Velocity C;
    input SI.SpecificHeatCapacity Cp;
    output SI.Temperature T0;
  algorithm
    T0 := T + C ^ 2 / (2 * Cp);
  end tempTotale;

  function pressTotale
    input SI.Pressure P;
    input SI.Density Rho;
    input SI.Velocity C;
    output SI.Pressure P0;
  algorithm
    P0 := P + 0.5 * Rho * C ^ 2;
  end pressTotale;
end compresseurCentrifuge;