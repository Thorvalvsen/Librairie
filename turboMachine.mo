package turboMachine
  import Cons = Modelica.Constants;
  import SI = Modelica.SIunits;
  import Math = Modelica.Math;
  //extends Modelica.Icons.VariantsPackage;
  // Compresseur

  model compCaract "Compresseur centrifuge sans courbe"
    extends turboMachinePartielle;
    extends arbreCompresseur;
    //parameter SI.AngularVelocity N(displayUnit = "rpm");
    parameter SI.Area section1;
    parameter SI.Length diam;
    parameter SI.Angle Beta(displayUnit = "°");
    //parameter SI.Power Pmeca;
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
    SI.Velocity sound = Medium.velocityOfSound(etat1);
    Medium.Density rho1 = Medium.density(etat1);
    Medium.Density rho2 = Medium.density(etat2);
    // Variable du compresseur
    SI.Velocity C1, Co, Cr, C2;
    SI.Velocity U;
  equation
    U = SI.Conversions.from_rpm(N) * diam / 2;
    C1 = mt_flow / rho1 / section1;
    T01 = T1 + 0.5 * C1 ^ 2 / Cp;
    P01 = port_a.p + 0.5 * rho1 * C1 ^ 2;
    T02 = T01 + Pmeca / mt_flow / Cp;
    P02 = P01 * (T02 / T01) ^ (Poly * (Cp / Cv) / (Cp / Cv - 1));
    Co = Pmeca / mt_flow / U;
    if U > Co then
      Cr = (U - Co) / Math.tan(SI.Conversions.from_deg(Beta));
    else
      Cr = 0 "A reflechir";
    end if;
    C2 = (Co ^ 2 + Cr ^ 2) ^ 0.5;
    T2 = T02 - 0.5 / Cp * C2 ^ 2;
    port_b.p = P02 - 0.5 * rho2 * C2 ^ 2;
    port_b.h_outflow = Medium.specificEnthalpy(etat2);
    //Warning
    assert(U < Co, "Regime instable vitesse trop faible ou puissance trop forte");
    //assert(C2 > Medium.velocityOfSound(etat1), "Overload, vitesse fluide > Vitesse du son");
  end compCaract;

  model compresseurAxial "Compresseur axial sans courbe"
    extends turboMachinePartielle;
    extends arbreCompresseur;
    //parameter SI.Conversions.NonSIunits.AngularVelocity_rpm N;
    parameter SI.Diameter diamMobile;
    parameter SI.Diameter diamHub;
    parameter SI.Angle angleGuide;
    parameter SI.Angle angleMobile;
    // parameter SI.Power Pmeca;
    parameter SI.Efficiency nu = 0.9;
    // etat du fluide
    Medium.ThermodynamicState etat1 = Medium.setState_phX(port_a.p, inStream(port_a.h_outflow), inStream(port_a.Xi_outflow)) "Equation etat entree compresseur";
    Medium.ThermodynamicState etat2 = Medium.setState_phX(port_b.p, port_b.h_outflow, inStream(port_a.Xi_outflow)) "Equation etat sortie compresseur";
    // variable
    SI.Area section;
    SI.Velocity Cx, U, C1, C1y, W1, W1y, W2, W2y, C2, C2y;
    SI.Enthalpy h00, h01, h02, h01rel;
    SI.Enthalpy h1, h2;
    SI.Temperature T01, T02;
    SI.Pressure P01, P02;
    SI.Energy dW;
    Real gamma = Cp / Cv;
    Medium.SpecificHeatCapacity Cp = Medium.specificHeatCapacityCp(etat1);
    Medium.SpecificHeatCapacity Cv = Medium.specificHeatCapacityCv(etat1);
    Medium.SpecificHeatCapacity Cp2 = Medium.specificHeatCapacityCp(etat2);
    Medium.Density rho1 = Medium.density(etat1);
    Medium.Density rho2 = Medium.density(etat2);
  equation
    section = Cons.pi * (diamMobile ^ 2 - diamHub ^ 2) / 4;
    Cx = mt_flow / rho1 / section;
    U = SI.Conversions.from_rpm(N) * diamMobile / 2;
    h00 = port_a.h_outflow + 0.5 * Cx ^ 2;
    C1 = Cx / Math.cos(SI.Conversions.from_deg(angleGuide));
    h1 = h00 - 0.5 * C1 ^ 2;
    C1y = C1 * Math.sin(SI.Conversions.from_deg(angleGuide));
    W1y = U - C1y;
    W1 = (Cx ^ 2 + W1y ^ 2) ^ 0.5;
    h01rel = h1 + 0.5 * W1 ^ 2;
    h01 = h1 + 0.5 * C1 ^ 2;
    dW = Pmeca / mt_flow;
    C2y = dW / U + C1y;
    W2y = U - C2y;
    W2 = W2y / Math.sin(SI.Conversions.from_deg(angleMobile));
    h2 = h01rel - 0.5 * W2 ^ 2;
    h02 = h01 + dW;
    C2 = 2 * (h02 - h2) ^ 0.5;
    port_b.h_outflow = h02 - 0.5 * C1 ^ 2;
    T01 = h01 / Cp;
    T02 = h02 / Cp2;
    P01 = port_a.p + 0.5 * rho1 * Cx ^ 2 "Manque rendement Poly le long du guide";
    P02 = P01 * (T02 / T01) ^ (nu * gamma / (gamma - 1));
    port_b.p = P02 - 0.5 * rho2 * C1 ^ 2;
  end compresseurAxial;

  // Turbine

  model turbineAxiale
    extends turboMachinePartielle;
    extends arbreCompresseur;
    parameter SI.Length diamMobile "B11";
    parameter SI.Length diamHub "B12";
    parameter SI.Angle angleGuide(displayUnit = "°") "B13, angle de sortie des guides";
    parameter SI.Angle angleMobile(displayUnit = "°") "B14, angle de sortie du mobile";
    //parameter SI.Conversions.NonSIunits.AngularVelocity_rpm N "B15, vitesse de rotation de la turbine";
    parameter SI.Efficiency nu = 0.9 "Rendement Polytropic";
    // etat du fluide
    Medium.ThermodynamicState etat1 = Medium.setState_phX(port_a.p, inStream(port_a.h_outflow), inStream(port_a.Xi_outflow)) "Equation etat entree turbine";
    Medium.ThermodynamicState etat2 = Medium.setState_phX(port_b.p, port_b.h_outflow, inStream(port_a.Xi_outflow)) "Equation etat sortie turbine";
    // variable (Etage 1 : Guide, Etage 2:Entree rotor, Etage 3:sortie rotor)
    //SI.Temperature T1 = Medium.temperature(etat1);
    //SI.Temperature T2;
    //SI.Temperature T3;
    SI.Temperature T02, T03;
    SI.Pressure P01, P03;
    //SI.Enthalpy h01, h02rel, h03;
    SI.Enthalpy h01, h03;
    SI.Enthalpy h2, h3;
    SI.Velocity U, Cx(start = 100);
    SI.Velocity C2, C2y, W2;
    SI.Velocity C3, C3y, W3;
    SI.Enthalpy W;
    //SI.Power Pmeca;
    SI.Area section;
    //SI.VolumeFlowRate Q;
    Medium.SpecificHeatCapacity Cp_e = Medium.specificHeatCapacityCp(etat1);
    Medium.SpecificHeatCapacity Cv_e = Medium.specificHeatCapacityCv(etat1);
    Medium.SpecificHeatCapacity Cp_s = Medium.specificHeatCapacityCp(etat2);
    Medium.Density rho_e = Medium.density(etat1);
    Medium.Density rho_s = Medium.density(etat2);
    Real gamma = Cp_e / Cv_e;
    // Variable de la turbine
  equation
    section = Cons.pi / 4 * (diamMobile ^ 2 - diamHub ^ 2);
    Cx = mt_flow / rho_e / section;
    C2 = Cx / cos(SI.Conversions.from_deg(angleGuide));
    U = SI.Conversions.from_rpm(N) * diamMobile / 2;
    C2y = C2 * sin(SI.Conversions.from_deg(angleGuide));
    W2 = ((C2y - U) ^ 2 + Cx ^ 2) ^ 0.5;
    h01 = port_a.h_outflow + 0.5 * Cx ^ 2;
    h2 = h01 - 0.5 * C2 ^ 2 "manque rendement polytropique le long des guides";
    //h02rel = h2 + 0.5 * W2 ^ 2;
    W3 = Cx / cos(SI.Conversions.from_deg(angleMobile));
    C3y = W3 * sin(SI.Conversions.from_deg(angleMobile)) - U;
    W = U * (C2y - C3y) "Verifier signe de la relation";
    h03 = h01 - W;
    C3 = (C3y ^ 2 + Cx ^ 2) ^ 0.5;
    h3 = h03 - 0.5 * C3 ^ 2;
    T02 = h01 / Cp_e;
    T03 = h03 / Cp_s;
    //T03 = h03 / Cp_e;
    P01 = port_a.p + 0.5 * rho_e * Cx ^ 2;
    P03 = P01 * (T03 / T02) ^ (gamma / nu / (gamma - 1));
    port_b.p = P03 - 0.5 * rho_s * C3 ^ 2;
    port_b.h_outflow = h3;
    //inStream(port_b.h_outflow) = h3;
    Pmeca = -mt_flow * W;
  end turbineAxiale;

  model turbineCentrifuge
    extends turboMachinePartielle;
    extends arbreCompresseur;
    parameter SI.Length diamEntree "B10";
    parameter SI.Length diamSortie "B12";
    parameter SI.Area epaisseurEntree "B11";
    parameter SI.Angle angleEntree(displayUnit = "°") "B13";
    parameter SI.Angle angleSortie(displayUnit = "°") "B14";
    // etat du fluide
    Medium.ThermodynamicState etat1 = Medium.setState_phX(port_a.p, inStream(port_a.h_outflow), inStream(port_a.Xi_outflow)) "Equation etat entree turbine";
    Medium.ThermodynamicState etat2 = Medium.setState_pTX(port_b.p, T3, inStream(port_a.Xi_outflow)) "Equation etat sortie turbine";
    // variable (Etage 1 : Inducer, Etage 2:Entree rotor, Etage 3:sortie rotor)
    //SI.AngularVelocity N(displayUnit = "rpm");
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
    //SI.Power Pmeca;
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
    Pmeca = -mt_flow * W;
    port_b.h_outflow = Medium.specificEnthalpy(etat2);
    //Warning
    // assert(C2 > Medium.velocityOfSound(etat1), "Overload, vitesse fluide > Vitesse du son");
  end turbineCentrifuge;

  // Pompe centrifuge

  model centrifugeGeometrie "Pompe centrifuge avec definition de la geometrie"
    extends pompeCentrifuge;
    extends arbreCompresseur;
    // Parametres
    parameter SI.Angle Beta;
    parameter SI.Diameter diam;
    parameter Real nombreAillettes;
    parameter SI.Efficiency nu = 0.7;
    parameter SI.AngularVelocity Wiscli = 3 "Coeff de Wisclisius entre 3 et 4 rad/s";
    //parameter Real N "Vitesse de rotation en tr/min";
    parameter SI.Pressure pv "Pression de vapeur";
    //Variable
    SI.VolumeFlowRate Q;
    SI.AngularVelocity omega = SI.Conversions.from_rpm(N);
    //SI.Power P;
    Real phi;
    Real slipping;
    SI.Velocity U2 = bladeSpeed(diam, omega);
    SI.Length NPSHr;
    SI.Length NPSHa;
  equation
    Q = mt_flow / rho;
    phi = flowCoefficient(diam, omega, Q);
    slipping = slipFactor(Beta, nombreAillettes, phi);
    Head = nu * slipping * U2 ^ 2 * (1 - phi * tan(Beta)) / Cons.g_n;
    Pmeca = rho * (N / 60) ^ 3 * diam ^ 5 * nu;
    NPSHr = NPSHWiscli(Q, omega, Wiscli);
    NPSHa = (port_b.p - pv) / (rho * Cons.g_n);
  end centrifugeGeometrie;

  model centrifugeCourbe "Pompe centrifuge avec courbe"
    extends pompeCentrifuge;
    extends arbreCompresseur;
    // Parametres
    parameter SI.Diameter diam;
    parameter SI.AngularVelocity Wiscli;
    //parameter Real N "Vitesse de rotation en tr/min";
    parameter SI.Pressure pv "Pression de vapeur fluide";
    // Variables
    SI.VolumeFlowRate Q;
    SI.AngularVelocity omega = SI.Conversions.from_rpm(N);
    Real phi;
    Real psi;
    SI.Length NPSHr;
    SI.Length NPSHa;
    Real nu;
    //SI.Power P;
  equation
    Q = mt_flow / rho;
    phi = flowCoefficient(diam, omega, Q);
    nu = fnu(phi);
    NPSHr = NPSHWiscli(Q, omega, Wiscli);
    psi = fpsi(phi);
    Head = psi * (N / 60 * diam) ^ 2 / Cons.g_n;
    Pmeca = rho * (N / 60) ^ 3 * diam ^ 5 * nu;
    NPSHa = (port_b.p - pv) / (rho * Cons.g_n);
  end centrifugeCourbe;

  //Turbine hydraulique

  model Pelton01 "Turbine de pelton avec determination de la vitesse en fonction de la meilleure efficacite"
    extends liquidIncomp;
    extends arbreCompresseur;
    // Definition des parametres
    parameter SI.Angle Beta "angle de sortie des aubes >90";
    parameter SI.Diameter diamMobile "Diametre du mobile";
    parameter SI.Diameter diamNozzle "diametre des injecteurs";
    parameter Real nombreNozzle "Nombre d'injecteur";
    parameter SI.Efficiency nu = 0.9 "Rendement total de la turbine";
    //Variables internes
    SI.Velocity C1 "Vitesse du fluide en sortie d'injecteur";
    SI.Velocity U "Vitesse en bout d'aube";
    //SI.AngularVelocity N "Vitesse de rotation du mobile";
    //SI.Power Pmeca "Puissance mecanique";
  equation
    C1 = 0.98 * (2 * Cons.g_n * Head) ^ 0.5;
    U = C1 / 2 "Point de meilleur rendement";
    mt_flow = nombreNozzle * C1 * Cons.pi * diamNozzle ^ 2 / 4 * rho;
    Pmeca = -mt_flow * U * (C1 - U) * (1 - nu * cos(SI.Conversions.from_deg(Beta)));
    N = U / (diamMobile / 2) / Cons.pi * 30;
  end Pelton01;

  model Pelton02 "Turbine de pelton avec calcul du rendement en fonction de la vitesse de rotation"
    extends liquidIncomp;
    extends arbreCompresseur;
    // Definition des parametres
    parameter SI.Angle Beta "angle de sortie des aubes >90";
    parameter SI.Diameter diamMobile "Diametre du mobile";
    parameter SI.Diameter diamNozzle "diametre des injecteurs";
    parameter Real nombreNozzle "Nombre d'injecteur";
    //parameter SI.AngularVelocity N "Vitesse de rotation";
    //Variables internes
    SI.Velocity C1 "Vitesse du fluide en sortie d'injecteur";
    SI.Velocity U "Vitesse en bout d'aube";
    SI.Efficiency nu "Rendement de la turbine";
    //SI.Power Pmeca "Puissance mecanique";
  equation
    C1 = 0.98 * (2 * Cons.g_n * Head) ^ 0.5;
    U = SI.Conversions.from_rpm(N) * diamMobile / 2;
    nu = rendPelton(U, C1);
    mt_flow = nombreNozzle * C1 * Cons.pi * diamNozzle ^ 2 / 4 * rho;
    Pmeca = -nu * mt_flow * U * (C1 - U) * (1 - 0.9 * cos(SI.Conversions.from_deg(Beta)));
  end Pelton02;

  model Francis
    extends liquidIncomp;
    extends arbreCompresseur;
    // Definition des parametres
    parameter SI.Angle Beta "angle des guides";
    parameter SI.Length guide "hauteur des guides";
    parameter SI.Diameter diamMobile "Diametre de la turbine";
    parameter SI.Efficiency nu "Rendement hydraulique";
    //parameter SI.AngularVelocity N "Vitesse de rotation";
    parameter SI.Pressure psat = 1000 "Pression de vapeur saturante, evolution -> Monitoring fluide biphasique";
    // Definition des variables
    SI.Velocity U "Vitesses des aubes";
    SI.Velocity C2 "Vitesse du fluide";
    SI.Energy W "Travail unitaire";
    //SI.Power Pmeca "Puissance mécanique";
    SI.AngularVelocity PSS "Power Specific Speed";
    Real coeffCavit "Coefficient de cavitation";
    Real coeffThoma "Coefficient de Thoma";
  equation
    U = SI.Conversions.from_rpm(N) * diamMobile / 2;
    W = nu * Head * Cons.g_n;
    C2 = W / U / sin(SI.Conversions.from_deg(Beta));
    mt_flow = Cons.pi * diamMobile * guide * C2 * cos(SI.Conversions.from_deg(Beta)) * rho;
    Pmeca = -mt_flow * W;
    PSS = SI.Conversions.from_rpm(N) * (Pmeca / rho) ^ 0.5 / (Cons.g_n * Head) ^ (5 / 4);
    coeffCavit = coeffCavitation(PSS);
    coeffThoma = (port_b.p - psat) / (rho * Cons.g_n) / Head;
    //   assert(coeffThoma < coeffCavit, "Cavitation");
  end Francis;

  // Base Classe //

  partial model turboMachinePartielle
    //extends Modelica.Fluid.Interfaces.PartialTwoPort(port_a(p(start = p_a_start), m_flow(start = m_flow_start, min = 0)), port_b(p(start = p_b_start), m_flow(start = -m_flow_start, max = 0)));
    extends Modelica.Fluid.Interfaces.PartialTwoPort(port_a(p(start = p_a_start), m_flow(start = m_flow_start, min = 0)), port_b(p(start = p_b_start), m_flow(start = -m_flow_start)));
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

  partial model arbreCompresseur
    SI.Power Pmeca "Puisance mecanique de arbre";
    SI.Conversions.NonSIunits.AngularVelocity_rpm N;
    SI.Angle phi;
    Modelica.Mechanics.Rotational.Interfaces.Flange_a shaft;
  equation
    phi = shaft.phi;
    N = SI.Conversions.to_rpm(der(phi));
    Pmeca = SI.Conversions.from_rpm(N) * shaft.tau;
  end arbreCompresseur;

  partial model pompeCentrifuge
    extends Modelica.Fluid.Interfaces.PartialTwoPort(port_a(m_flow(start = m_flow_start)));
    // Initialization
    //parameter Medium.AbsolutePressure p_a_start=system.p_start "Guess value for inlet pressure";
    //parameter Medium.AbsolutePressure p_b_start=p_a_start "Guess value for outlet pressure";
    parameter Medium.MassFlowRate m_flow_start = system.m_flow_start "Guess value of m_flow = port_a.m_flow";
    // Variables
    Medium.MassFlowRate mt_flow "Mass flow rate";
    Medium.ThermodynamicState etat = Medium.setState_phX(port_a.p, inStream(port_a.h_outflow), inStream(port_a.Xi_outflow)) "Equation etat";
    Medium.Density rho = Medium.density(etat) "Masse volumique";
    SI.Length Head "Hauteur manometrique";
    // Assumptions
    // parameter Boolean checkValve=true "to prevent reverse flow";
  equation
    // Conservation masse
    mt_flow = port_a.m_flow;
    port_a.m_flow + port_b.m_flow = 0;
    // Prise de hauteur
    Head = (port_b.p - port_a.p) / (rho * Cons.g_n);
    // Enthalpie
    port_a.h_outflow = inStream(port_b.h_outflow);
    port_b.h_outflow = inStream(port_a.h_outflow);
    // Transport de substances
    port_a.Xi_outflow = inStream(port_b.Xi_outflow);
    port_b.Xi_outflow = inStream(port_a.Xi_outflow);
    port_a.C_outflow = inStream(port_b.C_outflow);
    port_b.C_outflow = inStream(port_a.C_outflow);
  end pompeCentrifuge;

  partial model liquidIncomp
    extends Modelica.Fluid.Interfaces.PartialTwoPort(port_a(m_flow(start = m_flow_start)));
    // Initialization
    //parameter Medium.AbsolutePressure p_a_start=system.p_start "Guess value for inlet pressure";
    //parameter Medium.AbsolutePressure p_b_start=p_a_start "Guess value for outlet pressure";
    parameter Medium.MassFlowRate m_flow_start = system.m_flow_start "Guess value of m_flow = port_a.m_flow";
    // Variables
    Medium.MassFlowRate mt_flow "Mass flow rate";
    Medium.ThermodynamicState etat = Medium.setState_phX(port_a.p, inStream(port_a.h_outflow), inStream(port_a.Xi_outflow)) "Equation etat";
    Medium.Density rho = Medium.density(etat) "Masse volumique";
    SI.Length Head "Hauteur manometrique";
    // Assumptions
    // parameter Boolean checkValve=true "to prevent reverse flow";
  equation
    // Conservation masse
    mt_flow = port_a.m_flow;
    port_a.m_flow + port_b.m_flow = 0;
    // Prise de hauteur
    Head = (port_a.p - port_b.p) / (rho * Cons.g_n);
    // Enthalpie
    port_a.h_outflow = inStream(port_b.h_outflow);
    port_b.h_outflow = inStream(port_a.h_outflow);
    // Transport de substances
    port_a.Xi_outflow = inStream(port_b.Xi_outflow);
    port_b.Xi_outflow = inStream(port_a.Xi_outflow);
    port_a.C_outflow = inStream(port_b.C_outflow);
    port_b.C_outflow = inStream(port_a.C_outflow);
  end liquidIncomp;

  function slipFactor "facteur de glissement"
    input SI.Angle Beta "Angle de sortie du fluide";
    input Real nombreAilettes;
    input Real flowCoef "rapport entre la vitesse axiale et la vitesse de rotation";
    output Real slip;
  algorithm
    if SI.Conversions.to_deg(Beta) < 50 then
      slip := 1 - 0.63 * Cons.pi / nombreAilettes / (1 - flowCoef * tan(Beta));
    else
      slip := 1 - Cons.pi / nombreAilettes * cos(Beta) / (1 - flowCoef * tan(Beta));
    end if;
  end slipFactor;

  function flowCoefficient
    input SI.Diameter diam "diametre de la roue";
    input SI.AngularVelocity omega "Vitesse de rotation en rad/s";
    input SI.VolumeFlowRate Q "Debit Volumique";
    output Real phi;
  algorithm
    phi := Q / (omega / (2 * Cons.pi) * diam ^ 3);
  end flowCoefficient;

  function bladeSpeed "Vitesse de rotation de la roue"
    input SI.Diameter diam "diametre de la roue";
    input SI.AngularVelocity omega "Vitesse angulaire de la roue";
    output SI.Velocity U;
  algorithm
    U := omega * 2 * Cons.pi * diam / 2;
  end bladeSpeed;

  function NPSHWiscli "Calcul du NPSH selon Wisclisius"
    input SI.VolumeFlowRate Q "Debit Volumique";
    input SI.AngularVelocity omega "Vitesse de rotation angulaire";
    input SI.AngularVelocity Wiscli = 3;
    output SI.Length NPSH;
  algorithm
    NPSH := (omega / Wiscli) ^ (4 / 3) * Q ^ (2 / 3) / Cons.g_n;
  end NPSHWiscli;

  function fpsi "premiere caracteristique de pompe"
    input Real phi "Rapport de Q/ND3";
    output Real psi "gH/(ND)² fonction de phi";
  protected
    constant Real a = 280;
    constant Real b = 6.68;
    constant Real c = 6.0755;
  algorithm
    psi := a * phi ^ 2 + b * phi + c;
  end fpsi;

  function fnu "seconde caracteristique de pompe"
    input Real phi;
    output Real nu "rapport de P/(N³Q⁵)";
  protected
    constant Real a = 3.56;
    constant Real b = 0.1505;
  algorithm
    nu := a * phi + b;
  end fnu;

  function rendPelton "Rendement de la turbine de pelton en fonction du rapport des vitesses"
    input SI.Velocity U;
    input SI.Velocity C1;
    output SI.Efficiency nu;
  algorithm
    nu := 3.767 * (U / C1) * (1 - U / C1);
  end rendPelton;

  function coeffCavitation "Calcul de la valeur limite du coeff de thoma en fonction de la PSS"
    input SI.AngularVelocity PSS;
    output Real coeffCavit;
  algorithm
    coeffCavit := 0.03 * exp(0.8 * PSS);
  end coeffCavitation;
end turboMachine;