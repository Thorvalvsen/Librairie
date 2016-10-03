package PompeCentrifuge

import Cons=Modelica.Constants;
import SI=Modelica.SIunits;

function slipFactor "facteur de glissement"
  input SI.Angle Beta "Angle de sortie du fluide";
  input Real nombreAilettes;
  input Real flowCoef "rapport entre la vitesse axiale et la vitesse de rotation";
  output Real slip;
algorithm
  if SI.Conversions.to_deg(Beta)<50 then
    slip:=1-0.63*Cons.pi/nombreAilettes/(1-flowCoef*tan(Beta));
    else 
      slip:=1-Cons.pi/nombreAilettes*cos(Beta)/(1-flowCoef*tan(Beta));
end if;  
end slipFactor;
  
function flowCoefficient
  input SI.Diameter diam "diametre de la roue";
  input SI.AngularVelocity omega "Vitesse de rotation en rad/s";
  input SI.VolumeFlowRate Q "Debit Volumique" ;
  output Real phi;
algorithm
  phi:=Q/(omega/(2*Cons.pi)*diam^3);
end flowCoefficient;

function bladeSpeed "Vitesse de rotation de la roue"
  input SI.Diameter diam "diametre de la roue";
  input SI.AngularVelocity omega "Vitesse angulaire de la roue";
  output SI.Velocity U;
algorithm
  U:=omega*2*Cons.pi*diam/2;
end bladeSpeed;

function NPSHWiscli "Calcul du NPSH selon Wisclisius"
  input SI.VolumeFlowRate Q "Debit Volumique";
  input SI.AngularVelocity omega "Vitesse de rotation angulaire";
  input SI.AngularVelocity Wiscli=3;
  output SI.Length NPSH;
algorithm
  NPSH:=(omega/Wiscli)^(4/3)*Q^(2/3)/Cons.g_n;
end NPSHWiscli;
  
function fpsi "premiere caracteristique de pompe"
  input Real phi "Rapport de Q/ND3";
  output Real psi "gH/(ND)² fonction de phi";    
  
  protected
  constant Real a=280;
  constant Real b=6.68;
  constant Real c=6.0755;
  
algorithm
  psi:=a*phi^2+b*phi+c;
  
end fpsi;
  
function fnu "seconde caracteristique de pompe"
  input Real phi;
  output Real nu "rapport de P/(N³Q⁵)";
  
  protected
  constant Real a=3.56;
  constant Real b=0.1505;
  
algorithm
  nu:=a*phi+b;
  
end fnu;
 
 
partial model pompeCentrifuge
  
  extends Modelica.Fluid.Interfaces.PartialTwoPort(port_a(m_flow(start = m_flow_start)));
   
  // Initialization
  //parameter Medium.AbsolutePressure p_a_start=system.p_start "Guess value for inlet pressure";
  //parameter Medium.AbsolutePressure p_b_start=p_a_start "Guess value for outlet pressure";
  parameter Medium.MassFlowRate m_flow_start = system.m_flow_start "Guess value of m_flow = port_a.m_flow";
 
  // Variables

  Medium.MassFlowRate mt_flow "Mass flow rate";
  Medium.ThermodynamicState etat=Medium.setState_phX(port_a.p,inStream(port_a.h_outflow),inStream(port_a.Xi_outflow)) "Equation etat";
  Medium.Density rho=Medium.density(etat) "Masse volumique";
  SI.Length Head "Hauteur manometrique";
   // Assumptions
   // parameter Boolean checkValve=true "to prevent reverse flow";
  
equation
  
  
  // Conservation masse
  mt_flow=port_a.m_flow;
  port_a.m_flow+port_b.m_flow=0;
  // Prise de hauteur
  Head=(port_b.p-port_a.p)/(rho*Cons.g_n);
  // Enthalpie
  port_a.h_outflow = inStream(port_b.h_outflow);
  port_b.h_outflow = inStream(port_a.h_outflow);
  // Transport de substances
  port_a.Xi_outflow = inStream(port_b.Xi_outflow);
  port_b.Xi_outflow = inStream(port_a.Xi_outflow);
  port_a.C_outflow = inStream(port_b.C_outflow);
  port_b.C_outflow = inStream(port_a.C_outflow);

end pompeCentrifuge;


model centrifugeGeometrie
   
  extends pompeCentrifuge;
  // Parametres
  parameter SI.Angle Beta;
  parameter SI.Diameter diam;
  parameter Real nombreAillettes;
  parameter SI.Efficiency nu=0.7;
  parameter SI.AngularVelocity Wiscli=3 "Coeff de Wisclisius entre 3 et 4 rad/s";
  parameter Real N "Vitesse de rotation en tr/min";
  parameter SI.Pressure pv "Pression de vapeur";
  //Variable
  SI.VolumeFlowRate Q;
  SI.AngularVelocity omega=N*2*Cons.pi/60;
  SI.Power P;
  Real phi;
  Real slipping;
  SI.Velocity U2=bladeSpeed(diam,omega);
  SI.Length NPSHr;
  SI.Length NPSHa;
  
equation
  Q=mt_flow/rho;
  phi=flowCoefficient(diam,omega,Q);
  slipping=slipFactor(Beta,nombreAillettes,phi);
  Head=nu*slipping*U2^2*(1-phi*tan(Beta))/Cons.g_n;
  P=rho*(N/60)^3*(diam)^5*nu;
  NPSHr=NPSHWiscli(Q,omega,Wiscli);
  NPSHa=(port_b.p-pv)/(rho*Cons.g_n);
    
end centrifugeGeometrie;
  
model centrifugeCourbe
  
  extends pompeCentrifuge;
  // Parametres
  parameter SI.Diameter diam;
  parameter SI.AngularVelocity Wiscli;
  parameter Real N "Vitesse de rotation en tr/min";
  parameter SI.Pressure pv "Pression de vapeur fluide";
  // Variables
  SI.VolumeFlowRate Q;
  SI.AngularVelocity omega=N*2*Cons.pi/60;
  Real phi;
  Real psi;
  SI.Length NPSHr;
  SI.Length NPSHa;  
  Real nu;
  SI.Power P;   
    
equation
  Q=mt_flow/rho;
  phi=flowCoefficient(diam,omega,Q);
  nu=fnu(phi);
  NPSHr=NPSHWiscli(Q,omega,Wiscli);
  psi=fpsi(phi);
  Head=psi*(N/60*diam)^2/Cons.g_n;
  P=rho*(N/60)^3*(diam)^5*nu;
  NPSHa=(port_b.p-pv)/(rho*Cons.g_n);
    
end centrifugeCourbe;
  
end PompeCentrifuge;
    
  
  
  
 