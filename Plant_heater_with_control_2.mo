within AES_project_2021_2022;

model Plant_heater_with_control_2
  /* Simulation starts at hs:ms:ss*/
  parameter Real hs = 8;
  parameter Real ms = 0;
  parameter Real ss = 0;
  Real P_loss = Qheat.Q + Hsupz1.Q_flow + Hsupz2.Q_flow;
  AES.ProcessComponents.Thermal.Liquid.Pressuriser pressuriser annotation(
    Placement(visible = true, transformation(origin = {-100, -74}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.Tube tubeh1(L = 50) annotation(
    Placement(visible = true, transformation(origin = {16, -4}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.Tube tubeh2(L = 50) annotation(
    Placement(visible = true, transformation(origin = {230, -4}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.Tube tubec1(L = 50) annotation(
    Placement(visible = true, transformation(origin = {230, -70}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.Tube tubec2(L = 50) annotation(
    Placement(visible = true, transformation(origin = {10, -80}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.Valve_linear vh1(dpnom = 300000, wnom = 0.2) annotation(
    Placement(visible = true, transformation(origin = {142, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  AES.ProcessComponents.Thermal.Liquid.Valve_linear vh2(dpnom = 300000, wnom = 2.2) annotation(
    Placement(visible = true, transformation(origin = {386, 22}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  AES.ProcessComponents.Thermal.Liquid.Tube h1(Di = 0.02, L = 5) annotation(
    Placement(visible = true, transformation(origin = {162, 46}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.Tube h2(L = 50) annotation(
    Placement(visible = true, transformation(origin = {398, 46}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.Pump_centrifugal pump(dp0 = 600000, w0 = 30) annotation(
    Placement(visible = true, transformation(origin = {-34, -4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.Tube closure(L = 1000) annotation(
    Placement(visible = true, transformation(origin = {482, -36}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  AES.ProcessComponents.Thermal.Liquid.Tube tubeheat annotation(
    Placement(visible = true, transformation(origin = {-144, -32}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  AES.ProcessComponents.Thermal.Liquid.DiffPressureSensor sDp annotation(
    Placement(visible = true, transformation(origin = {-54, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.Tsensor sTh annotation(
    Placement(visible = true, transformation(origin = {-176, -4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.surfQcond_prescribed Qheat annotation(
    Placement(visible = true, transformation(origin = {-180, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor Cz1(C = 1e4, T(displayUnit = "K")) annotation(
    Placement(visible = true, transformation(origin = {204, 86}, extent = {{-20, -20}, {20, 20}}, rotation = -90)));
  Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor sTz1 annotation(
    Placement(visible = true, transformation(origin = {90, 86}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  Modelica.Thermal.HeatTransfer.Components.ThermalConductor Gloss1(G = 80) annotation(
    Placement(visible = true, transformation(origin = {162, 116}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow Hsupz1 annotation(
    Placement(visible = true, transformation(origin = {122, 66}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  Modelica.Blocks.Math.Gain Psupz1(k = 500) annotation(
    Placement(visible = true, transformation(origin = {118, 40}, extent = {{-4, -4}, {4, 4}}, rotation = 0)));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature pTa annotation(
    Placement(visible = true, transformation(origin = {-2, 146}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.VectorHPtoHP_conductor convz1 annotation(
    Placement(visible = true, transformation(origin = {162, 66}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow Hsupz2 annotation(
    Placement(visible = true, transformation(origin = {358, 68}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor sTz2 annotation(
    Placement(visible = true, transformation(origin = {280, 126}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor Cz2(C = 1e4, T(displayUnit = "K")) annotation(
    Placement(visible = true, transformation(origin = {446, 88}, extent = {{-20, -20}, {20, 20}}, rotation = -90)));
  AES.ProcessComponents.Thermal.Liquid.VectorHPtoHP_conductor convz2 annotation(
    Placement(visible = true, transformation(origin = {398, 68}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
  Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor(G = 80) annotation(
    Placement(visible = true, transformation(origin = {398, 116}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Modelica.Blocks.Math.Gain Psupz2(k = 500) annotation(
    Placement(visible = true, transformation(origin = {366, 44}, extent = {{-4, -4}, {4, 4}}, rotation = 90)));
  inner AES.ProcessComponents.Thermal.System_settings.System_liquid system(ro(displayUnit = "kg/m3")) annotation(
    Placement(visible = true, transformation(origin = {-180, 148}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.CombiTimeTable sp_Tz(extrapolation = Modelica.Blocks.Types.Extrapolation.Periodic, offset = {273.15}, smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, table = [0, 7; 7, 7; 8, 20; 17, 20; 20, 14; 22, 10; 24, 10], timeScale = 3600) annotation(
    Placement(visible = true, transformation(origin = {-168, 116}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression Pressure_difference(y = 2.027 * 10 ^ 5) annotation(
    Placement(visible = true, transformation(origin = {-219, 60}, extent = {{-25, -18}, {25, 18}}, rotation = 0)));
  AES.ControlBlocks.AnalogueControllers.PI_awfb_basic PI_z1(CSmax = 1, CSmin = 0, K = 0.07692, Ti = 76.92) annotation(
    Placement(visible = true, transformation(origin = {-6, 122}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ControlBlocks.AnalogueControllers.PI_awfb_basic PI_z2(CSmax = 1, CSmin = 0, K = 0.3846, Ti = 76.92) annotation(
    Placement(visible = true, transformation(origin = {276, 72}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ControlBlocks.ActuationSchemes.DaisyChain_uniform daisyChain_z1 annotation(
    Placement(visible = true, transformation(origin = {48, 122}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ControlBlocks.ActuationSchemes.DaisyChain_uniform daisyChain_z2 annotation(
    Placement(visible = true, transformation(origin = {320, 72}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Gain gain(k = 10000) annotation(
    Placement(visible = true, transformation(origin = {-224, -40}, extent = {{-4, -4}, {4, 4}}, rotation = 0)));
  Modelica.Blocks.Sources.CombiTimeTable Tamb(extrapolation = Modelica.Blocks.Types.Extrapolation.Periodic, offset = {273.15}, smoothness = Modelica.Blocks.Types.Smoothness.ContinuousDerivative, table = [0, 0; 4, -2; 8, 8; 12, 10; 15, 10; 18, 3; 20, 1; 22, 0; 24, 0], timeScale = 3600) annotation(
    Placement(visible = true, transformation(origin = {-124, 146}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ControlBlocks.AnalogueControllers.PI_awfb_basic PI_Heater(CSmax = 1, CSmin = 0, K = 0.0047, Ti = 39.2699) annotation(
    Placement(visible = true, transformation(origin = {-312, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression Heater_T_Max(y = 45 + 273.15) annotation(
    Placement(visible = true, transformation(origin = {-443, -34}, extent = {{-25, -18}, {25, 18}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression P_Loss(y = P_loss) annotation(
    Placement(visible = true, transformation(origin = {-401, 138}, extent = {{-19, -10}, {19, 10}}, rotation = 0)));
  Modelica.Blocks.Continuous.Integrator E_loss(initType = Modelica.Blocks.Types.Init.NoInit, use_reset = false) annotation(
    Placement(visible = true, transformation(origin = {-300, 138}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ControlBlocks.AnalogueControllers.PI_awfb_basic PI_Pressure(CSmax = 1, CSmin = 0, K = 0.001, Ti = 1000) annotation(
    Placement(visible = true, transformation(origin = {-42, 32}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ControlBlocks.AnalogueControllers.PI_awfb_basic pI_awfb_basic1(CSmax = 10, CSmin = 0, K = 0.005, Ti = 500) annotation(
    Placement(visible = true, transformation(origin = {268, 14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression LO_limit(y = 5 + 273.15) annotation(
    Placement(visible = true, transformation(origin = {-28, 68}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Logical.Switch switch1 annotation(
    Placement(visible = true, transformation(origin = {76, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Logical.Switch switch2 annotation(
    Placement(visible = true, transformation(origin = {318, 22}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.BooleanExpression booleanExpression(y = if AES_project_2021_2022.NotteGiorno(time) == 0 then false else true)  annotation(
    Placement(visible = true, transformation(origin = {170, 168}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ControlBlocks.AnalogueControllers.PI_awfb_basic pI_awfb_basic(CSmax = 1, CSmin = 0, K = 0.005, Ti = 500) annotation(
    Placement(visible = true, transformation(origin = {18, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
equation
  connect(pump.pwh_b, tubeh1.pwh_a) annotation(
    Line(points = {{-22, -4}, {4, -4}}, color = {46, 52, 54}));
  connect(tubeh1.pwh_b, vh1.pwh_a) annotation(
    Line(points = {{28, -4}, {142, -4}, {142, 6}}, color = {46, 52, 54}));
  connect(vh1.pwh_b, h1.pwh_a) annotation(
    Line(points = {{142, 30}, {142, 46}, {150, 46}}, color = {46, 52, 54}));
  connect(h1.pwh_b, tubec2.pwh_a) annotation(
    Line(points = {{174, 46}, {186, 46}, {186, -80}, {22, -80}}, color = {46, 52, 54}));
  connect(tubeh1.pwh_b, tubeh2.pwh_a) annotation(
    Line(points = {{28, -4}, {218, -4}}, color = {46, 52, 54}));
  connect(tubeh2.pwh_b, vh2.pwh_a) annotation(
    Line(points = {{242, -4}, {386, -4}, {386, 10}}, color = {46, 52, 54}));
  connect(vh2.pwh_b, h2.pwh_a) annotation(
    Line(points = {{386, 34}, {386, 46}}, color = {46, 52, 54}));
  connect(tubec1.pwh_b, tubec2.pwh_a) annotation(
    Line(points = {{218, -70}, {120, -70}, {120, -80}, {22, -80}}, color = {46, 52, 54}));
  connect(h2.pwh_b, tubec1.pwh_a) annotation(
    Line(points = {{410, 46}, {426, 46}, {426, -70}, {242, -70}}, color = {46, 52, 54}));
  connect(pressuriser.pwh_b, tubec2.pwh_b) annotation(
    Line(points = {{-88, -80}, {-2, -80}}, color = {46, 52, 54}));
  connect(tubeh2.pwh_b, closure.pwh_a) annotation(
    Line(points = {{242, -4}, {482, -4}, {482, -24}}, color = {46, 52, 54}));
  connect(closure.pwh_b, tubec1.pwh_a) annotation(
    Line(points = {{482, -48}, {482, -70}, {242, -70}}, color = {46, 52, 54}));
  connect(tubeheat.pwh_b, pump.pwh_a) annotation(
    Line(points = {{-144, -20}, {-144, -4}, {-46, -4}}, color = {46, 52, 54}));
  connect(pressuriser.pwh_a, tubeheat.pwh_a) annotation(
    Line(points = {{-112, -80}, {-144, -80}, {-144, -44}}, color = {46, 52, 54}));
  connect(pump.pwh_b, sDp.pwh_hi) annotation(
    Line(points = {{-22, -4}, {-22, -19}, {-42, -19}, {-42, -32}}, color = {46, 52, 54}));
  connect(sDp.pwh_lo, pressuriser.pwh_b) annotation(
    Line(points = {{-42, -44}, {-42, -59}, {-88, -59}, {-88, -80}}, color = {46, 52, 54}));
  connect(sTh.pwh_a, tubeheat.pwh_b) annotation(
    Line(points = {{-164, -4}, {-144, -4}, {-144, -20}}, color = {46, 52, 54}));
  connect(Qheat.surf, tubeheat.surf) annotation(
    Line(points = {{-180, -48.6667}, {-162, -48.6667}, {-162, -32.3334}, {-150, -32.3334}}, color = {144, 5, 5}));
  connect(Cz1.port, Gloss1.port_b) annotation(
    Line(points = {{184, 86}, {162, 86}, {162, 106}}, color = {191, 0, 0}));
  connect(Psupz1.y, Hsupz1.Q_flow) annotation(
    Line(points = {{122, 40}, {122, 56.4}}, color = {0, 0, 127}));
  connect(Cz1.port, Hsupz1.port) annotation(
    Line(points = {{184, 86}, {122, 86}, {122, 76}}, color = {191, 0, 0}));
  connect(pTa.port, Gloss1.port_a) annotation(
    Line(points = {{8, 146}, {162, 146}, {162, 126}}, color = {191, 0, 0}));
  connect(Cz1.port, convz1.HP) annotation(
    Line(points = {{184, 86}, {162, 86}, {162, 70}}, color = {191, 0, 0}));
  connect(convz1.vectorHP, h1.surf) annotation(
    Line(points = {{162, 62}, {162, 52}}, color = {144, 5, 5}));
  connect(sTz1.port, Cz1.port) annotation(
    Line(points = {{100, 86}, {184, 86}}, color = {191, 0, 0}));
  connect(Psupz2.y, Hsupz2.Q_flow) annotation(
    Line(points = {{366, 48.4}, {366, 52.4}, {358, 52.4}, {358, 58.4}}, color = {0, 0, 127}));
  connect(h2.surf, convz2.vectorHP) annotation(
    Line(points = {{398, 51.4}, {398, 63.4}}, color = {144, 5, 5}));
  connect(convz2.HP, Cz2.port) annotation(
    Line(points = {{398, 72}, {398, 88}, {426, 88}}, color = {191, 0, 0}));
  connect(Cz2.port, Hsupz2.port) annotation(
    Line(points = {{426, 88}, {358, 88}, {358, 78}}, color = {191, 0, 0}));
  connect(sTz2.port, Cz2.port) annotation(
    Line(points = {{290, 126}, {364, 126}, {364, 88}, {426, 88}}, color = {191, 0, 0}));
  connect(pTa.port, thermalConductor.port_a) annotation(
    Line(points = {{8, 146}, {398, 146}, {398, 126}}, color = {191, 0, 0}));
  connect(thermalConductor.port_b, Cz2.port) annotation(
    Line(points = {{398, 106}, {398, 88}, {426, 88}}, color = {191, 0, 0}));
  connect(sTz1.T, PI_z1.PV) annotation(
    Line(points = {{80, 86}, {-18, 86}, {-18, 118}}, color = {0, 0, 127}));
  connect(PI_z1.CS, daisyChain_z1.CSi01) annotation(
    Line(points = {{6, 122}, {36, 122}}, color = {0, 0, 127}));
  connect(sTz2.T, PI_z2.PV) annotation(
    Line(points = {{270, 126}, {270, 102}, {264, 102}, {264, 68}}, color = {0, 0, 127}));
  connect(daisyChain_z2.CSo01[1], Psupz2.u) annotation(
    Line(points = {{332, 72}, {332, 23.5}, {366, 23.5}, {366, 39}}, color = {0, 0, 127}));
  connect(gain.y, Qheat.Q) annotation(
    Line(points = {{-220, -40}, {-187.6, -40}}, color = {0, 0, 127}));
  connect(sp_Tz.y[1], PI_z1.SP) annotation(
    Line(points = {{-157, 116}, {-89.5, 116}, {-89.5, 128}, {-18, 128}}, color = {0, 0, 127}));
  connect(sp_Tz.y[1], PI_z2.SP) annotation(
    Line(points = {{-157, 116}, {-65, 116}, {-65, 180}, {245, 180}, {245, 78}, {264, 78}}, color = {0, 0, 127}));
  connect(Tamb.y[1], pTa.T) annotation(
    Line(points = {{-113, 146}, {-14, 146}}, color = {0, 0, 127}));
  connect(Heater_T_Max.y, PI_Heater.SP) annotation(
    Line(points = {{-415.5, -34}, {-324, -34}}, color = {0, 0, 127}));
  connect(P_Loss.y, E_loss.u) annotation(
    Line(points = {{-380.1, 138}, {-313.1, 138}}, color = {0, 0, 127}));
  connect(PI_Heater.CS, gain.u) annotation(
    Line(points = {{-300, -40}, {-229, -40}}, color = {0, 0, 127}));
  connect(sTh.oT, PI_Heater.PV) annotation(
    Line(points = {{-174, -4}, {-394, -4}, {-394, -44}, {-324, -44}}, color = {0, 0, 127}));
  connect(daisyChain_z1.CSo01[2], Psupz1.u) annotation(
    Line(points = {{60, 122}, {74, 122}, {74, 40}, {113, 40}}, color = {0, 0, 127}));
  connect(PI_Pressure.CS, pump.cmd) annotation(
    Line(points = {{-30, 32}, {-30, 12}, {-34, 12}, {-34, 4}}, color = {0, 0, 127}));
  connect(PI_Pressure.PV, sDp.oDp) annotation(
    Line(points = {{-54, 28}, {-54, -5.2}, {-52, -5.2}, {-52, -38.4}}, color = {0, 0, 127}));
  connect(Pressure_difference.y, PI_Pressure.SP) annotation(
    Line(points = {{-191.5, 60}, {-96.75, 60}, {-96.75, 38}, {-54, 38}}, color = {0, 0, 127}));
  connect(PI_z2.CS, daisyChain_z2.CSi01) annotation(
    Line(points = {{288, 72}, {308, 72}}, color = {0, 0, 127}));
  connect(sTz2.T, pI_awfb_basic1.PV) annotation(
    Line(points = {{270, 126}, {234, 126}, {234, 10}, {256, 10}}, color = {0, 0, 127}));
  connect(switch1.y, vh1.x) annotation(
    Line(points = {{87, 18}, {132, 18}}, color = {0, 0, 127}));
  connect(daisyChain_z1.CSo01[1], switch1.u1) annotation(
    Line(points = {{60, 122}, {52, 122}, {52, 26}, {64, 26}}, color = {0, 0, 127}));
  connect(daisyChain_z2.CSo01[1], switch2.u1) annotation(
    Line(points = {{332, 72}, {290, 72}, {290, 30}, {306, 30}}, color = {0, 0, 127}));
  connect(pI_awfb_basic1.CS, switch2.u3) annotation(
    Line(points = {{280, 14}, {306, 14}}, color = {0, 0, 127}));
  connect(switch2.y, vh2.x) annotation(
    Line(points = {{329, 22}, {376, 22}}, color = {0, 0, 127}));
  connect(booleanExpression.y, switch1.u2) annotation(
    Line(points = {{181, 168}, {64, 168}, {64, 18}}, color = {255, 0, 255}));
  connect(booleanExpression.y, switch2.u2) annotation(
    Line(points = {{181, 168}, {294, 168}, {294, 22}, {306, 22}}, color = {255, 0, 255}));
  connect(pI_awfb_basic.CS, switch1.u3) annotation(
    Line(points = {{30, 18}, {41, 18}, {41, 10}, {64, 10}}, color = {0, 0, 127}));
  connect(sTz1.T, pI_awfb_basic.PV) annotation(
    Line(points = {{80, 86}, {70, 86}, {70, 60}, {-12, 60}, {-12, 14}, {6, 14}}, color = {0, 0, 127}));
  connect(LO_limit.y, pI_awfb_basic.SP) annotation(
    Line(points = {{-17, 68}, {-26, 68}, {-26, 24}, {6, 24}}, color = {0, 0, 127}));
  connect(LO_limit.y, pI_awfb_basic1.SP) annotation(
    Line(points = {{-17, 68}, {246, 68}, {246, 20}, {256, 20}}, color = {0, 0, 127}));
  annotation(
    Diagram(coordinateSystem(extent = {{-480, 180}, {500, -100}})),
    experiment(StartTime = 0, StopTime = 864000, Tolerance = 1e-6, Interval = 86.4),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian",
    __OpenModelica_simulationFlags(lv = "LOG_STATS", s = "dassl"));
end Plant_heater_with_control_2;
