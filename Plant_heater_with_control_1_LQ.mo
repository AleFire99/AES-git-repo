within AES_project_2021_2022;

model Plant_heater_with_control_1_LQ
  /* Simulation starts at hs:ms:ss*/
  parameter Real hs = 8;
  parameter Real ms = 0;
  parameter Real ss = 0;
  Real P_loss = Qheat.Q + Hsupz1.Q_flow + Hsupz2.Q_flow;
 
  AES.ProcessComponents.Thermal.Liquid.Pressuriser pressuriser annotation(
    Placement(visible = true, transformation(origin = {-108, -64}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.Tube tubeh1(L = 50) annotation(
    Placement(visible = true, transformation(origin = {10, -4}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.Tube tubeh2(L = 50) annotation(
    Placement(visible = true, transformation(origin = {230, -4}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.Tube tubec1(L = 50) annotation(
    Placement(visible = true, transformation(origin = {230, -70}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.Tube tubec2(L = 50) annotation(
    Placement(visible = true, transformation(origin = {10, -70}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.Valve_linear vh1(dpnom = 300000, wnom = 0.2) annotation(
    Placement(visible = true, transformation(origin = {142, 26}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  AES.ProcessComponents.Thermal.Liquid.Valve_linear vh2(dpnom = 300000, wnom = 2.2) annotation(
    Placement(visible = true, transformation(origin = {386, 26}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
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
    Placement(visible = true, transformation(origin = {192, 86}, extent = {{-20, -20}, {20, 20}}, rotation = -90)));
  Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor sTz1 annotation(
    Placement(visible = true, transformation(origin = {90, 86}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  Modelica.Thermal.HeatTransfer.Components.ThermalConductor Gloss1(G = 80) annotation(
    Placement(visible = true, transformation(origin = {162, 116}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow Hsupz1 annotation(
    Placement(visible = true, transformation(origin = {122, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  Modelica.Blocks.Math.Gain Psupz1(k = 500) annotation(
    Placement(visible = true, transformation(origin = {122, 52}, extent = {{-4, -4}, {4, 4}}, rotation = 90)));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature pTa annotation(
    Placement(visible = true, transformation(origin = {-2, 146}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.VectorHPtoHP_conductor convz1 annotation(
    Placement(visible = true, transformation(origin = {162, 66}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow Hsupz2 annotation(
    Placement(visible = true, transformation(origin = {358, 72}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor sTz2 annotation(
    Placement(visible = true, transformation(origin = {292, 88}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor Cz2(C = 1e4, T(displayUnit = "K")) annotation(
    Placement(visible = true, transformation(origin = {446, 88}, extent = {{-20, -20}, {20, 20}}, rotation = -90)));
  AES.ProcessComponents.Thermal.Liquid.VectorHPtoHP_conductor convz2 annotation(
    Placement(visible = true, transformation(origin = {398, 68}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
  Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor(G = 80) annotation(
    Placement(visible = true, transformation(origin = {398, 116}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Modelica.Blocks.Math.Gain Psupz2(k = 500) annotation(
    Placement(visible = true, transformation(origin = {358, 54}, extent = {{-4, -4}, {4, 4}}, rotation = 90)));
  inner AES.ProcessComponents.Thermal.System_settings.System_liquid system annotation(
    Placement(visible = true, transformation(origin = {-182, 148}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ControlBlocks.AnalogueControllers.PI_awfb_basic PI_Pressure(CSmax = 1, CSmin = 0, K = 100)  annotation(
    Placement(visible = true, transformation(origin = {-36, 54}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression Pressure_difference(y = 2 * 10 ^ 5)  annotation(
    Placement(visible = true, transformation(origin = {-219, 60}, extent = {{-25, -18}, {25, 18}}, rotation = 0)));
  Modelica.Blocks.Math.Gain gain(k = 10000) annotation(
    Placement(visible = true, transformation(origin = {-212, -40}, extent = {{-4, -4}, {4, 4}}, rotation = 0)));
  Modelica.Blocks.Sources.CombiTimeTable Tamb(extrapolation = Modelica.Blocks.Types.Extrapolation.Periodic, offset = {273.15}, smoothness = Modelica.Blocks.Types.Smoothness.ContinuousDerivative, table = [0, 0; 4, -2; 8, 8; 12, 10; 15, 10; 18, 3; 20, 1; 22, 0; 24, 0], timeScale = 3600) annotation(
    Placement(visible = true, transformation(origin = {-124, 146}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
 Modelica.Blocks.Sources.RealExpression P_Loss(y = P_loss)  annotation(
    Placement(visible = true, transformation(origin = {-401, 138}, extent = {{-19, -10}, {19, 10}}, rotation = 0)));
 Modelica.Blocks.Continuous.Integrator E_loss(initType = Modelica.Blocks.Types.Init.NoInit, use_reset = false)  annotation(
    Placement(visible = true, transformation(origin = {-300, 138}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
 Modelica.Blocks.Sources.CombiTimeTable sp_Tz1(extrapolation = Modelica.Blocks.Types.Extrapolation.Periodic, offset = {273.15}, smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, table = [0, 7; 7, 7; 8, 20; 17, 20; 20, 14; 22, 10; 24, 10], timeScale = 3600) annotation(
    Placement(visible = true, transformation(origin = {-28, 104}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
 Modelica.Blocks.Sources.CombiTimeTable sp_Tz2(extrapolation = Modelica.Blocks.Types.Extrapolation.Periodic, offset = {273.15}, smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, table = [0, 7; 7, 7; 8, 20; 17, 20; 20, 14; 22, 10; 24, 10], timeScale = 3600) annotation(
    Placement(visible = true, transformation(origin = {218, 118}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
 Modelica.Blocks.Sources.RealExpression realExpression(y = 45 + 273.15) annotation(
    Placement(visible = true, transformation(origin = {-445, -34}, extent = {{-25, -18}, {25, 18}}, rotation = 0)));
 AES.ControlBlocks.AnalogueControllers.PI_awfb_basic PI_Heater(CSmax = 1, CSmin = 0, K = 0.0047, Ti = 39.2699) annotation(
    Placement(visible = true, transformation(origin = {-304, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
 Modelica.Blocks.Math.Add ref_z1(k2 = -1)  annotation(
    Placement(visible = true, transformation(origin = {2, 44}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
 Modelica.Blocks.Math.Add ref_z2(k2 = -1) annotation(
    Placement(visible = true, transformation(origin = {238, 46}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
 Modelica.Blocks.Continuous.Integrator integrator_z2 annotation(
    Placement(visible = true, transformation(origin = {258, 26}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
 Modelica.Blocks.Continuous.Integrator integrator_z1 annotation(
    Placement(visible = true, transformation(origin = {16, 24}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
 Modelica.Blocks.Math.MatrixGain matrixGain_z2(K = [-102.2973, 99.9993; -0.0120, 0.0117]) annotation(
    Placement(visible = true, transformation(origin = {320, 44}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
 Modelica.Blocks.Routing.Multiplex2 mux_z2 annotation(
    Placement(visible = true, transformation(origin = {288, 44}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
 Modelica.Blocks.Math.MatrixGain matrixGain_z1(K = [-102.2973, 99.9993; -0.0120, 0.0117]) annotation(
    Placement(visible = true, transformation(origin = {84, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
 Modelica.Blocks.Routing.Multiplex2 mux_z1 annotation(
    Placement(visible = true, transformation(origin = {54, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
 Modelica.Blocks.Nonlinear.Limiter sat2_z1(limitsAtInit = true, uMax = 1, uMin = 0)  annotation(
    Placement(visible = true, transformation(origin = {108, 42}, extent = {{-4, -4}, {4, 4}}, rotation = 0)));
 Modelica.Blocks.Nonlinear.Limiter sat1_z1(limitsAtInit = true, uMax = 1, uMin = 0) annotation(
    Placement(visible = true, transformation(origin = {108, 26}, extent = {{-4, -4}, {4, 4}}, rotation = 0)));
 Modelica.Blocks.Nonlinear.Limiter sat2_z2(limitsAtInit = true, uMax = 1, uMin = 0) annotation(
    Placement(visible = true, transformation(origin = {344, 44}, extent = {{-4, -4}, {4, 4}}, rotation = 0)));
 Modelica.Blocks.Nonlinear.Limiter sat1_z2(limitsAtInit = true, uMax = 1, uMin = 0) annotation(
    Placement(visible = true, transformation(origin = {346, 26}, extent = {{-4, -4}, {4, 4}}, rotation = 0)));
equation
  connect(pump.pwh_b, tubeh1.pwh_a) annotation(
    Line(points = {{-22, -4}, {-2, -4}}, color = {46, 52, 54}));
  connect(tubeh1.pwh_b, vh1.pwh_a) annotation(
    Line(points = {{22, -4}, {142, -4}, {142, 14}}, color = {46, 52, 54}));
  connect(vh1.pwh_b, h1.pwh_a) annotation(
    Line(points = {{142, 38}, {142, 46}, {150, 46}}, color = {46, 52, 54}));
  connect(h1.pwh_b, tubec2.pwh_a) annotation(
    Line(points = {{174, 46}, {186, 46}, {186, -70}, {22, -70}}, color = {46, 52, 54}));
  connect(tubeh1.pwh_b, tubeh2.pwh_a) annotation(
    Line(points = {{22, -4}, {218, -4}}, color = {46, 52, 54}));
  connect(tubeh2.pwh_b, vh2.pwh_a) annotation(
    Line(points = {{242, -4}, {386, -4}, {386, 14}}, color = {46, 52, 54}));
  connect(vh2.pwh_b, h2.pwh_a) annotation(
    Line(points = {{386, 38}, {386, 46}}, color = {46, 52, 54}));
  connect(tubec1.pwh_b, tubec2.pwh_a) annotation(
    Line(points = {{218, -70}, {22, -70}}, color = {46, 52, 54}));
  connect(h2.pwh_b, tubec1.pwh_a) annotation(
    Line(points = {{410, 46}, {426, 46}, {426, -70}, {242, -70}}, color = {46, 52, 54}));
  connect(pressuriser.pwh_b, tubec2.pwh_b) annotation(
    Line(points = {{-96, -70}, {-2, -70}}, color = {46, 52, 54}));
  connect(tubeh2.pwh_b, closure.pwh_a) annotation(
    Line(points = {{242, -4}, {482, -4}, {482, -24}}, color = {46, 52, 54}));
  connect(closure.pwh_b, tubec1.pwh_a) annotation(
    Line(points = {{482, -48}, {482, -70}, {242, -70}}, color = {46, 52, 54}));
  connect(tubeheat.pwh_b, pump.pwh_a) annotation(
    Line(points = {{-144, -20}, {-144, -4}, {-46, -4}}, color = {46, 52, 54}));
  connect(pressuriser.pwh_a, tubeheat.pwh_a) annotation(
    Line(points = {{-120, -70}, {-144, -70}, {-144, -44}}, color = {46, 52, 54}));
  connect(pump.pwh_b, sDp.pwh_hi) annotation(
    Line(points = {{-22, -4}, {-22, -19}, {-42, -19}, {-42, -32}}, color = {46, 52, 54}));
  connect(sDp.pwh_lo, pressuriser.pwh_b) annotation(
    Line(points = {{-42, -44}, {-42, -70}, {-96, -70}}, color = {46, 52, 54}));
  connect(sTh.pwh_a, tubeheat.pwh_b) annotation(
    Line(points = {{-164, -4}, {-144, -4}, {-144, -20}}, color = {46, 52, 54}));
  connect(Qheat.surf, tubeheat.surf) annotation(
    Line(points = {{-180, -48.6667}, {-162, -48.6667}, {-162, -32.3334}, {-150, -32.3334}}, color = {144, 5, 5}));
  connect(Cz1.port, Gloss1.port_b) annotation(
    Line(points = {{172, 86}, {162, 86}, {162, 106}}, color = {191, 0, 0}));
 connect(Cz1.port, Hsupz1.port) annotation(
    Line(points = {{172, 86}, {122, 86}, {122, 80}}, color = {191, 0, 0}));
  connect(pTa.port, Gloss1.port_a) annotation(
    Line(points = {{8, 146}, {162, 146}, {162, 126}}, color = {191, 0, 0}));
  connect(Cz1.port, convz1.HP) annotation(
    Line(points = {{172, 86}, {162, 86}, {162, 70}}, color = {191, 0, 0}));
  connect(convz1.vectorHP, h1.surf) annotation(
    Line(points = {{162, 62}, {162, 52}}, color = {144, 5, 5}));
  connect(sTz1.port, Cz1.port) annotation(
    Line(points = {{100, 86}, {172, 86}}, color = {191, 0, 0}));
  connect(h2.surf, convz2.vectorHP) annotation(
    Line(points = {{398, 51.4}, {398, 63.4}}, color = {144, 5, 5}));
  connect(convz2.HP, Cz2.port) annotation(
    Line(points = {{398, 72}, {398, 88}, {426, 88}}, color = {191, 0, 0}));
 connect(Cz2.port, Hsupz2.port) annotation(
    Line(points = {{426, 88}, {358, 88}, {358, 82}}, color = {191, 0, 0}));
  connect(sTz2.port, Cz2.port) annotation(
    Line(points = {{302, 88}, {426, 88}}, color = {191, 0, 0}));
  connect(pTa.port, thermalConductor.port_a) annotation(
    Line(points = {{8, 146}, {398, 146}, {398, 126}}, color = {191, 0, 0}));
  connect(thermalConductor.port_b, Cz2.port) annotation(
    Line(points = {{398, 106}, {398, 88}, {426, 88}}, color = {191, 0, 0}));
  connect(Pressure_difference.y, PI_Pressure.SP) annotation(
    Line(points = {{-191.5, 60}, {-48, 60}}, color = {0, 0, 127}));
  connect(PI_Pressure.PV, sDp.oDp) annotation(
    Line(points = {{-48, 50}, {-52, 50}, {-52, -38.4}}, color = {0, 0, 127}));
  connect(PI_Pressure.CS, pump.cmd) annotation(
    Line(points = {{-24, 54}, {-24, 15.5}, {-34, 15.5}, {-34, 4}}, color = {0, 0, 127}));
  connect(gain.y, Qheat.Q) annotation(
    Line(points = {{-207.6, -40}, {-187.6, -40}}, color = {0, 0, 127}));
  connect(Tamb.y[1], pTa.T) annotation(
    Line(points = {{-113, 146}, {-14, 146}}, color = {0, 0, 127}));
  connect(P_Loss.y, E_loss.u) annotation(
    Line(points = {{-380.1, 138}, {-313.1, 138}}, color = {0, 0, 127}));
  connect(realExpression.y, PI_Heater.SP) annotation(
    Line(points = {{-417.5, -34}, {-316, -34}}, color = {0, 0, 127}));
  connect(PI_Heater.CS, gain.u) annotation(
    Line(points = {{-292, -40}, {-216, -40}}, color = {0, 0, 127}));
  connect(sTh.oT, PI_Heater.PV) annotation(
    Line(points = {{-174, -4}, {-372, -4}, {-372, -44}, {-316, -44}}, color = {0, 0, 127}));
 connect(Psupz2.y, Hsupz2.Q_flow) annotation(
    Line(points = {{358, 58}, {358, 62}}, color = {0, 0, 127}));
 connect(Psupz1.y, Hsupz1.Q_flow) annotation(
    Line(points = {{122, 56}, {122, 60}}, color = {0, 0, 127}));
 connect(sp_Tz1.y[1], ref_z1.u1) annotation(
    Line(points = {{-17, 104}, {7, 104}, {7, 56}, {8, 56}}, color = {0, 0, 127}));
 connect(sTz1.T, ref_z1.u2) annotation(
    Line(points = {{80, 86}, {-4, 86}, {-4, 56}}, color = {0, 0, 127}));
 connect(ref_z2.y, integrator_z2.u) annotation(
    Line(points = {{238, 35}, {238, 24.5}, {251, 24.5}, {251, 26}}, color = {0, 0, 127}));
 connect(sTz2.T, ref_z2.u2) annotation(
    Line(points = {{282, 88}, {232, 88}, {232, 58}}, color = {0, 0, 127}));
 connect(ref_z1.y, integrator_z1.u) annotation(
    Line(points = {{2, 34}, {2, 23.5}, {9, 23.5}, {9, 24}}, color = {0, 0, 127}));
 connect(sp_Tz2.y[1], ref_z2.u1) annotation(
    Line(points = {{230, 118}, {244, 118}, {244, 58}}, color = {0, 0, 127}));
 connect(integrator_z2.y, mux_z2.u2[1]) annotation(
    Line(points = {{264, 26}, {270, 26}, {270, 38}, {276, 38}}, color = {0, 0, 127}));
 connect(sTz2.T, mux_z2.u1[1]) annotation(
    Line(points = {{282, 88}, {266, 88}, {266, 50}, {276, 50}}, color = {0, 0, 127}));
 connect(integrator_z1.y, mux_z1.u2[1]) annotation(
    Line(points = {{22, 24}, {38, 24}, {38, 36}, {42, 36}}, color = {0, 0, 127}));
 connect(sTz1.T, mux_z1.u1[1]) annotation(
    Line(points = {{80, 86}, {28, 86}, {28, 48}, {42, 48}}, color = {0, 0, 127}));
 connect(mux_z1.y, matrixGain_z1.u) annotation(
    Line(points = {{66, 42}, {72, 42}}, color = {0, 0, 127}, thickness = 0.5));
 connect(mux_z2.y, matrixGain_z2.u) annotation(
    Line(points = {{300, 44}, {308, 44}}, color = {0, 0, 127}, thickness = 0.5));
 connect(matrixGain_z1.y[2], sat2_z1.u) annotation(
    Line(points = {{96, 42}, {104, 42}}, color = {0, 0, 127}));
 connect(sat2_z1.y, Psupz1.u) annotation(
    Line(points = {{112, 42}, {122, 42}, {122, 48}}, color = {0, 0, 127}));
 connect(matrixGain_z1.y[1], sat1_z1.u) annotation(
    Line(points = {{96, 42}, {98, 42}, {98, 26}, {103, 26}}, color = {0, 0, 127}));
 connect(sat1_z1.y, vh1.x) annotation(
    Line(points = {{112, 26}, {132, 26}}, color = {0, 0, 127}));
 connect(matrixGain_z2.y[2], sat2_z2.u) annotation(
    Line(points = {{332, 44}, {340, 44}}, color = {0, 0, 127}));
 connect(sat2_z2.y, Psupz2.u) annotation(
    Line(points = {{348, 44}, {358, 44}, {358, 50}}, color = {0, 0, 127}));
 connect(matrixGain_z2.y[1], sat1_z2.u) annotation(
    Line(points = {{332, 44}, {334, 44}, {334, 26}, {342, 26}}, color = {0, 0, 127}));
 connect(sat1_z2.y, vh2.x) annotation(
    Line(points = {{350, 26}, {376, 26}}, color = {0, 0, 127}));
  annotation(
    Diagram(coordinateSystem(extent = {{-500, -200}, {500, 200}})),
    experiment(StartTime = 0, StopTime = 864000, Tolerance = 1e-6, Interval = 86.4),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian",
    __OpenModelica_simulationFlags(lv = "LOG_STATS", s = "dassl"));
end Plant_heater_with_control_1_LQ;
