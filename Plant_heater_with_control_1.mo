within AES_project_2021_2022;

model Plant_heater_with_control_1
  /* Simulation starts at hs:ms:ss*/
  parameter Real hs = 8;
  parameter Real ms = 0;
  parameter Real ss = 0;
  Real P_heater = Qheat.Q;
  Real P_zones = Hsupz1.Q_flow + Hsupz2.Q_flow;
  Real P_ambient = Gloss1.G * (sTz1.T - pTa.T) + thermalConductor.G * (sTz2.T - pTa.T);
  Real P_pump = pump.pwh_a.w * (pump.pwh_b.p - pump.pwh_a.p) / system.ro;
  Real P_tot = P_heater + P_pump + P_zones;
  Real P_average = E_tot.y/time;
  
 
  
  Real eta = P_tot/(P_tot+P_ambient);
  Real eta_avg = eta/time;
  
  AES.ProcessComponents.Thermal.Liquid.Pressuriser pressuriser annotation(
    Placement(visible = true, transformation(origin = {-64, -64}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
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
    Placement(visible = true, transformation(origin = {398, 48}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.Pump_centrifugal pump(dp0 = 600000, w0 = 30) annotation(
    Placement(visible = true, transformation(origin = {-34, -4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.Tube closure(L = 1000) annotation(
    Placement(visible = true, transformation(origin = {482, -36}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  AES.ProcessComponents.Thermal.Liquid.Tube tubeheat annotation(
    Placement(visible = true, transformation(origin = {-100, -32}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  AES.ProcessComponents.Thermal.Liquid.DiffPressureSensor sDp annotation(
    Placement(visible = true, transformation(origin = {-54, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.Tsensor sTh annotation(
    Placement(visible = true, transformation(origin = {-132, -4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.surfQcond_prescribed Qheat annotation(
    Placement(visible = true, transformation(origin = {-136, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor Cz1(C = 1e4, T(displayUnit = "K")) annotation(
    Placement(visible = true, transformation(origin = {192, 86}, extent = {{-20, -20}, {20, 20}}, rotation = -90)));
  Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor sTz1 annotation(
    Placement(visible = true, transformation(origin = {90, 86}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  Modelica.Thermal.HeatTransfer.Components.ThermalConductor Gloss1(G = 80) annotation(
    Placement(visible = true, transformation(origin = {162, 116}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow Hsupz1 annotation(
    Placement(visible = true, transformation(origin = {122, 66}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  Modelica.Blocks.Math.Gain Psupz1(k = 500) annotation(
    Placement(visible = true, transformation(origin = {108, 42}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature pTa annotation(
    Placement(visible = true, transformation(origin = {-2, 146}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.VectorHPtoHP_conductor convz1 annotation(
    Placement(visible = true, transformation(origin = {162, 66}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow Hsupz2 annotation(
    Placement(visible = true, transformation(origin = {366, 68}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor sTz2 annotation(
    Placement(visible = true, transformation(origin = {292, 88}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor Cz2(C = 1e4, T(displayUnit = "K")) annotation(
    Placement(visible = true, transformation(origin = {446, 88}, extent = {{-20, -20}, {20, 20}}, rotation = -90)));
  AES.ProcessComponents.Thermal.Liquid.VectorHPtoHP_conductor convz2 annotation(
    Placement(visible = true, transformation(origin = {398, 68}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
  Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor(G = 80) annotation(
    Placement(visible = true, transformation(origin = {398, 116}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Modelica.Blocks.Math.Gain Psupz2(k = 500) annotation(
    Placement(visible = true, transformation(origin = {366, 42}, extent = {{-6, -6}, {6, 6}}, rotation = 90)));
  inner AES.ProcessComponents.Thermal.System_settings.System_liquid system(ro(displayUnit = "kg/m3")) annotation(
    Placement(visible = true, transformation(origin = {-168, 150}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.CombiTimeTable sp_Tz(extrapolation = Modelica.Blocks.Types.Extrapolation.Periodic, offset = {273.15}, smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, table = [0, 7; 7, 7; 8, 20; 17, 20; 20, 14; 22, 10; 24, 10], timeScale = 3600) annotation(
    Placement(visible = true, transformation(origin = {-168, 116}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ControlBlocks.AnalogueControllers.PI_awfb_basic PI_Pressure(CSmax = 1, CSmin = 0, K = 0.001, Ti = 1000) annotation(
    Placement(visible = true, transformation(origin = {6, 54}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression Pressure_difference(y = 2.027 * 10 ^ 5) annotation(
    Placement(visible = true, transformation(origin = {-219, 60}, extent = {{-25, -18}, {25, 18}}, rotation = 0)));
  AES.ControlBlocks.AnalogueControllers.PI_awfb_basic PI_z1(CSmax = 1, CSmin = 0, K = 0.0317, Ti = 76.9231) annotation(
    Placement(visible = true, transformation(origin = {-6, 110}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ControlBlocks.AnalogueControllers.PI_awfb_basic PI_z2(CSmax = 1, CSmin = 0, K = 0.1587, Ti = 76.92) annotation(
    Placement(visible = true, transformation(origin = {298, 52}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ControlBlocks.ActuationSchemes.DaisyChain_uniform daisyChain_z1 annotation(
    Placement(visible = true, transformation(origin = {49, 111}, extent = {{-15, -13}, {15, 13}}, rotation = 0)));
  AES.ControlBlocks.ActuationSchemes.DaisyChain_uniform daisyChain_z2 annotation(
    Placement(visible = true, transformation(origin = {341, 26}, extent = {{-15, -14}, {15, 14}}, rotation = 0)));
  Modelica.Blocks.Math.Gain gain(k = 10000) annotation(
    Placement(visible = true, transformation(origin = {-177, -41}, extent = {{-9, -9}, {9, 9}}, rotation = 0)));
  Modelica.Blocks.Sources.CombiTimeTable Tamb(extrapolation = Modelica.Blocks.Types.Extrapolation.Periodic, offset = {273.15}, smoothness = Modelica.Blocks.Types.Smoothness.ContinuousDerivative, table = [0, 0; 4, -2; 8, 8; 12, 10; 15, 10; 18, 3; 20, 1; 22, 0; 24, 0], timeScale = 3600) annotation(
    Placement(visible = true, transformation(origin = {-102, 146}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ControlBlocks.AnalogueControllers.PI_awfb_basic PI_Heater(CSmax = 1, CSmin = 0, K = 8.2877, Ti = 1.1446 * 10 ^ 3) annotation(
    Placement(visible = true, transformation(origin = {-224, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression Heater_T_Max(y = 40 + 273.15) annotation(
    Placement(visible = true, transformation(origin = {-291, -34}, extent = {{-25, -18}, {25, 18}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression P_tot_loss(y = P_tot) annotation(
    Placement(visible = true, transformation(origin = {-275, 114}, extent = {{-19, -10}, {19, 10}}, rotation = 0)));
  Modelica.Blocks.Continuous.Integrator E_tot(initType = Modelica.Blocks.Types.Init.NoInit, use_reset = false) annotation(
    Placement(visible = true, transformation(origin = {-212, 114}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
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
    Line(points = {{386, 38}, {386, 48}}, color = {46, 52, 54}));
  connect(tubec1.pwh_b, tubec2.pwh_a) annotation(
    Line(points = {{218, -70}, {22, -70}}, color = {46, 52, 54}));
  connect(h2.pwh_b, tubec1.pwh_a) annotation(
    Line(points = {{410, 48}, {426, 48}, {426, -70}, {242, -70}}, color = {46, 52, 54}));
  connect(pressuriser.pwh_b, tubec2.pwh_b) annotation(
    Line(points = {{-52, -70}, {-2, -70}}, color = {46, 52, 54}));
  connect(tubeh2.pwh_b, closure.pwh_a) annotation(
    Line(points = {{242, -4}, {482, -4}, {482, -24}}, color = {46, 52, 54}));
  connect(closure.pwh_b, tubec1.pwh_a) annotation(
    Line(points = {{482, -48}, {482, -70}, {242, -70}}, color = {46, 52, 54}));
  connect(tubeheat.pwh_b, pump.pwh_a) annotation(
    Line(points = {{-100, -20}, {-100, -4}, {-46, -4}}, color = {46, 52, 54}));
  connect(pressuriser.pwh_a, tubeheat.pwh_a) annotation(
    Line(points = {{-76, -70}, {-100, -70}, {-100, -44}}, color = {46, 52, 54}));
  connect(pump.pwh_b, sDp.pwh_hi) annotation(
    Line(points = {{-22, -4}, {-22, -19}, {-42, -19}, {-42, -32}}, color = {46, 52, 54}));
  connect(sDp.pwh_lo, pressuriser.pwh_b) annotation(
    Line(points = {{-42, -44}, {-42, -70}, {-52, -70}}, color = {46, 52, 54}));
  connect(sTh.pwh_a, tubeheat.pwh_b) annotation(
    Line(points = {{-120, -4}, {-100, -4}, {-100, -20}}, color = {46, 52, 54}));
  connect(Qheat.surf, tubeheat.surf) annotation(
    Line(points = {{-136, -48.6667}, {-118, -48.6667}, {-118, -32.3334}, {-106, -32.3334}}, color = {144, 5, 5}));
  connect(Cz1.port, Gloss1.port_b) annotation(
    Line(points = {{172, 86}, {162, 86}, {162, 106}}, color = {191, 0, 0}));
  connect(Cz1.port, Hsupz1.port) annotation(
    Line(points = {{172, 86}, {122, 86}, {122, 76}}, color = {191, 0, 0}));
  connect(pTa.port, Gloss1.port_a) annotation(
    Line(points = {{8, 146}, {162, 146}, {162, 126}}, color = {191, 0, 0}));
  connect(Cz1.port, convz1.HP) annotation(
    Line(points = {{172, 86}, {162, 86}, {162, 70}}, color = {191, 0, 0}));
  connect(convz1.vectorHP, h1.surf) annotation(
    Line(points = {{162, 62}, {162, 52}}, color = {144, 5, 5}));
  connect(sTz1.port, Cz1.port) annotation(
    Line(points = {{100, 86}, {172, 86}}, color = {191, 0, 0}));
  connect(Psupz2.y, Hsupz2.Q_flow) annotation(
    Line(points = {{366, 49}, {366, 58}}, color = {0, 0, 127}));
  connect(h2.surf, convz2.vectorHP) annotation(
    Line(points = {{398, 53}, {398, 63.4}}, color = {144, 5, 5}));
  connect(convz2.HP, Cz2.port) annotation(
    Line(points = {{398, 72}, {398, 88}, {426, 88}}, color = {191, 0, 0}));
  connect(Cz2.port, Hsupz2.port) annotation(
    Line(points = {{426, 88}, {366, 88}, {366, 78}}, color = {191, 0, 0}));
  connect(sTz2.port, Cz2.port) annotation(
    Line(points = {{302, 88}, {426, 88}}, color = {191, 0, 0}));
  connect(pTa.port, thermalConductor.port_a) annotation(
    Line(points = {{8, 146}, {398, 146}, {398, 126}}, color = {191, 0, 0}));
  connect(thermalConductor.port_b, Cz2.port) annotation(
    Line(points = {{398, 106}, {398, 88}, {426, 88}}, color = {191, 0, 0}));
  connect(Pressure_difference.y, PI_Pressure.SP) annotation(
    Line(points = {{-191.5, 60}, {-6, 60}}, color = {0, 0, 127}));
  connect(sTz1.T, PI_z1.PV) annotation(
    Line(points = {{80, 86}, {-66, 86}, {-66, 105}, {-18, 105}, {-18, 106}}, color = {0, 0, 127}));
  connect(sTz2.T, PI_z2.PV) annotation(
    Line(points = {{282, 88}, {282, 47}, {286, 47}, {286, 48}}, color = {0, 0, 127}));
  connect(PI_Pressure.CS, pump.cmd) annotation(
    Line(points = {{18, 54}, {18, 15.5}, {-34, 15.5}, {-34, 4}}, color = {0, 0, 127}));
  connect(sp_Tz.y[1], PI_z1.SP) annotation(
    Line(points = {{-157, 116}, {-18, 116}}, color = {0, 0, 127}));
  connect(sp_Tz.y[1], PI_z2.SP) annotation(
    Line(points = {{-157, 116}, {-65, 116}, {-65, 180}, {245, 180}, {245, 58}, {286, 58}}, color = {0, 0, 127}));
  connect(Heater_T_Max.y, PI_Heater.SP) annotation(
    Line(points = {{-263.5, -34}, {-236, -34}}, color = {0, 0, 127}));
  connect(P_tot_loss.y, E_tot.u) annotation(
    Line(points = {{-254.1, 114}, {-224.1, 114}}, color = {0, 0, 127}));
  connect(sTh.oT, PI_Heater.PV) annotation(
    Line(points = {{-130, -4}, {-252, -4}, {-252, -44}, {-236, -44}}, color = {0, 0, 127}));
  connect(daisyChain_z1.CSo01[1], vh1.x) annotation(
    Line(points = {{67, 111}, {74, 111}, {74, 26}, {132, 26}}, color = {0, 0, 127}));
  connect(daisyChain_z1.CSo01[2], Psupz1.u) annotation(
    Line(points = {{67, 111}, {74, 111}, {74, 42}, {98, 42}}, color = {0, 0, 127}));
  connect(Tamb.y[1], pTa.T) annotation(
    Line(points = {{-91, 146}, {-14, 146}}, color = {0, 0, 127}));
  connect(PI_z1.CS, daisyChain_z1.CSi01) annotation(
    Line(points = {{6, 110}, {30, 110}}, color = {0, 0, 127}));
  connect(daisyChain_z2.CSo01[1], vh2.x) annotation(
    Line(points = {{359, 26}, {376, 26}}, color = {0, 0, 127}));
  connect(daisyChain_z2.CSo01[2], Psupz2.u) annotation(
    Line(points = {{359, 26}, {366, 26}, {366, 35}}, color = {0, 0, 127}));
  connect(PI_z2.CS, daisyChain_z2.CSi01) annotation(
    Line(points = {{310, 52}, {314, 52}, {314, 26}, {324, 26}}, color = {0, 0, 127}));
  connect(gain.y, Qheat.Q) annotation(
    Line(points = {{-167.1, -41}, {-143.1, -41}}, color = {0, 0, 127}));
  connect(PI_Heater.CS, gain.u) annotation(
    Line(points = {{-212, -40}, {-188, -40}}, color = {0, 0, 127}));
  connect(sDp.oDp, PI_Pressure.PV) annotation(
    Line(points = {{-52, -38}, {-74, -38}, {-74, 50}, {-6, 50}}, color = {0, 0, 127}));
  connect(Psupz1.y, Hsupz1.Q_flow) annotation(
    Line(points = {{116, 42}, {122, 42}, {122, 56}}, color = {0, 0, 127}));
  annotation(
    Diagram(coordinateSystem(extent = {{-500, -200}, {500, 200}})),
    experiment(StartTime = 0, StopTime = 864000, Tolerance = 1e-6, Interval = 86.4),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian",
    __OpenModelica_simulationFlags(lv = "LOG_STATS", s = "dassl"));
end Plant_heater_with_control_1;
