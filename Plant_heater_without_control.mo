within AES_project_2021_2022;

model Plant_heater_without_control
  /* Simulation starts at hs:ms:ss*/
  parameter Real hs = 8;
  parameter Real ms = 0;
  parameter Real ss = 0;
  AES.ProcessComponents.Thermal.Liquid.Pressuriser pressuriser annotation(
    Placement(visible = true, transformation(origin = {-296, -62}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.Tube tubeh1(L = 50) annotation(
    Placement(visible = true, transformation(origin = {-178, -2}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.Tube tubeh2(L = 50) annotation(
    Placement(visible = true, transformation(origin = {42, -2}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.Tube tubec1(L = 50) annotation(
    Placement(visible = true, transformation(origin = {42, -68}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.Tube tubec2(L = 50) annotation(
    Placement(visible = true, transformation(origin = {-178, -68}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.Valve_linear vh1(dpnom = 300000, wnom = 0.2) annotation(
    Placement(visible = true, transformation(origin = {-48, 26}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  AES.ProcessComponents.Thermal.Liquid.Valve_linear vh2(dpnom = 300000, wnom = 2.2) annotation(
    Placement(visible = true, transformation(origin = {198, 28}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  AES.ProcessComponents.Thermal.Liquid.Tube h1(Di = 0.02, L = 5) annotation(
    Placement(visible = true, transformation(origin = {-26, 48}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.Tube h2(L = 50) annotation(
    Placement(visible = true, transformation(origin = {210, 48}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.Pump_centrifugal pump(dp0 = 600000, w0 = 30) annotation(
    Placement(visible = true, transformation(origin = {-222, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.Tube closure(L = 1000) annotation(
    Placement(visible = true, transformation(origin = {294, -34}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  AES.ProcessComponents.Thermal.Liquid.Tube tubeheat annotation(
    Placement(visible = true, transformation(origin = {-332, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  AES.ProcessComponents.Thermal.Liquid.DiffPressureSensor sDp annotation(
    Placement(visible = true, transformation(origin = {-242, -36}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.Tsensor sTh annotation(
    Placement(visible = true, transformation(origin = {-364, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.surfQcond_prescribed Qheat annotation(
    Placement(visible = true, transformation(origin = {-366, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor Cz1(C = 1e4, T(displayUnit = "K")) annotation(
    Placement(visible = true, transformation(origin = {4, 88}, extent = {{-20, -20}, {20, 20}}, rotation = -90)));
  Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor sTz1 annotation(
    Placement(visible = true, transformation(origin = {-100, 88}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  Modelica.Thermal.HeatTransfer.Components.ThermalConductor Gloss1(G = 80) annotation(
    Placement(visible = true, transformation(origin = {-26, 118}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow Hsupz1 annotation(
    Placement(visible = true, transformation(origin = {-66, 68}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  Modelica.Blocks.Math.Gain Psupz1(k = 500) annotation(
    Placement(visible = true, transformation(origin = {-66, 46}, extent = {{-4, -4}, {4, 4}}, rotation = 0)));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature pTa annotation(
    Placement(visible = true, transformation(origin = {-190, 148}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.VectorHPtoHP_conductor convz1 annotation(
    Placement(visible = true, transformation(origin = {-26, 68}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow Hsupz2 annotation(
    Placement(visible = true, transformation(origin = {170, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor sTz2 annotation(
    Placement(visible = true, transformation(origin = {104, 90}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor Cz2(C = 1e4, T(displayUnit = "K")) annotation(
    Placement(visible = true, transformation(origin = {258, 90}, extent = {{-20, -20}, {20, 20}}, rotation = -90)));
  AES.ProcessComponents.Thermal.Liquid.VectorHPtoHP_conductor convz2 annotation(
    Placement(visible = true, transformation(origin = {210, 70}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
  Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor(G = 80) annotation(
    Placement(visible = true, transformation(origin = {210, 118}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Modelica.Blocks.Math.Gain Psupz2(k = 500) annotation(
    Placement(visible = true, transformation(origin = {178, 46}, extent = {{-4, -4}, {4, 4}}, rotation = 90)));
  inner AES.ProcessComponents.Thermal.System_settings.System_liquid system annotation(
    Placement(visible = true, transformation(origin = {-370, 150}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.CombiTimeTable Tamb(extrapolation = Modelica.Blocks.Types.Extrapolation.Periodic, offset = {273.15}, smoothness = Modelica.Blocks.Types.Smoothness.ContinuousDerivative, table = [0, 0; 4, -2; 8, 8; 12, 10; 15, 10; 18, 3; 20, 1; 22, 0; 24, 0], timeScale = 3600) annotation(
    Placement(visible = true, transformation(origin = {-312, 148}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.CombiTimeTable sp_Tz1(extrapolation = Modelica.Blocks.Types.Extrapolation.Periodic, offset = {273.15}, smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, table = [0, 7; 7, 7; 8, 20; 17, 20; 20, 14; 22, 10; 24, 10], timeScale = 3600) annotation(
    Placement(visible = true, transformation(origin = {-356, 118}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression Heater_T_ref(y = 45 + 273.15)  annotation(
    Placement(visible = true, transformation(origin = {-409, -96}, extent = {{-25, -18}, {25, 18}}, rotation = 180)));
  AES.ControlBlocks.AnalogueControllers.PI_awfb_basic PI_Heater(CSmax = 1, CSmin = 0, K = 100)  annotation(
    Placement(visible = true, transformation(origin = {-436, -32}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ControlBlocks.AnalogueControllers.PI_awfb_basic PI_Pressure(CSmax = 1, CSmin = 0, K = 10000)  annotation(
    Placement(visible = true, transformation(origin = {-218, 54}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression Pressure_difference(y = 2 * 10 ^ 5)  annotation(
    Placement(visible = true, transformation(origin = {-407, 62}, extent = {{-25, -18}, {25, 18}}, rotation = 0)));
  Modelica.Blocks.Sources.CombiTimeTable sp_z2(extrapolation = Modelica.Blocks.Types.Extrapolation.Periodic, offset = {273.15}, smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, table = [0, 7; 7, 7; 8, 20; 17, 20; 20, 14; 22, 10; 24, 10], timeScale = 3600) annotation(
    Placement(visible = true, transformation(origin = {40, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ControlBlocks.AnalogueControllers.PI_awfb_basic PI_z2(CSmax = 2, CSmin = 0)  annotation(
    Placement(visible = true, transformation(origin = {118, 54}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ControlBlocks.ActuationSchemes.DaisyChain_uniform daisyChain_z2 annotation(
    Placement(visible = true, transformation(origin = {160, 14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Gain gain1(k = -9.9812)  annotation(
    Placement(visible = true, transformation(origin = {-130, 28}, extent = {{8, 8}, {-8, -8}}, rotation = 180)));
  Modelica.Blocks.Math.Gain gain2(k = -0.0012) annotation(
    Placement(visible = true, transformation(origin = {-133, 57}, extent = {{-9, -9}, {9, 9}}, rotation = 0)));
  Modelica.Blocks.Math.Add add2 annotation(
    Placement(visible = true, transformation(origin = {-104, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Gain gain(k = 10000) annotation(
    Placement(visible = true, transformation(origin = {-396, -30}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  Modelica.Blocks.Math.Add add3 annotation(
    Placement(visible = true, transformation(origin = {-90, 22}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
equation
  connect(pump.pwh_b, tubeh1.pwh_a) annotation(
    Line(points = {{-210, -2}, {-190, -2}}, color = {46, 52, 54}));
  connect(tubeh1.pwh_b, vh1.pwh_a) annotation(
    Line(points = {{-166, -2}, {-48, -2}, {-48, 14}}, color = {46, 52, 54}));
  connect(vh1.pwh_b, h1.pwh_a) annotation(
    Line(points = {{-48, 38}, {-48, 48}, {-38, 48}}, color = {46, 52, 54}));
  connect(h1.pwh_b, tubec2.pwh_a) annotation(
    Line(points = {{-14, 48}, {-2, 48}, {-2, -68}, {-166, -68}}, color = {46, 52, 54}));
  connect(tubeh1.pwh_b, tubeh2.pwh_a) annotation(
    Line(points = {{-166, -2}, {30, -2}}, color = {46, 52, 54}));
  connect(tubeh2.pwh_b, vh2.pwh_a) annotation(
    Line(points = {{54, -2}, {198, -2}, {198, 16}}, color = {46, 52, 54}));
  connect(vh2.pwh_b, h2.pwh_a) annotation(
    Line(points = {{198, 40}, {198, 48}}, color = {46, 52, 54}));
  connect(tubec1.pwh_b, tubec2.pwh_a) annotation(
    Line(points = {{30, -68}, {-166, -68}}, color = {46, 52, 54}));
  connect(h2.pwh_b, tubec1.pwh_a) annotation(
    Line(points = {{222, 48}, {238, 48}, {238, -68}, {54, -68}}, color = {46, 52, 54}));
  connect(pressuriser.pwh_b, tubec2.pwh_b) annotation(
    Line(points = {{-284, -68}, {-190, -68}}, color = {46, 52, 54}));
  connect(tubeh2.pwh_b, closure.pwh_a) annotation(
    Line(points = {{54, -2}, {294, -2}, {294, -22}}, color = {46, 52, 54}));
  connect(closure.pwh_b, tubec1.pwh_a) annotation(
    Line(points = {{294, -46}, {294, -68}, {54, -68}}, color = {46, 52, 54}));
  connect(tubeheat.pwh_b, pump.pwh_a) annotation(
    Line(points = {{-332, -18}, {-332, -2}, {-234, -2}}, color = {46, 52, 54}));
  connect(pressuriser.pwh_a, tubeheat.pwh_a) annotation(
    Line(points = {{-308, -68}, {-332, -68}, {-332, -42}}, color = {46, 52, 54}));
  connect(pump.pwh_b, sDp.pwh_hi) annotation(
    Line(points = {{-210, -2}, {-210, -17}, {-230, -17}, {-230, -30}}, color = {46, 52, 54}));
  connect(sDp.pwh_lo, pressuriser.pwh_b) annotation(
    Line(points = {{-230, -42}, {-230, -68}, {-284, -68}}, color = {46, 52, 54}));
  connect(sTh.pwh_a, tubeheat.pwh_b) annotation(
    Line(points = {{-352, -2}, {-332, -2}, {-332, -18}}, color = {46, 52, 54}));
  connect(Qheat.surf, tubeheat.surf) annotation(
    Line(points = {{-366, -39}, {-350, -39}, {-350, -30.6667}, {-338, -30.6667}}, color = {144, 5, 5}));
  connect(Cz1.port, Gloss1.port_b) annotation(
    Line(points = {{-16, 88}, {-26, 88}, {-26, 108}}, color = {191, 0, 0}));
  connect(Psupz1.y, Hsupz1.Q_flow) annotation(
    Line(points = {{-62, 46}, {-62, 52.2}, {-66, 52.2}, {-66, 58.4}}, color = {0, 0, 127}));
  connect(Cz1.port, Hsupz1.port) annotation(
    Line(points = {{-16, 88}, {-66, 88}, {-66, 78}}, color = {191, 0, 0}));
  connect(pTa.port, Gloss1.port_a) annotation(
    Line(points = {{-180, 148}, {-26, 148}, {-26, 128}}, color = {191, 0, 0}));
  connect(Cz1.port, convz1.HP) annotation(
    Line(points = {{-16, 88}, {-26, 88}, {-26, 72}}, color = {191, 0, 0}));
  connect(convz1.vectorHP, h1.surf) annotation(
    Line(points = {{-26, 64}, {-26, 54}}, color = {144, 5, 5}));
  connect(sTz1.port, Cz1.port) annotation(
    Line(points = {{-90, 88}, {-16, 88}}, color = {191, 0, 0}));
  connect(Psupz2.y, Hsupz2.Q_flow) annotation(
    Line(points = {{178, 50}, {178, 54}, {170, 54}, {170, 60}}, color = {0, 0, 127}));
  connect(h2.surf, convz2.vectorHP) annotation(
    Line(points = {{210, 54}, {210, 66}}, color = {144, 5, 5}));
  connect(convz2.HP, Cz2.port) annotation(
    Line(points = {{210, 74}, {210, 90}, {238, 90}}, color = {191, 0, 0}));
  connect(Cz2.port, Hsupz2.port) annotation(
    Line(points = {{238, 90}, {170, 90}, {170, 80}}, color = {191, 0, 0}));
  connect(sTz2.port, Cz2.port) annotation(
    Line(points = {{114, 90}, {238, 90}}, color = {191, 0, 0}));
  connect(pTa.port, thermalConductor.port_a) annotation(
    Line(points = {{-180, 148}, {210, 148}, {210, 128}}, color = {191, 0, 0}));
  connect(thermalConductor.port_b, Cz2.port) annotation(
    Line(points = {{210, 108}, {210, 90}, {238, 90}}, color = {191, 0, 0}));
  connect(Tamb.y[1], pTa.T) annotation(
    Line(points = {{-301, 148}, {-202, 148}}, color = {0, 0, 127}));
  connect(Pressure_difference.y, PI_Pressure.SP) annotation(
    Line(points = {{-379.5, 62}, {-285.75, 62}, {-285.75, 60}, {-230, 60}}, color = {0, 0, 127}));
  connect(sp_z2.y[1], PI_z2.SP) annotation(
    Line(points = {{51, 60}, {106, 60}}, color = {0, 0, 127}));
  connect(PI_z2.CS, daisyChain_z2.CSi01) annotation(
    Line(points = {{130, 54}, {130, 13}, {148, 13}, {148, 14}}, color = {0, 0, 127}));
  connect(daisyChain_z2.CSo01[2], vh2.x) annotation(
    Line(points = {{172, 14}, {172, 21}, {188, 21}, {188, 28}}, color = {0, 0, 127}));
  connect(PI_Pressure.PV, sDp.oDp) annotation(
    Line(points = {{-230, 50}, {-230, 6}, {-240, 6}, {-240, -36}}, color = {0, 0, 127}));
  connect(sTz2.T, PI_z2.PV) annotation(
    Line(points = {{94, 90}, {94, 49}, {106, 49}, {106, 50}}, color = {0, 0, 127}));
  connect(sTh.oT, PI_Heater.PV) annotation(
    Line(points = {{-362, -2}, {-456, -2}, {-456, -36}, {-448, -36}}, color = {0, 0, 127}));
  connect(Heater_T_ref.y, PI_Heater.SP) annotation(
    Line(points = {{-436.5, -96}, {-461.75, -96}, {-461.75, -26}, {-448, -26}}, color = {0, 0, 127}));
  connect(daisyChain_z2.CSo01[1], Psupz2.u) annotation(
    Line(points = {{172, 14}, {172, 25.5}, {178, 25.5}, {178, 41}}, color = {0, 0, 127}));
  connect(PI_Pressure.CS, pump.cmd) annotation(
    Line(points = {{-206, 54}, {-206, 29}, {-222, 29}, {-222, 6}}, color = {0, 0, 127}));
  connect(PI_Heater.CS, gain.u) annotation(
    Line(points = {{-424, -32}, {-413.5, -32}, {-413.5, -30}, {-403, -30}}, color = {0, 0, 127}));
  connect(gain.y, Qheat.Q) annotation(
    Line(points = {{-390, -30}, {-374, -30}}, color = {0, 0, 127}));
  connect(add2.y, Psupz1.u) annotation(
    Line(points = {{-93, 50}, {-82, 50}, {-82, 46}, {-70, 46}}, color = {0, 0, 127}));
  connect(sTz1.T, gain2.u) annotation(
    Line(points = {{-110, 88}, {-160, 88}, {-160, 57}, {-144, 57}}, color = {0, 0, 127}));
  connect(gain2.y, add2.u1) annotation(
    Line(points = {{-123, 57}, {-120, 57}, {-120, 56}, {-116, 56}}, color = {0, 0, 127}));
  connect(sp_Tz1.y[1], add2.u2) annotation(
    Line(points = {{-344, 118}, {-174, 118}, {-174, 44}, {-116, 44}}, color = {0, 0, 127}));
  connect(add3.y, vh1.x) annotation(
    Line(points = {{-78, 22}, {-72, 22}, {-72, 26}, {-58, 26}}, color = {0, 0, 127}));
  connect(gain1.y, add3.u1) annotation(
    Line(points = {{-122, 28}, {-102, 28}}, color = {0, 0, 127}));
  connect(sTz1.T, gain1.u) annotation(
    Line(points = {{-110, 88}, {-160, 88}, {-160, 28}, {-140, 28}}, color = {0, 0, 127}));
  connect(sp_Tz1.y[1], add3.u2) annotation(
    Line(points = {{-344, 118}, {-174, 118}, {-174, 16}, {-102, 16}}, color = {0, 0, 127}));
  annotation(
    Diagram(coordinateSystem(extent = {{-500, -200}, {500, 200}})),
    experiment(StartTime = 0, StopTime = 864000, Tolerance = 1e-6, Interval = 86.4),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian",
    __OpenModelica_simulationFlags(lv = "LOG_STATS", s = "dassl"));
end Plant_heater_without_control;
