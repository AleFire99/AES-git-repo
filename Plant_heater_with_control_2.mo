within AES_project_2021_2022;

model Plant_heater_with_control_2
  /* Simulation starts at hs:ms:ss*/
  parameter Real hs = 8;
  parameter Real ms = 0;
  parameter Real ss = 0;
  Real P_heaters = Qheat.Q + Hsupz1.Q_flow + Hsupz2.Q_flow;
  Real P_ambient = Gloss1.G * (sTz1.T - pTa.T) + thermalConductor.G * (sTz2.T - pTa.T);
  Real P_pump = pump.pwh_a.w * (pump.pwh_b.p - pump.pwh_a.p) / system.ro;
  Real P_loss = P_heaters + P_ambient + P_pump;
  AES.ProcessComponents.Thermal.Liquid.Pressuriser pressuriser annotation(
    Placement(visible = true, transformation(origin = {-32, -86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.Tube tubeh1(L = 50) annotation(
    Placement(visible = true, transformation(origin = {18, -30}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.Tube tubeh2(L = 50) annotation(
    Placement(visible = true, transformation(origin = {232, -30}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.Tube tubec1(L = 50) annotation(
    Placement(visible = true, transformation(origin = {250, -92}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.Tube tubec2(L = 50) annotation(
    Placement(visible = true, transformation(origin = {54, -92}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.Valve_linear vh1(dpnom = 300000, wnom = 0.2) annotation(
    Placement(visible = true, transformation(origin = {142, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  AES.ProcessComponents.Thermal.Liquid.Valve_linear vh2(dpnom = 300000, wnom = 2.2) annotation(
    Placement(visible = true, transformation(origin = {386, 22}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  AES.ProcessComponents.Thermal.Liquid.Tube h1(Di = 0.02, L = 5) annotation(
    Placement(visible = true, transformation(origin = {162, 46}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.Tube h2(L = 50) annotation(
    Placement(visible = true, transformation(origin = {458, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.Pump_centrifugal pump(dp0 = 600000, w0 = 30) annotation(
    Placement(visible = true, transformation(origin = {-32, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.Tube closure(L = 1000) annotation(
    Placement(visible = true, transformation(origin = {514, -64}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  AES.ProcessComponents.Thermal.Liquid.Tube tubeheat annotation(
    Placement(visible = true, transformation(origin = {-156, -48}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  AES.ProcessComponents.Thermal.Liquid.DiffPressureSensor sDp annotation(
    Placement(visible = true, transformation(origin = {-32, -58}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.Tsensor sTh annotation(
    Placement(visible = true, transformation(origin = {-184, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.surfQcond_prescribed Qheat annotation(
    Placement(visible = true, transformation(origin = {-216, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor Cz1(C = 1e4, T(displayUnit = "K")) annotation(
    Placement(visible = true, transformation(origin = {204, 86}, extent = {{-20, -20}, {20, 20}}, rotation = -90)));
  Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor sTz1 annotation(
    Placement(visible = true, transformation(origin = {94, 86}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  Modelica.Thermal.HeatTransfer.Components.ThermalConductor Gloss1(G = 80) annotation(
    Placement(visible = true, transformation(origin = {162, 116}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow Hsupz1 annotation(
    Placement(visible = true, transformation(origin = {122, 66}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  Modelica.Blocks.Math.Gain Psupz1(k = 500) annotation(
    Placement(visible = true, transformation(origin = {118, 36}, extent = {{-4, -4}, {4, 4}}, rotation = 0)));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature pTa annotation(
    Placement(visible = true, transformation(origin = {-2, 158}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.VectorHPtoHP_conductor convz1 annotation(
    Placement(visible = true, transformation(origin = {162, 66}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow Hsupz2 annotation(
    Placement(visible = true, transformation(origin = {436, 78}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor sTz2 annotation(
    Placement(visible = true, transformation(origin = {374, 102}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor Cz2(C = 1e4, T(displayUnit = "K")) annotation(
    Placement(visible = true, transformation(origin = {492, 102}, extent = {{-20, -20}, {20, 20}}, rotation = -90)));
  AES.ProcessComponents.Thermal.Liquid.VectorHPtoHP_conductor convz2 annotation(
    Placement(visible = true, transformation(origin = {456, 68}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
  Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor(G = 80) annotation(
    Placement(visible = true, transformation(origin = {456, 132}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Modelica.Blocks.Math.Gain Psupz2(k = 500) annotation(
    Placement(visible = true, transformation(origin = {424, 56}, extent = {{-4, -4}, {4, 4}}, rotation = 0)));
  inner AES.ProcessComponents.Thermal.System_settings.System_liquid system(ro(displayUnit = "kg/m3")) annotation(
    Placement(visible = true, transformation(origin = {-226, 164}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.CombiTimeTable sp_Tz(extrapolation = Modelica.Blocks.Types.Extrapolation.Periodic, offset = {273.15}, smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, table = [0, 7; 7, 7; 8, 20; 17, 20; 20, 14; 22, 10; 24, 10], timeScale = 3600) annotation(
    Placement(visible = true, transformation(origin = {-170, 140}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression Pressure_difference(y = 2.027 * 10 ^ 5) annotation(
    Placement(visible = true, transformation(origin = {-233, 6}, extent = {{-25, -18}, {25, 18}}, rotation = 0)));
  AES.ControlBlocks.AnalogueControllers.PI_awfb_basic PI_z1(CSmax = 1, CSmin = 0, K = 0.07692, Ti = 76.92) annotation(
    Placement(visible = true, transformation(origin = {-120, 124}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ControlBlocks.AnalogueControllers.PI_awfb_basic PI_z2(CSmax = 1, CSmin = 0, K = 0.3846, Ti = 76.92) annotation(
    Placement(visible = true, transformation(origin = {262, 142}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ControlBlocks.ActuationSchemes.DaisyChain_uniform daisyChain_z1 annotation(
    Placement(visible = true, transformation(origin = {-82, 124}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ControlBlocks.ActuationSchemes.DaisyChain_uniform daisyChain_z2 annotation(
    Placement(visible = true, transformation(origin = {316, 142}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Gain gain(k = 10000) annotation(
    Placement(visible = true, transformation(origin = {-260, -40}, extent = {{-4, -4}, {4, 4}}, rotation = 0)));
  Modelica.Blocks.Sources.CombiTimeTable Tamb(extrapolation = Modelica.Blocks.Types.Extrapolation.Periodic, offset = {273.15}, smoothness = Modelica.Blocks.Types.Smoothness.ContinuousDerivative, table = [0, 0; 4, -2; 8, 8; 12, 10; 15, 10; 18, 3; 20, 1; 22, 0; 24, 0], timeScale = 3600) annotation(
    Placement(visible = true, transformation(origin = {-196, 158}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ControlBlocks.AnalogueControllers.PI_awfb_basic PI_Heater(CSmax = 1, CSmin = 0, K = 0.0047, Ti = 39.2699) annotation(
    Placement(visible = true, transformation(origin = {-312, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression Heater_T_Max(y = 45 + 273.15) annotation(
    Placement(visible = true, transformation(origin = {-443, -34}, extent = {{-25, -18}, {25, 18}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression P_Loss(y = P_loss) annotation(
    Placement(visible = true, transformation(origin = {-401, 138}, extent = {{-19, -10}, {19, 10}}, rotation = 0)));
  Modelica.Blocks.Continuous.Integrator E_loss(initType = Modelica.Blocks.Types.Init.NoInit, use_reset = false) annotation(
    Placement(visible = true, transformation(origin = {-300, 138}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ControlBlocks.AnalogueControllers.PI_awfb_basic night_PI_z2(CSmax = 10, CSmin = 0, K = 0.07692, Ti = 100) annotation(
    Placement(visible = true, transformation(origin = {268, 14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression LO_limit(y = 5 + 273.15) annotation(
    Placement(visible = true, transformation(origin = {-16, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Logical.Switch switch2 annotation(
    Placement(visible = true, transformation(origin = {350, 22}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ControlBlocks.AnalogueControllers.PI_awfb_basic night_PI_z1(CSmax = 1, CSmin = 0, K = 0.07692, Ti = 100) annotation(
    Placement(visible = true, transformation(origin = {-2, 12}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.CombiTimeTable hours_switch(extrapolation = Modelica.Blocks.Types.Extrapolation.Periodic, smoothness = Modelica.Blocks.Types.Smoothness.ConstantSegments, table = [0, 0; 8, 1; 22, 0; 24, 0], tableOnFile = false, timeScale = 3600)  annotation(
    Placement(visible = true, transformation(origin = {-110, 216}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Logical.Switch switch1 annotation(
    Placement(visible = true, transformation(origin = {60, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Logical.GreaterThreshold greaterThreshold(threshold = 0.5)  annotation(
    Placement(visible = true, transformation(origin = {-40, 216}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
  Modelica.Blocks.Logical.Switch switch annotation(
    Placement(visible = true, transformation(origin = {58, 116}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression realExpression(y = 0) annotation(
    Placement(visible = true, transformation(origin = {15, 104}, extent = {{-7, -6}, {7, 6}}, rotation = 0)));
  Modelica.Blocks.Logical.Switch switch3 annotation(
    Placement(visible = true, transformation(origin = {394, 134}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression realExpression1(y = 0) annotation(
    Placement(visible = true, transformation(origin = {349, 126}, extent = {{-7, -6}, {7, 6}}, rotation = 0)));
  AES.ControlBlocks.AnalogueControllers.PI_awfb_basic PI_Pressure(CSmax = 1, CSmin = 0, K = 0.001, Ti = 500) annotation(
    Placement(visible = true, transformation(origin = {-113, -1}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
equation
  connect(pump.pwh_b, tubeh1.pwh_a) annotation(
    Line(points = {{-20, -30}, {6, -30}}, color = {46, 52, 54}));
  connect(tubeh1.pwh_b, vh1.pwh_a) annotation(
    Line(points = {{30, -30}, {142, -30}, {142, 6}}, color = {46, 52, 54}));
  connect(vh1.pwh_b, h1.pwh_a) annotation(
    Line(points = {{142, 30}, {142, 46}, {150, 46}}, color = {46, 52, 54}));
  connect(h1.pwh_b, tubec2.pwh_a) annotation(
    Line(points = {{174, 46}, {186, 46}, {186, -92}, {66, -92}}, color = {46, 52, 54}));
  connect(tubeh1.pwh_b, tubeh2.pwh_a) annotation(
    Line(points = {{30, -30}, {220, -30}}, color = {46, 52, 54}));
  connect(tubeh2.pwh_b, vh2.pwh_a) annotation(
    Line(points = {{244, -30}, {386, -30}, {386, 10}}, color = {46, 52, 54}));
  connect(vh2.pwh_b, h2.pwh_a) annotation(
    Line(points = {{386, 34}, {386, 38}, {446, 38}}, color = {46, 52, 54}));
  connect(tubec1.pwh_b, tubec2.pwh_a) annotation(
    Line(points = {{238, -92}, {66, -92}}, color = {46, 52, 54}));
  connect(h2.pwh_b, tubec1.pwh_a) annotation(
    Line(points = {{470, 38}, {470, -92}, {262, -92}}, color = {46, 52, 54}));
  connect(pressuriser.pwh_b, tubec2.pwh_b) annotation(
    Line(points = {{-20, -92}, {42, -92}}, color = {46, 52, 54}));
  connect(tubeh2.pwh_b, closure.pwh_a) annotation(
    Line(points = {{244, -30}, {514, -30}, {514, -52}}, color = {46, 52, 54}));
  connect(closure.pwh_b, tubec1.pwh_a) annotation(
    Line(points = {{514, -76}, {514, -92}, {262, -92}}, color = {46, 52, 54}));
  connect(tubeheat.pwh_b, pump.pwh_a) annotation(
    Line(points = {{-156, -36}, {-156, -30}, {-44, -30}}, color = {46, 52, 54}));
  connect(pressuriser.pwh_a, tubeheat.pwh_a) annotation(
    Line(points = {{-44, -92}, {-156, -92}, {-156, -60}}, color = {46, 52, 54}));
  connect(pump.pwh_b, sDp.pwh_hi) annotation(
    Line(points = {{-20, -30}, {-20, -52}}, color = {46, 52, 54}));
  connect(sDp.pwh_lo, pressuriser.pwh_b) annotation(
    Line(points = {{-20, -64}, {-20, -92}}, color = {46, 52, 54}));
  connect(sTh.pwh_a, tubeheat.pwh_b) annotation(
    Line(points = {{-172, -30}, {-156, -30}, {-156, -36}}, color = {46, 52, 54}));
  connect(Qheat.surf, tubeheat.surf) annotation(
    Line(points = {{-216, -49}, {-162, -49}, {-162, -48}, {-161, -48}}, color = {144, 5, 5}));
  connect(Cz1.port, Gloss1.port_b) annotation(
    Line(points = {{184, 86}, {162, 86}, {162, 106}}, color = {191, 0, 0}));
  connect(Psupz1.y, Hsupz1.Q_flow) annotation(
    Line(points = {{122, 36}, {122, 56.4}}, color = {0, 0, 127}));
  connect(Cz1.port, Hsupz1.port) annotation(
    Line(points = {{184, 86}, {122, 86}, {122, 76}}, color = {191, 0, 0}));
  connect(pTa.port, Gloss1.port_a) annotation(
    Line(points = {{8, 158}, {162, 158}, {162, 126}}, color = {191, 0, 0}));
  connect(Cz1.port, convz1.HP) annotation(
    Line(points = {{184, 86}, {162, 86}, {162, 70}}, color = {191, 0, 0}));
  connect(convz1.vectorHP, h1.surf) annotation(
    Line(points = {{162, 62}, {162, 52}}, color = {144, 5, 5}));
  connect(sTz1.port, Cz1.port) annotation(
    Line(points = {{104, 86}, {184, 86}}, color = {191, 0, 0}));
  connect(h2.surf, convz2.vectorHP) annotation(
    Line(points = {{458, 43}, {458, 53.5}, {456, 53.5}, {456, 64}}, color = {144, 5, 5}));
  connect(convz2.HP, Cz2.port) annotation(
    Line(points = {{456, 72}, {456, 102}, {472, 102}}, color = {191, 0, 0}));
  connect(Cz2.port, Hsupz2.port) annotation(
    Line(points = {{472, 102}, {436, 102}, {436, 88}}, color = {191, 0, 0}));
  connect(sTz2.port, Cz2.port) annotation(
    Line(points = {{384, 102}, {472, 102}}, color = {191, 0, 0}));
  connect(pTa.port, thermalConductor.port_a) annotation(
    Line(points = {{8, 158}, {456, 158}, {456, 142}}, color = {191, 0, 0}));
  connect(thermalConductor.port_b, Cz2.port) annotation(
    Line(points = {{456, 122}, {456, 102}, {472, 102}}, color = {191, 0, 0}));
  connect(PI_z1.CS, daisyChain_z1.CSi01) annotation(
    Line(points = {{-108, 124}, {-94, 124}}, color = {0, 0, 127}));
  connect(gain.y, Qheat.Q) annotation(
    Line(points = {{-256, -40}, {-224, -40}}, color = {0, 0, 127}));
  connect(sp_Tz.y[1], PI_z1.SP) annotation(
    Line(points = {{-159, 140}, {-155.5, 140}, {-155.5, 130}, {-132, 130}}, color = {0, 0, 127}));
  connect(Tamb.y[1], pTa.T) annotation(
    Line(points = {{-185, 158}, {-14, 158}}, color = {0, 0, 127}));
  connect(Heater_T_Max.y, PI_Heater.SP) annotation(
    Line(points = {{-415.5, -34}, {-324, -34}}, color = {0, 0, 127}));
  connect(P_Loss.y, E_loss.u) annotation(
    Line(points = {{-380.1, 138}, {-313.1, 138}}, color = {0, 0, 127}));
  connect(PI_Heater.CS, gain.u) annotation(
    Line(points = {{-300, -40}, {-265, -40}}, color = {0, 0, 127}));
  connect(sTh.oT, PI_Heater.PV) annotation(
    Line(points = {{-182, -30}, {-394, -30}, {-394, -44}, {-324, -44}}, color = {0, 0, 127}));
  connect(sTz2.T, night_PI_z2.PV) annotation(
    Line(points = {{364, 102}, {234, 102}, {234, 10}, {256, 10}}, color = {0, 0, 127}));
  connect(night_PI_z2.CS, switch2.u3) annotation(
    Line(points = {{280, 14}, {338, 14}}, color = {0, 0, 127}));
  connect(switch2.y, vh2.x) annotation(
    Line(points = {{361, 22}, {376, 22}}, color = {0, 0, 127}));
  connect(LO_limit.y, night_PI_z2.SP) annotation(
    Line(points = {{-5, 60}, {246, 60}, {246, 20}, {256, 20}}, color = {0, 0, 127}));
  connect(sTz1.T, night_PI_z1.PV) annotation(
    Line(points = {{84, 86}, {-84, 86}, {-84, 8}, {-14, 8}}, color = {0, 0, 127}));
  connect(hours_switch.y[1], greaterThreshold.u) annotation(
    Line(points = {{-98, 216}, {-50, 216}}, color = {0, 0, 127}));
  connect(greaterThreshold.y, switch1.u2) annotation(
    Line(points = {{-32, 216}, {28, 216}, {28, 20}, {48, 20}}, color = {255, 0, 255}));
  connect(greaterThreshold.y, switch2.u2) annotation(
    Line(points = {{-32, 216}, {282, 216}, {282, 22}, {338, 22}}, color = {255, 0, 255}));
  connect(PI_z2.PV, sTz2.T) annotation(
    Line(points = {{250, 138}, {234, 138}, {234, 102}, {364, 102}}, color = {0, 0, 127}));
  connect(LO_limit.y, night_PI_z1.SP) annotation(
    Line(points = {{-5, 60}, {-74, 60}, {-74, 18}, {-14, 18}}, color = {0, 0, 127}));
  connect(night_PI_z1.CS, switch1.u3) annotation(
    Line(points = {{10, 12}, {48, 12}}, color = {0, 0, 127}));
  connect(switch1.y, vh1.x) annotation(
    Line(points = {{72, 20}, {132, 20}, {132, 18}}, color = {0, 0, 127}));
  connect(PI_z1.PV, sTz1.T) annotation(
    Line(points = {{-132, 120}, {-132, 86}, {84, 86}}, color = {0, 0, 127}));
  connect(daisyChain_z2.CSo01[1], switch2.u1) annotation(
    Line(points = {{328, 142}, {328, 30}, {338, 30}}, color = {0, 0, 127}));
  connect(PI_z2.CS, daisyChain_z2.CSi01) annotation(
    Line(points = {{274, 142}, {304, 142}}, color = {0, 0, 127}));
  connect(Psupz2.y, Hsupz2.Q_flow) annotation(
    Line(points = {{428, 56}, {436, 56}, {436, 68}}, color = {0, 0, 127}));
  connect(daisyChain_z1.CSo01[1], switch1.u1) annotation(
    Line(points = {{-70, 124}, {-61, 124}, {-61, 114}, {2, 114}, {2, 28}, {48, 28}}, color = {0, 0, 127}));
  connect(daisyChain_z1.CSo01[2], switch.u1) annotation(
    Line(points = {{-70, 124}, {46, 124}}, color = {0, 0, 127}));
  connect(switch.u2, greaterThreshold.y) annotation(
    Line(points = {{46, 116}, {28, 116}, {28, 216}, {-32, 216}}, color = {255, 0, 255}));
  connect(realExpression.y, switch.u3) annotation(
    Line(points = {{23, 104}, {46, 104}, {46, 108}}, color = {0, 0, 127}));
  connect(switch.y, Psupz1.u) annotation(
    Line(points = {{70, 116}, {78, 116}, {78, 36}, {114, 36}}, color = {0, 0, 127}));
  connect(sp_Tz.y[1], PI_z2.SP) annotation(
    Line(points = {{-159, 140}, {47.5, 140}, {47.5, 148}, {250, 148}}, color = {0, 0, 127}));
  connect(realExpression1.y, switch3.u3) annotation(
    Line(points = {{356, 126}, {382, 126}}, color = {0, 0, 127}));
  connect(greaterThreshold.y, switch3.u2) annotation(
    Line(points = {{-32, 216}, {378, 216}, {378, 134}, {382, 134}}, color = {255, 0, 255}));
  connect(daisyChain_z2.CSo01[2], switch3.u1) annotation(
    Line(points = {{328, 142}, {382, 142}}, color = {0, 0, 127}));
  connect(switch3.y, Psupz2.u) annotation(
    Line(points = {{406, 134}, {408, 134}, {408, 56}, {419, 56}}, color = {0, 0, 127}));
  connect(Pressure_difference.y, PI_Pressure.SP) annotation(
    Line(points = {{-206, 6}, {-126, 6}}, color = {0, 0, 127}));
  connect(sDp.oDp, PI_Pressure.PV) annotation(
    Line(points = {{-30, -58}, {-126, -58}, {-126, -6}}, color = {0, 0, 127}));
  connect(PI_Pressure.CS, pump.cmd) annotation(
    Line(points = {{-100, -1}, {-100, 0}, {-32, 0}, {-32, -22}}, color = {0, 0, 127}));
  annotation(
    Diagram(coordinateSystem(extent = {{-480, 180}, {500, -100}})),
    experiment(StartTime = 0, StopTime = 864000, Tolerance = 1e-6, Interval = 86.4),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian",
    __OpenModelica_simulationFlags(lv = "LOG_STATS", s = "dassl"));
end Plant_heater_with_control_2;
