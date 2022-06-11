within AES_project_2021_2022;

model Plant_heater_with_control_2
  /* Simulation starts at hs:ms:ss*/
  parameter Real hs = 8;
  parameter Real ms = 0;
  parameter Real ss = 0;
  Real P_heater = Qheat.Q;
  Real P_zones = Hsupz1.Q_flow + Hsupz2.Q_flow;
  Real P_env = Gloss1.G * (sTz1.T - pTa.T) + thermalConductor.G * (sTz2.T - pTa.T);
  Real P_pump = pump.pwh_a.w * (pump.pwh_b.p - pump.pwh_a.p) / system.ro;
  Real P_tot = P_heater + P_zones + P_pump;
  Real P_average = E_tot.y/time;
  
  Real eta = P_tot/(P_tot+P_env);
  
  AES.ProcessComponents.Thermal.Liquid.Pressuriser pressuriser annotation(
    Placement(visible = true, transformation(origin = {-24, -148}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.Tube tubeh1(L = 50) annotation(
    Placement(visible = true, transformation(origin = {26, -92}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.Tube tubeh2(L = 50) annotation(
    Placement(visible = true, transformation(origin = {288, -92}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.Tube tubec1(L = 50) annotation(
    Placement(visible = true, transformation(origin = {306, -154}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.Tube tubec2(L = 50) annotation(
    Placement(visible = true, transformation(origin = {62, -154}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.Valve_linear vh1(dpnom = 300000, wnom = 0.2) annotation(
    Placement(visible = true, transformation(origin = {142, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  AES.ProcessComponents.Thermal.Liquid.Valve_linear vh2(dpnom = 300000, wnom = 2.2) annotation(
    Placement(visible = true, transformation(origin = {460, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  AES.ProcessComponents.Thermal.Liquid.Tube h1(Di = 0.02, L = 5) annotation(
    Placement(visible = true, transformation(origin = {162, 32}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.Tube h2(L = 50) annotation(
    Placement(visible = true, transformation(origin = {530, 36}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.Pump_centrifugal pump(dp0 = 600000, w0 = 30) annotation(
    Placement(visible = true, transformation(origin = {-24, -92}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.Tube closure(L = 1000) annotation(
    Placement(visible = true, transformation(origin = {592, -126}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  AES.ProcessComponents.Thermal.Liquid.Tube tubeheat annotation(
    Placement(visible = true, transformation(origin = {-148, -114}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  AES.ProcessComponents.Thermal.Liquid.DiffPressureSensor sDp annotation(
    Placement(visible = true, transformation(origin = {-24, -120}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.Tsensor sTh annotation(
    Placement(visible = true, transformation(origin = {-176, -92}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.surfQcond_prescribed Qheat annotation(
    Placement(visible = true, transformation(origin = {-208, -106}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor Cz1(C = 1e4, T(displayUnit = "K", fixed = true)) annotation(
    Placement(visible = true, transformation(origin = {202, 86}, extent = {{-20, -20}, {20, 20}}, rotation = -90)));
  Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor sTz1 annotation(
    Placement(visible = true, transformation(origin = {94, 86}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  Modelica.Thermal.HeatTransfer.Components.ThermalConductor Gloss1(G = 80) annotation(
    Placement(visible = true, transformation(origin = {162, 116}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow Hsupz1 annotation(
    Placement(visible = true, transformation(origin = {132, 66}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  Modelica.Blocks.Math.Gain Psupz1(k = 500) annotation(
    Placement(visible = true, transformation(origin = {128, 30}, extent = {{-4, -4}, {4, 4}}, rotation = 0)));
  AES.ProcessComponents.Thermal.Liquid.VectorHPtoHP_conductor convz1 annotation(
    Placement(visible = true, transformation(origin = {162, 66}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow Hsupz2 annotation(
    Placement(visible = true, transformation(origin = {510, 76}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor sTz2 annotation(
    Placement(visible = true, transformation(origin = {448, 100}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor Cz2(C = 1e4, T(displayUnit = "K", fixed = true)) annotation(
    Placement(visible = true, transformation(origin = {612, 100}, extent = {{-20, -20}, {20, 20}}, rotation = -90)));
  AES.ProcessComponents.Thermal.Liquid.VectorHPtoHP_conductor convz2 annotation(
    Placement(visible = true, transformation(origin = {530, 66}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
  Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor(G = 80) annotation(
    Placement(visible = true, transformation(origin = {524, 130}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Modelica.Blocks.Math.Gain Psupz2(k = 500) annotation(
    Placement(visible = true, transformation(origin = {498, 54}, extent = {{-4, -4}, {4, 4}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression Pressure_difference(y = 2.027 * 10 ^ 5) annotation(
    Placement(visible = true, transformation(origin = {-341, -56}, extent = {{-25, -18}, {25, 18}}, rotation = 0)));
  AES.ControlBlocks.ActuationSchemes.DaisyChain_uniform daisyChain_z2 annotation(
    Placement(visible = true, transformation(origin = {414, 74}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Gain gain(k = 10000) annotation(
    Placement(visible = true, transformation(origin = {-250, -106}, extent = {{-4, -4}, {4, 4}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression Heater_T_Max(y = 45 + 273.15) annotation(
    Placement(visible = true, transformation(origin = {-389, -98}, extent = {{-25, -18}, {25, 18}}, rotation = 0)));
  Modelica.Blocks.Logical.Switch switch2 annotation(
    Placement(visible = true, transformation(origin = {380, 74}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.CombiTimeTable hours_switch(extrapolation = Modelica.Blocks.Types.Extrapolation.Periodic, smoothness = Modelica.Blocks.Types.Smoothness.ConstantSegments, table = [0, 0; 7.5, 1; 22, 0; 24, 0], tableOnFile = false, timeScale = 3600) annotation(
    Placement(visible = true, transformation(origin = {-232, 234}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Logical.GreaterThreshold greaterThreshold(threshold = 0.5) annotation(
    Placement(visible = true, transformation(origin = {-158, 234}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
  AES.ControlBlocks.AnalogueControllers.PI_awfb_basic PI_Pressure(CSmax = 1, CSmin = 0, K = 0.001, Ti = 500) annotation(
    Placement(visible = true, transformation(origin = {-99, -63}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
  Modelica.Blocks.Logical.Not not1 annotation(
    Placement(visible = true, transformation(origin = {-72, 216}, extent = {{-6, -6}, {6, 6}}, rotation = -90)));
  AES.ControlBlocks.AnalogueControllers.PI_awfb_full pI_z2D(CSmax = 2, CSmin = 0, K = 0.3846, Ti = 76.92, hasTracking = true) annotation(
    Placement(visible = true, transformation(origin = {322, 126}, extent = {{-10, -20}, {10, 20}}, rotation = 0)));
  AES.ControlBlocks.AnalogueControllers.PI_awfb_full pI_z2N(CSmax = 1, CSmin = 0, K = 0.3846, Ti = 76.92, hasTracking = true) annotation(
    Placement(visible = true, transformation(origin = {326, -4}, extent = {{-10, -20}, {10, 20}}, rotation = 0)));
  AES.ControlBlocks.ActuationSchemes.DaisyChain_uniform daisyChain_z1 annotation(
    Placement(visible = true, transformation(origin = {58, 64}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ControlBlocks.AnalogueControllers.PI_awfb_full pI_z1N(CSmax = 1, CSmin = 0, K = 0.07692, Ti = 76.92, hasTracking = true) annotation(
    Placement(visible = true, transformation(origin = {-36, 2}, extent = {{-10, -20}, {10, 20}}, rotation = 0)));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature pTa annotation(
    Placement(visible = true, transformation(origin = {8, 154}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Logical.Switch switch4 annotation(
    Placement(visible = true, transformation(origin = {12, 68}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression LO_limit(y = 5 + 273.15) annotation(
    Placement(visible = true, transformation(origin = {-58, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  AES.ControlBlocks.AnalogueControllers.PI_awfb_full pI_z1D(CSmax = 2, CSmin = 0, K = 0.07692, Ti = 76.92, hasTracking = true) annotation(
    Placement(visible = true, transformation(origin = {-40, 112}, extent = {{-10, -20}, {10, 20}}, rotation = 0)));
  Modelica.Blocks.Sources.CombiTimeTable sp_Tz(extrapolation = Modelica.Blocks.Types.Extrapolation.Periodic, offset = {273.15}, smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, table = [0, 7; 7, 7; 8, 20; 17, 20; 20, 14; 22, 10; 24, 10], timeScale = 3600) annotation(
    Placement(visible = true, transformation(origin = {-148, 138}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.CombiTimeTable Tamb(extrapolation = Modelica.Blocks.Types.Extrapolation.Periodic, offset = {273.15}, smoothness = Modelica.Blocks.Types.Smoothness.ContinuousDerivative, table = [0, 0; 4, -2; 8, 8; 12, 10; 15, 10; 18, 3; 20, 1; 22, 0; 24, 0], timeScale = 3600) annotation(
    Placement(visible = true, transformation(origin = {-186, 154}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  inner AES.ProcessComponents.Thermal.System_settings.System_liquid system(ro(displayUnit = "kg/m3")) annotation(
    Placement(visible = true, transformation(origin = {-216, 160}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Continuous.Integrator E_tot(initType = Modelica.Blocks.Types.Init.NoInit, use_reset = false) annotation(
    Placement(visible = true, transformation(origin = {-268, 106}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression P_Tot(y = P_tot) annotation(
    Placement(visible = true, transformation(origin = {-345, 106}, extent = {{-19, -10}, {19, 10}}, rotation = 0)));
  Modelica.Blocks.Logical.Not not2 annotation(
    Placement(visible = true, transformation(origin = {292, 208}, extent = {{-6, -6}, {6, 6}}, rotation = -90)));
  AES.ControlBlocks.AnalogueControllers.PI_awfb_basic pI_awfb_basic(CSmax = 1, CSmin = 0, K = 8.2877, Ti = 1.145 * 10 ^ 3) annotation(
    Placement(visible = true, transformation(origin = {-307, -105}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
equation
  connect(pump.pwh_b, tubeh1.pwh_a) annotation(
    Line(points = {{-12, -92}, {14, -92}}, color = {46, 52, 54}));
  connect(tubeh1.pwh_b, vh1.pwh_a) annotation(
    Line(points = {{38, -92}, {142, -92}, {142, -8}}, color = {46, 52, 54}));
  connect(vh1.pwh_b, h1.pwh_a) annotation(
    Line(points = {{142, 16}, {142, 32}, {150, 32}}, color = {46, 52, 54}));
  connect(h1.pwh_b, tubec2.pwh_a) annotation(
    Line(points = {{174, 32}, {186, 32}, {186, -154}, {74, -154}}, color = {46, 52, 54}));
  connect(tubeh1.pwh_b, tubeh2.pwh_a) annotation(
    Line(points = {{38, -92}, {276, -92}}, color = {46, 52, 54}));
  connect(tubeh2.pwh_b, vh2.pwh_a) annotation(
    Line(points = {{300, -92}, {460, -92}, {460, 8}}, color = {46, 52, 54}));
  connect(vh2.pwh_b, h2.pwh_a) annotation(
    Line(points = {{460, 32}, {460, 36}, {518, 36}}, color = {46, 52, 54}));
  connect(tubec1.pwh_b, tubec2.pwh_a) annotation(
    Line(points = {{294, -154}, {74, -154}}, color = {46, 52, 54}));
  connect(h2.pwh_b, tubec1.pwh_a) annotation(
    Line(points = {{542, 36}, {542, -154}, {318, -154}}, color = {46, 52, 54}));
  connect(pressuriser.pwh_b, tubec2.pwh_b) annotation(
    Line(points = {{-12, -154}, {50, -154}}, color = {46, 52, 54}));
  connect(tubeh2.pwh_b, closure.pwh_a) annotation(
    Line(points = {{300, -92}, {592, -92}, {592, -114}}, color = {46, 52, 54}));
  connect(closure.pwh_b, tubec1.pwh_a) annotation(
    Line(points = {{592, -138}, {592, -154}, {318, -154}}, color = {46, 52, 54}));
  connect(tubeheat.pwh_b, pump.pwh_a) annotation(
    Line(points = {{-148, -102}, {-148, -92}, {-36, -92}}, color = {46, 52, 54}));
  connect(pressuriser.pwh_a, tubeheat.pwh_a) annotation(
    Line(points = {{-36, -154}, {-148, -154}, {-148, -126}}, color = {46, 52, 54}));
  connect(pump.pwh_b, sDp.pwh_hi) annotation(
    Line(points = {{-12, -92}, {-12, -114}}, color = {46, 52, 54}));
  connect(sDp.pwh_lo, pressuriser.pwh_b) annotation(
    Line(points = {{-12, -126}, {-12, -154}}, color = {46, 52, 54}));
  connect(sTh.pwh_a, tubeheat.pwh_b) annotation(
    Line(points = {{-164, -92}, {-148, -92}, {-148, -102}}, color = {46, 52, 54}));
  connect(Qheat.surf, tubeheat.surf) annotation(
    Line(points = {{-208, -115}, {-154, -115}, {-154, -114}, {-153, -114}}, color = {144, 5, 5}));
  connect(Cz1.port, Gloss1.port_b) annotation(
    Line(points = {{182, 86}, {162, 86}, {162, 106}}, color = {191, 0, 0}));
  connect(Psupz1.y, Hsupz1.Q_flow) annotation(
    Line(points = {{132, 30}, {132, 56}}, color = {0, 0, 127}));
  connect(Cz1.port, Hsupz1.port) annotation(
    Line(points = {{182, 86}, {132, 86}, {132, 76}}, color = {191, 0, 0}));
  connect(Cz1.port, convz1.HP) annotation(
    Line(points = {{182, 86}, {162, 86}, {162, 70}}, color = {191, 0, 0}));
  connect(convz1.vectorHP, h1.surf) annotation(
    Line(points = {{162, 62}, {162, 37}}, color = {144, 5, 5}));
  connect(sTz1.port, Cz1.port) annotation(
    Line(points = {{104, 86}, {182, 86}}, color = {191, 0, 0}));
  connect(h2.surf, convz2.vectorHP) annotation(
    Line(points = {{530, 41}, {530, 62.4}}, color = {144, 5, 5}));
  connect(convz2.HP, Cz2.port) annotation(
    Line(points = {{530, 70}, {530, 100}, {592, 100}}, color = {191, 0, 0}));
  connect(Cz2.port, Hsupz2.port) annotation(
    Line(points = {{592, 100}, {510, 100}, {510, 86}}, color = {191, 0, 0}));
  connect(sTz2.port, Cz2.port) annotation(
    Line(points = {{458, 100}, {592, 100}}, color = {191, 0, 0}));
  connect(thermalConductor.port_b, Cz2.port) annotation(
    Line(points = {{524, 120}, {524, 100}, {592, 100}}, color = {191, 0, 0}));
  connect(gain.y, Qheat.Q) annotation(
    Line(points = {{-246, -106}, {-216, -106}}, color = {0, 0, 127}));
  connect(hours_switch.y[1], greaterThreshold.u) annotation(
    Line(points = {{-221, 234}, {-168, 234}}, color = {0, 0, 127}));
  connect(greaterThreshold.y, switch2.u2) annotation(
    Line(points = {{-149, 234}, {352, 234}, {352, 74}, {368, 74}}, color = {255, 0, 255}));
  connect(Psupz2.y, Hsupz2.Q_flow) annotation(
    Line(points = {{502.4, 54}, {510.4, 54}, {510.4, 66}}, color = {0, 0, 127}));
  connect(Pressure_difference.y, PI_Pressure.SP) annotation(
    Line(points = {{-313.5, -56}, {-112, -56}}, color = {0, 0, 127}));
  connect(sDp.oDp, PI_Pressure.PV) annotation(
    Line(points = {{-22, -120}, {-112, -120}, {-112, -68}}, color = {0, 0, 127}));
  connect(PI_Pressure.CS, pump.cmd) annotation(
    Line(points = {{-86, -63}, {-86, -62}, {-23.8, -62}, {-23.8, -84}}, color = {0, 0, 127}));
  connect(greaterThreshold.y, not1.u) annotation(
    Line(points = {{-149.2, 234}, {-72, 234}, {-72, 223}}, color = {255, 0, 255}));
  connect(switch2.y, daisyChain_z2.CSi01) annotation(
    Line(points = {{391, 74}, {402, 74}}, color = {0, 0, 127}));
  connect(sTz2.T, pI_z2D.PV) annotation(
    Line(points = {{438, 100}, {298, 100}, {298, 136}, {312, 136}}, color = {0, 0, 127}));
  connect(sTz2.T, pI_z2N.PV) annotation(
    Line(points = {{438, 100}, {298, 100}, {298, 6}, {316, 6}}, color = {0, 0, 127}));
  connect(pI_z2N.CS, switch2.u3) annotation(
    Line(points = {{336, 10}, {350, 10}, {350, 66}, {368, 66}}, color = {0, 0, 127}));
  connect(pI_z2D.CS, switch2.u1) annotation(
    Line(points = {{332, 140}, {348, 140}, {348, 82}, {368, 82}}, color = {0, 0, 127}));
  connect(daisyChain_z2.CSo01[2], Psupz2.u) annotation(
    Line(points = {{426, 74}, {466, 74}, {466, 54}, {493, 54}}, color = {0, 0, 127}));
  connect(daisyChain_z2.CSo01[1], vh2.x) annotation(
    Line(points = {{426, 74}, {444, 74}, {444, 20}, {450, 20}}, color = {0, 0, 127}));
  connect(daisyChain_z1.CSo01[2], Psupz1.u) annotation(
    Line(points = {{70, 64}, {70, 30}, {123, 30}}, color = {0, 0, 127}));
  connect(daisyChain_z1.CSo01[1], vh1.x) annotation(
    Line(points = {{70, 64}, {70, 4}, {132, 4}}, color = {0, 0, 127}));
  connect(pTa.port, thermalConductor.port_a) annotation(
    Line(points = {{18, 154}, {524, 154}, {524, 140}}, color = {191, 0, 0}));
  connect(pTa.port, Gloss1.port_a) annotation(
    Line(points = {{18, 154}, {162, 154}, {162, 126}}, color = {191, 0, 0}));
  connect(switch4.y, pI_z1N.TR) annotation(
    Line(points = {{23, 68}, {33, 68}, {33, -24}, {-159, -24}, {-159, 4}, {-46, 4}}, color = {0, 0, 127}));
  connect(switch4.y, daisyChain_z1.CSi01) annotation(
    Line(points = {{23, 68}, {32.5, 68}, {32.5, 64}, {46, 64}}, color = {0, 0, 127}));
  connect(greaterThreshold.y, switch4.u2) annotation(
    Line(points = {{-149, 234}, {0, 234}, {0, 68}}, color = {255, 0, 255}));
  connect(pI_z1N.CS, switch4.u3) annotation(
    Line(points = {{-26, 16}, {-16, 16}, {-16, 60}, {0, 60}}, color = {0, 0, 127}));
  connect(LO_limit.y, pI_z2N.SP) annotation(
    Line(points = {{-47, 50}, {232, 50}, {232, 10}, {316, 10}}, color = {0, 0, 127}));
  connect(LO_limit.y, pI_z1N.SP) annotation(
    Line(points = {{-47, 50}, {-111, 50}, {-111, 16}, {-46, 16}}, color = {0, 0, 127}));
  connect(sTz1.T, pI_z1N.PV) annotation(
    Line(points = {{84, 86}, {-134, 86}, {-134, 12}, {-46, 12}}, color = {0, 0, 127}));
  connect(not1.y, pI_z1D.TS) annotation(
    Line(points = {{-72, 209}, {-72, 118}, {-50, 118}}, color = {255, 0, 255}));
  connect(switch4.y, pI_z1D.TR) annotation(
    Line(points = {{23, 68}, {33, 68}, {33, -24}, {-159, -24}, {-159, 114}, {-50, 114}}, color = {0, 0, 127}));
  connect(pI_z1D.CS, switch4.u1) annotation(
    Line(points = {{-30, 126}, {-10, 126}, {-10, 76}, {0, 76}}, color = {0, 0, 127}));
  connect(sTz1.T, pI_z1D.PV) annotation(
    Line(points = {{84, 86}, {-134, 86}, {-134, 122}, {-50, 122}}, color = {0, 0, 127}));
  connect(sp_Tz.y[1], pI_z2D.SP) annotation(
    Line(points = {{-137, 138}, {312, 138}, {312, 140}}, color = {0, 0, 127}));
  connect(sp_Tz.y[1], pI_z1D.SP) annotation(
    Line(points = {{-137, 138}, {-119.5, 138}, {-119.5, 126}, {-50, 126}}, color = {0, 0, 127}));
  connect(Tamb.y[1], pTa.T) annotation(
    Line(points = {{-175, 154}, {-4, 154}}, color = {0, 0, 127}));
  connect(P_Tot.y, E_tot.u) annotation(
    Line(points = {{-324.1, 106}, {-280.1, 106}}, color = {0, 0, 127}));
  connect(sTz1.T, pI_z1N.PV) annotation(
    Line(points = {{84, 86}, {-134, 86}, {-134, 18}, {-86, 18}}, color = {0, 0, 127}));
  connect(greaterThreshold.y, pI_z1N.TS) annotation(
    Line(points = {{-150, 234}, {-98, 234}, {-98, 8}, {-46, 8}}, color = {255, 0, 255}));
  connect(greaterThreshold.y, not2.u) annotation(
    Line(points = {{-150, 234}, {292, 234}, {292, 215}}, color = {255, 0, 255}));
  connect(not2.y, pI_z2D.TS) annotation(
    Line(points = {{292, 201.4}, {292, 132.4}, {312, 132.4}}, color = {255, 0, 255}));
  connect(greaterThreshold.y, pI_z2N.TS) annotation(
    Line(points = {{-150, 234}, {282, 234}, {282, 2}, {316, 2}}, color = {255, 0, 255}));
  connect(switch2.y, pI_z2N.TR) annotation(
    Line(points = {{392, 74}, {396, 74}, {396, -48}, {268, -48}, {268, -2}, {316, -2}}, color = {0, 0, 127}));
  connect(switch2.y, pI_z2D.TR) annotation(
    Line(points = {{392, 74}, {396, 74}, {396, -48}, {268, -48}, {268, 128}, {312, 128}}, color = {0, 0, 127}));
  connect(pI_awfb_basic.CS, gain.u) annotation(
    Line(points = {{-294, -104}, {-254, -104}, {-254, -106}}, color = {0, 0, 127}));
  connect(sTh.oT, pI_awfb_basic.PV) annotation(
    Line(points = {{-174, -92}, {-330, -92}, {-330, -110}, {-320, -110}}, color = {0, 0, 127}));
  connect(Heater_T_Max.y, pI_awfb_basic.SP) annotation(
    Line(points = {{-362, -98}, {-320, -98}}, color = {0, 0, 127}));
  annotation(
    Diagram(coordinateSystem(extent = {{-480, 180}, {500, -100}})),
    experiment(StartTime = 0, StopTime = 864000, Tolerance = 1e-6, Interval = 86.4),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian",
    __OpenModelica_simulationFlags(lv = "LOG_STATS", s = "dassl"));
end Plant_heater_with_control_2;
