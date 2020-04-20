# Anti-windup-Compensator
The windup problem in active vibration control (AVC) is investigated in details. Instead of reviewing a list of methods on various abstract simulations, a benchmark problem in AVC is defined. Then, the proposed methods are adapted to the output regulation problem in disturbance rejection control. The selected approaches are based on their fundamental contributions to the anti-windup compensation problem. Large attention is given to capture the similarities and differences of the methods in dealing with the windup problem. Therefore, instead of categorizing the methods to static and non-static methods or model recovery and direct linear anti-windup schemes, etc., a logical route is followed to highlight the significance of each method. The mathematical interpretations of the methods are provided for the vibration engineer while delivering forthright implementation algorithms for AVC. In this regards, the methods are unified for a state space representation that is commonly used in AVC modeling based on system identification. Practical issues that may raise for each technique are mentioned in the form of some remarks, and additionally, some guidelines are provided for tuning each algorithm. Finally, in order to investigate the compensated system's performance, detailed time-domain studies are carried out by separating the transient response of the systems to three modes: linear mode, where the actuation nonlinearity is inactive. The nonlinear mode, where the windup event is in progress, and finally, the output mismatch rejection mode, where the windup incident is over but performance degradation still continues to exist.
Key-words: Active vibration control; Actuation windup; Lyapunov-based methods; Mismatch; Disturbance; Linear matrix inequality. 


Refs.

[1]	E. Pereira, J. Becedas, I. Payo, F. Ramos, V. Feliu, Control of Flexible Manipulators. Theory and Practice, in: A. Jimenez, B.M. Al Hadithi (Eds.), Robot Manip. Trends, InTech, 2010: pp. 267–296.

[2]	R. Buckingham, Advanced Studies of Flexible Robotic Manipulators, Modeling, Design, Control and Applications, Ind. Robot An Int. J. 33 (2006) ir.2006.04933aae.001. doi:10.1108/ir.2006.04933aae.001.

[3]	D. Hirano, Y. Fujii, S. Abiko, R. Lampariello, K. Nagaoka, K. Yoshida, Vibration Suppression Control of a Space Robot with Flexible Appendage based on Simple Dynamic Model*, in: 2013 IEEE/RSJ Int. Conf. Intell. Robot. Syst., IEEE, Tokyo, 2013: pp. 789–794.

[4]	M. Benosman, G. Le Vey, Control of flexible manipulators: A survey, Robotica. 22 (2004) 533–545. doi:10.1017/S0263574703005642.

[5]	S.M. Hasheminejad, A. Oveisi, Active vibration control of an arbitrary thick smart cylindrical panel with optimally placed piezoelectric sensor/actuator pairs, Int. J. Mech. Mater. Des. 12 (2016) 1–16. doi:10.1007/s10999-015-9293-2.

[6]	A. Oveisi, T. Nestorović, Robust observer-based adaptive fuzzy sliding mode controller, Mech. Syst. Signal Process. 76–77 (2016) 58–71. doi:10.1016/j.ymssp.2016.01.015.

[7]	M. Nezami, B. Gholami, Active Flutter Control of a Supersonic Honeycomb Sandwich Beam Resting on Elastic Foundation with Piezoelectric Sensor/Actuator Pair, Int. J. Struct. Stab. Dyn. 15 (2015) 1450052. doi:10.1142/S0219455414500527.

[8]	E. Omidi, S.N. Mahmoodi, W.S. Shepard, Vibration reduction in aerospace structures via an optimized modified positive velocity feedback control, Aerosp. Sci. Technol. 45 (2015) 408–415. doi:10.1016/j.ast.2015.06.012.

[9]	A. Oveisi, T. Nestorovic, Robust nonfragile observer-based H2/H∞ controller, J. Vib. Control. (2016) 1077546316651548. doi:10.1177/1077546316651548.

[10]	W. Wu, S. Jayasuriya, A QFT design methodology for feedback systems with input rate or amplitude and rate saturation, in: Am. Control Conf., 2001: pp. 376–383. doi:10.1109/ACC.2001.945575.

[11]	V.R. Marcopoli, S.M. Phillips, Analysis and synthesis tools for a class of actuator-limited multivariable control systems: A linear matrix inequality approach, Int. J. Robust Nonlinear Control. 6 (1996) 1045–1063. doi:10.1002/(SICI)1099-1239(199611)6:9/10<1045::AID-RNC268>3.0.CO;2-S.

[12]	B.G. Romanchuk, Some comments on anti-windup synthesis using LMIs, Int. J. Robust Nonlinear Control. 9 (1999) 717–734. doi:10.1002/(SICI)1099-1239(199908)9:10<717::AID-RNC430>3.0.CO;2-F.

[13]	T. Nestorović, M. Trajkov, Optimal actuator and sensor placement based on balanced reduced models, Mech. Syst. Signal Process. 36 (2013) 271–289. doi:10.1016/j.ymssp.2012.12.008.

[14]	L. Bossi, C. Rottenbacher, G. Mimmi, L. Magni, Multivariable predictive control for vibrating structures: An application, Control Eng. Pract. 19 (2011) 1087–1098. doi:10.1016/j.conengprac.2011.05.003.

[15]	A. Wills, B.M. Ninness, S. Gibson, Maximum Likelihood Estimation of state space models from frequency domain data, IEEE Trans. Automat. Contr. 54 (2009) 19–33. doi:10.1109/TAC.2008.2009485.

[16]	B. Ninness, A. Wills, A. Mills, UNIT: A freely available system identification toolbox, Control Eng. Pract. 21 (2013) 631–644. doi:10.1016/j.conengprac.2012.10.007.

[17]	B. Sevim, A. Bayraktar, A.C. Altunisik, Finite element model calibration of berke arch dam using operational modal testing, J. Vib. Control. 17 (2011) 1065–1079. doi:10.1177/1077546310377912.

[18]	D. Tcherniak, S. Chauhan, M.H. Hansen, Applicability Limits of Operational Modal Analysis to Operational Wind Turbines, in: Springer New York, 2011: pp. 317–327. doi:10.1007/978-1-4419-9716-6_29.

[19]	R. (Rik) Pintelon, J. (Johan) Schoukens, Wiley InterScience (Online service), System identification : a frequency domain approach, Wiley, 2012.

[20]	V. Kapila, K.M. Grigoriadis, Actuator saturation control, M. Dekker, 2002.

[21]	R. Hanus, M. Kinnaert, J.L. Henrotte, Conditioning technique, a general anti-windup and bumpless transfer method, Automatica. 23 (1987) 729–739. doi:10.1016/0005-1098(87)90029-X.

[22]	D. Horla, On Directional Change and Anti-Windup Compensation in Multivariable Control Systems, Int. J. Appl. Math. Comput. Sci. 19 (2009) 281–289. doi:10.2478/v10006-009-0024-4.

[23]	A. Zheng, M. V. Kothare, M. Morari, Anti-windup design for internal model control, Int. J. Control. 60 (1994) 1015–1024. doi:10.1080/00207179408921506.

[24]	A.R. Teel, N. Kapoor, The L2 anti-winup problem: Its definition and solution, in: Control Conf. (ECC), 1997 Eur., IEEE, 1997. http://ieeexplore.ieee.org/document/7082381/.

[25]	M. Saeki, N. Wada, Design of anti-windup controller based on matrix inequalities, in: Proc. 35th IEEE Conf. Decis. Control, IEEE, 1996: pp. 261–262. doi:10.1109/CDC.1996.574310.

[26]	P.J. Campo, M. Morari, C.N. Nett, Multivariable Anti-Windup and Bumpless Transfer: A General Theory, in: Am. Control Conf. 1989, 1989: pp. 1706–1711.

[27]	M. V. Kothare, P.J. Campo, M. Morari, C.N. Nett, A unified framework for the study of anti-windup designs, Automatica. 30 (1994) 1869–1883. doi:10.1016/0005-1098(94)90048-5.

[28]	N. Kapoor, A.R. Teel, P. Daoutidis, An Anti-Windup Design for Linear Systems with Input Saturation, Automatica. 34 (1998) 559–574. doi:10.1016/S0005-1098(97)00194-5.

[29]	C.A. Desoer, M. (Mathukumalli) Vidyasagar, Feedback systems: input-output properties, Academic Press, 1975.

[30]	C. Edwards, I. Postlethwaite, Anti-windup and bumpless-transfer schemes, Automatica. 34 (1998) 199–210. 
doi:10.1016/S0005-1098(97)00165-9.

[31]	J.C. Doyle, R.S. Smith, D.F. Enns, Control of Plants with Input Saturation Nonlinearities, 1987 Am. Control Conf. (1987) 1034–1039.

[32]	P.F. Weston, I. Postlethwaite, Linear conditioning for systems containing saturating actuators, Automatica. 36 (2000) 1347–1354. doi:10.1016/S0005-1098(00)00044-3.

[33]	S. Galeani, S. Onori, L. Zaccarian, Nonlinear scheduled control for linear systems subject to saturation with application to anti-windup control, in: 2007 46th IEEE Conf. Decis. Control, IEEE, 2007: pp. 1168–1173. doi:10.1109/CDC.2007.4434692.

[34]	L. Lu, Z. Lin, Design of a Nonlinear Anti-Windup Gain by Using a Composite Quadratic Lyapunov Function, IEEE Trans. Automat. Contr. 56 (2011) 2997–3001. doi:10.1109/TAC.2011.2161832.

[35]	S. Sajjadi-Kia, F. Jabbari, Modified dynamic anti-windup through deferral of activation, Int. J. Robust Nonlinear Control. 22 (2012) 1661–1673. doi:10.1002/rnc.1772.

[36]	M.C. Turner, I. Postlethwaite, A new perspective on static and low order anti-windup synthesis, Int. J. Control. 77 (2004) 27–44. doi:10.1080/00207170310001640116.

[37]	L. Zaccarian, A.R. Teel, A common framework for anti-windup, bumpless transfer and reliable designs, Automatica. 38 (2002) 1735–1744. doi:10.1016/S0005-1098(02)00072-9.

[38]	M. V. Kothare, M. Morari, Multiplier theory for stability analysis of anti-windup control systems, Automatica. 35 (1999) 917–928. doi:10.1016/S0005-1098(98)00229-5.

[39]	A.H. Glattfelder, W. Schaufelberger, Stability of discrete override and cascade-limiter single-loop control systems, IEEE Trans. Automat. Contr. 33 (1988) 532–540. doi:10.1109/9.1248.

[40]	E.F. Mulder, M. V. Kothare, M. Morari, Multivariable anti-windup controller synthesis using linear matrix inequalities, Automatica. 37 (2001) 1407–1416. doi:10.1016/S0005-1098(01)00075-9.

[41]	A. Syaichu-Rohman, R.H. Middleton, M.M. Seron, A multivariable nonlinear algebraic loop as a QP with applications to MPC, Eur. Control Conf. ECC 2003. (2003) 1–6.

[42]	E.F. Mulder, M.V. Kothare, Synthesis of stabilizing anti-windup controllers using piecewise quadratic Lyapunov functions, in: Proc. 2000 Am. Control Conf. ACC (IEEE Cat. No.00CH36334), IEEE, 2000: pp. 3239–3243 vol.5. doi:10.1109/ACC.2000.879163.

[43]	G. Grimm, J. Hatfield, I. Postlethwaite, A.R. Teel, M.C. Turner, L. Zaccarian, Antiwindup for Stable Linear Systems with Input Saturation: An LMI-Based Synthesis, IEEE Trans. Automat. Contr. 48 (2003) 1509–1525. doi:10.1109/TAC.2003.816965.

[44]	G. Li, W.P. Heath, B. Lennox, Concise stability conditions for systems with static nonlinear feedback expressed by a quadratic program, IFAC Proc. Vol. 16 (2007) 263–268. doi:10.1049/iet-cta.

[45]	A.A. Adegbege, W.P. Heath, Multivariable algebraic loops in linear anti-windup implementations, in: 2015 23rd Mediterr. Conf. Control Autom., IEEE, 2015: pp. 514–519. doi:10.1109/MED.2015.7158799.

[46]	R. Fletcher, Practical Methods of Optimization, 1987. doi:10.1097/00000539-200101000-00069.

[47]	A.R. King-Hans, W.P. Heath, R. Alli-Oke, Two-stage multivariable IMC Antiwindup (TMIA) control of a Quadruple Tank process using a PLC, in: 2014 IEEE Conf. Control Appl., IEEE, 2014: pp. 1681–1686. doi:10.1109/CCA.2014.6981554.

[48]	A.A. Adegbege, W.P. Heath, Directionality compensation for linear multivariable anti-windup synthesis, Int. J. Control. 88 (2015) 2392–2402. doi:10.1080/00207179.2015.1045556.

[49]	Y. Peng, D. Vrančić, R. Hanus, S.S.R. Weller, Anti-Windup Designs for Multivariable Controllers, Automatica. 34 (1998) 1559–1565. doi:10.1016/S0005-1098(98)80009-5.

[50]	M. Soroush, N. Mehranbod, Optimal compensation for directionality in processes with a saturating actuator, Comput. Chem. Eng. 26 (2002) 1633–1641. doi:10.1016/S0098-1354(02)00145-X.

[51]	G. Li, G. Herrmann, D.P. Stoten, J. Tu, M.C. Turner, A novel robust disturbance rejection anti-windup framework, Int. J. Control. 84 (2011) 123–137. doi:10.1080/00207179.2010.542774.

[52]	G. Grimm, I. Postlethwaite, A.R. Teel, M.C. Turner, L. Zaccarian, Linear matrix inequalities for full and reduced order anti-windup synthesis, in: Proc. Am. Control Conf., 2001: pp. 4134–4139. doi:10.1109/ACC.2001.946384.

[53]	M.C. Turner, I. Postlethwaite, G. Herrmann, Further results on full-order anti-windup synthesis: Exploiting the stability multiplier - Google Scholar, in: Proc. 6th IFAC Nonlinear Control Syst. Des. Symp., Stuttgart, 2004. 

[54]	T. Hu, A.R. Teel, L. Zaccarian, Stability and Performance for Saturated Systems via Quadratic and Nonquadratic Lyapunov Functions, IEEE Trans. Automat. Contr. 51 (2006) 1770–1786. doi:10.1109/TAC.2006.884942.

[55]	A. Syaichu-Rohman, R.H. Middleton, On the robustness of multivariable algebraic loops with sector nonlinearities, in: Proc. 41st IEEE Conf. Decis. Control. 2002., IEEE, n.d.: pp. 1054–1059. doi:10.1109/CDC.2002.1184650.

[56]	C. Scherer, P. Gahinet, M. Chilali, Multiobjective output-feedback control via LMI optimization, IEEE Trans. Automat. Contr. 42 (1997) 896–911. doi:10.1109/9.599969.

[57]	M.V. Kothare, M. Morari, Multivariable anti-windup controller synthesis using multi-objective optimization, in: Proc. 1997 Am. Control Conf. (Cat. No.97CH36041), IEEE, 1997: pp. 3093–3097 vol.5. doi:10.1109/ACC.1997.612027.

[58]	E.F. Mulder, P.Y. Tiwari, M. V. Kothare, Simultaneous linear and anti-windup controller synthesis using multiobjective convex optimization, 2009. doi:10.1016/j.automatica.2008.10.019.

[59]	S. Galeani, S. Tarbouriech, M. Turner, L. Zaccarian, A Tutorial on Modern Anti-windup Design, Eur. J. Control. 15 (2009) 418–440. doi:10.3166/ejc.15.418-440.

[60]	S. Sastry, Nonlinear Systems: Analysis, Stability, and Control, Springer Verlag, New York, 1999. 

[61]	A.A. Adegbege, W.P. Heath, Internal model control design for input-constrained multivariable processes, AIChE J. 57 (2011) 3459–3472. doi:10.1002/aic.12540.

[62]	R. Gessing, Implementability of Regulation and Partial Decoupling of MIMO Plants, in: 2006 14th Mediterr. Conf. Control Autom., IEEE, 2006: pp. 1–6. doi:10.1109/MED.2006.328742.

[63]	M. Morari, E. Zafiriou, Robust process control, Prentice Hall, Englewood Cliffs, NJ, 1989.

[64]	D.E. Kirk, Optimal control theory; an introduction, Prentice-Hall, 1970.

[65]	H. Nobahari, S.A. Hosseini Kordkheili, S.S. Afshari, Hardware-in-the-loop optimization of an active vibration controller in a flexible beam structure using evolutionary algorithms, J. Intell. Mater. Syst. Struct. 25 (2014) 1211–1223. doi:10.1177/1045389X13502874.

[66]	L. Zaccarian, A.R. Teel, Nonlinear Scheduled Anti-Windup Design for Linear Systems, IEEE Trans. Automat. Contr. 49 (2004) 2055–2061. doi:10.1109/TAC.2004.837539.

[67]	M.C. Turner, J. Sofrony, G. Herrmann, An alternative approach to anti-windup in anticipation of actuator saturation, Int. J. Robust Nonlinear Control. (2016). doi:10.1002/rnc.3609.

[68]	A. Oveisi, T. Nestorović, Mu-Synthesis based active robust vibration control of an MRI inlet, FACTA Univ. Ser. Mech. Eng. 14 (2016) 37–53.

[69]	A.R. Teel, L. Zaccarian, J.J. Marcinkowski, An anti-windup strategy for active vibration isolation systems, Control Eng. Pract. 14 (2006) 17–27. doi:10.1016/j.conengprac.2004.12.018.

[70]	Z. Qiu, J. Han, X. Zhang, Y. Wang, Z. Wu, Active vibration control of a flexible beam using a non-collocated acceleration sensor and piezoelectric patch actuator, J. Sound Vib. 326 (2009) 438–455. doi:10.1016/j.jsv.2009.05.034.

[71]	J.M. Rodriguez-Fortun, J. Orus, J. Alfonso, F.B. Gimeno, J.A. Castellanos, Flatness-Based Active Vibration Control for Piezoelectric Actuators, IEEE/ASME Trans. Mechatronics. 18 (2013) 221–229. doi:10.1109/TMECH.2011.2166998.

[72]	R. a Horn, C.R. Johnson, Matrix Analysis, 1985. doi:10.1002/zamm.19870670330.

[73]	R.W. Cottle, J.-S. Pang, R.E. Stone, The Linear Complementarity Problem, Society for Industrial and Applied Mathematics, 2009. doi:10.1137/1.9780898719000.

[74]	S.P. Dirkse, M.C. Ferris, The path solver: a nommonotone stabilization scheme for mixed complementarity problems, Optim. Methods Softw. 5 (1995) 123–156. doi:10.1080/10556789508805606.
