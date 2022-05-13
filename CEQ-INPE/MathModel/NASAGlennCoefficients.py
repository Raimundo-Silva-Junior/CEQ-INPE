
#This python file contains all Nasa Glenn Coefficients for a range of species
#considering C, H, O and N atoms.
#This data can be used to find values for termodynamics functions of ideal and real gases

#MADE BY JOSÉ RAIMUNDO DA SILVA JUNIOR (FEG - UNESP) AS A PIBIC PROJECT FOR CNPq.
#ADVISED BY FERNANDO DE SOUZA COSTA (INPE - LCP).

import cmath

class Coefficients:

    """
    The class Coefficients is a class dedicated to lists of values for NASA-Glenn Coefficients for calculating
    thermodynamic properties of individual species.

    how to read the list:
     - cl(name of the specie capitalized) means that the Coefficients are made for the temperature range of 200 < T < 1000,
     - cm(name of the specie capitalized) means that the Coefficients are made for the temperature range of 1000 < T < 6000,
     - ch(name of the specie capitalized) means that the Coefficients are made for the temperature range of 6000 < T < 20000.

    Ex: clCO2, cmH2O or chO2.
    Ex: clN2 = (a1, a2, a3, a4, a5, a6, a7, b1, b2)
    """
    comb_list = (
    'CO2', 'CO', 'H2O', 'H2', 'N2', 'H', 'N', 'CH4', 'C2H2(ACETY)', 'C2H2(VINY)', 
    'C2H4', 'C2H6', 'C3H8', 'C4H8(ISOBUT)', 'C4H10(NBUT)', 'C4H10(ISOBUT)', 'C5H12(NPENT)',
    'C6H14(NHEX)', 'C7H16(NHEP)', 'C7H16(2METH)', 'C8H18(NOCT)', 'C8H18(ISOCT)', 'NH', 'NH2',
    'NH3', 'N2H2', 'N2H4', 'N3H', 'CH3NO2(L)', 'CH6N2(L)', 'C2H8N2(L)', 'CH4(L)'
    )

    #TUPLE OF OXIDANT SPECIES CHOSEN FOR CHON + HON PROBLEM IN ORDER

    oxid_list = (
    'N2', 'N', 'H2O', 'O', 'O2', 'NO', 'NO2', 'N2O', 'N2O3', 'N2O4', 'N2O5', 
    'OH', 'HO2', 'H2O2', 'O2(L)'
    )
    
    clCO2 = (4.943650540e+04, -6.264116010e+02, 5.301725240e+00, 2.503813816e-03, -2.127308728e-07, -7.689988780e-10, 2.849677801e-13,-4.528198460e+04,-7.048279440e+00)
    cmCO2 = (1.176962419e+05, -1.788791477e+03, 8.291523190e+00,-9.223156780e-05, 4.863676880e-09, -1.891053312e-12, 6.330036590e-16,-3.908350590e+04,-2.652669281e+01)
    chCO2 = (-1.544423287e+09, 1.016847056e+06, -2.561405230e+02, 3.369401080e-02, -2.181184337e-06, 6.991420840e-11, -8.842351500e-16, -8.043214510e+06, 2.254177493e+03)
    clCO = (1.489045326e+04, -2.922285939e+02, 5.724527170e+00,-8.176235030e-03, 1.456903469e-05,-1.087746302e-08, 3.027941827e-12,-1.303131878e+04,-7.859241350e+00)
    cmCO = (4.619197250e+05, -1.944704863e+03, 5.916714180e+00,-5.664282830e-04, 1.398814540e-07,-1.787680361e-11, 9.620935570e-16,-2.466261084e+03,-1.387413108e+01)
    chCO = (8.868662960e+08, -7.500377840e+05, 2.495474979e+02, -3.956351100e-02, 3.297772080e-06, -1.318409933e-10, 1.998937948e-15, 5.701421130e+06,-2.060704786e+03)
    clH2O = (-3.947960830e+04, 5.755731020e+02, 9.317826530e-01, 7.222712860e-03,-7.342557370e-06, 4.955043490e-09, -1.336933246e-12, -3.303974310e+04, 1.724205775e+01)
    cmH2O = (1.034972096e+06, -2.412698562e+03, 4.646110780e+00, 2.291998307e-03,-6.836830480e-07, 9.426468930e-11, -4.822380530e-15, -1.384286509e+04, -7.978148510e+00)
    clH2 = (4.078323210e+04, -8.009186040e+02, 8.214702010e+00, -1.269714457e-02, 1.753605076e-05, -1.202860270e-08,  3.368093490e-12, 2.682484665e+03, -3.043788844e+01)
    cmH2 = (5.608128010e+05, -8.371504740e+02, 2.975364532e+00, 1.252249124e-03, -3.740716190e-07, 5.936625200e-11, -3.606994100e-15, 5.339824410e+03, -2.202774769e+00)
    chH2 = (4.966884120e+08, -3.147547149e+05, 7.984121880e+01, -8.414789210e-03, 4.753248350e-07, -1.371873492e-11, 1.605461756e-16, 2.488433516e+06, -6.695728110e+02)
    clN2 = (2.210371497e+04, -3.818461820e+02, 6.082738360e+00, -8.530914410e-03, 1.384646189e-05, -9.625793620e-09, 2.519705809e-12, 7.108460860e+02, -1.076003744e+01)
    cmN2 = (5.877124060e+05, -2.239249073e+03, 6.066949220e+00, -6.139685500e-04, 1.491806679e-07, -1.923105485e-11, 1.061954386e-15, 1.283210415e+04, -1.586640027e+01)
    chN2 = (8.310139160e+08, -6.420733540e+05, 2.020264635e+02, -3.065092046e-02, 2.486903333e-06, -9.705954110e-11, 1.437538881e-15, 4.938707040e+06, -1.672099740e+03)
    clAr = (0.000000000e+00, 0.000000000e+00, 2.500000000e+00, 0.000000000e+00, 0.000000000e+00, 0.000000000e+00, 0.000000000e+00, -7.453750000e+02, 4.379674910e+00)
    cmAr = (2.010538475e+01, -5.992661070e-02, 2.500069401e+00, -3.992141160e-08, 1.205272140e-11, -1.819015576e-15, 1.078576636e-19, -7.449939610e+02, 4.379180110e+00)
    chAr = (-9.951265080e+08, 6.458887260e+05, -1.675894697e+02, 2.319933363e-02, -1.721080911e-06, 6.531938460e-11, -9.740147729e-16, -5.078300340e+06, 1.465298484e+03)
    clO2 = (-3.425563420e+04, 4.847000970e+02, 1.119010961e+00, 4.293889240e-03, -6.836300520e-07,-2.023372700e-09, 1.039040018e-12, -3.391454870e+03, 1.849699470e+01)
    cmO2 = (-1.037939022e+06, 2.344830282e+03, 1.819732036e+00, 1.267847582e-03, -2.188067988e-07, 2.053719572e-11, -8.193467050e-16, -1.689010929e+04, 1.738716506e+01)
    chO2 = (4.975294300e+08, -2.866106874e+05, 6.690352250e+01, -6.169959020e-03, 3.016396027e-07, -7.421416600e-12, 7.278175770e-17, 2.293554027e+06, -5.530621610e+02)
    clH1 = (0.000000000e+00, 0.000000000e+00, 2.500000000e+00, 0.000000000e+00, 0.000000000e+00, 0.000000000e+00, 0.000000000e+00, 2.547370801e+04, -4.466828530e-01)
    cmH1 = (6.078774250e+01,-1.819354417e-01, 2.500211817e+00, -1.226512864e-07, 3.732876330e-11, -5.687744560e-15, 3.410210197e-19, 2.547486398e+04, -4.481917770e-01)
    chH1 = (2.173757694e+08, -1.312035403e+05, 3.399174200e+01, -3.813999680e-03, 2.432854837e-07, -7.694275540e-12, 9.644105630e-17, 1.067638086e+06, -2.742301051e+02)
    clN1 = (0.000000000e+00, 0.000000000e+00, 2.500000000e+00, 0.000000000e+00, 0.000000000e+00, 0.000000000e+00, 0.000000000e+00, 5.610463780e+04, 4.193905036e+00)
    cmN1 = (8.876501380e+04, -1.071231500e+02, 2.362188287e+00, 2.916720081e-04, -1.729515100e-07, 4.012657880e-11, -2.677227571e-15, 5.697351330e+04, 4.865231506e+00)
    chN1 = (5.475181050e+08, -3.107574980e+05, 6.916782740e+01, -6.847988130e-03, 3.827572400e-07, -1.098367709e-11, 1.277986024e-16, 2.550585618e+06, -5.848769753e+02)
    clO1 = (-7.953611300e+03, 1.607177787e+02, 1.966226438e+00, 1.013670310e-03, -1.110415423e-06, 6.517507500e-10, -1.584779251e-13, 2.840362437e+04, 8.404241820e+00)
    cmO1 = (2.619020262e+05, -7.298722030e+02, 3.317177270e+00, -4.281334360e-04, 1.036104594e-07, 9.438304330e-12, -2.725038297e-16, 3.392428060e+04, -6.679585350e-01)
    chO1 = (1.779004264e+08, -1.082328257e+05, 2.810778365e+01, -2.975232262e-03, 1.854997534e-07, -5.796231540e-12, 7.191720164e-17, 8.890942630e+05, -2.181728151e+02)
    clNO = (-1.143916503e+04, 1.536467592e+02, 3.431468730e+00, -2.668592368e-03, 8.481399120e-06, -7.685111050e-09, 2.386797655e-12, 9.098214410e+03, 6.728725490e+00)
    cmNO = (2.239018716e+05, -1.289651623e+03, 5.433936030e+00, -3.656034900e-04, 9.880966450e-08, -1.416076856e-11, 9.380184620e-16, 1.750317656e+04, -8.501669090e+00)
    chNO = (-9.575303540e+08, 5.912434480e+05, -1.384566826e+02, 1.694339403e-02, -1.007351096e-06, 2.912584076e-11, -3.295109350e-16, -4.677501240e+06, 1.242081216e+03)
    clOH = (-1.998858990e+03, 9.300136160e+01, 3.050854229e+00, 1.529529288e-03, -3.157890998e-06, 3.315446180e-09, -1.138762683e-12, 2.991214235e+03, 4.674110790e+00)
    cmOH = (1.017393379e+06, -2.509957276e+03, 5.116547860e+00, 1.305299930e-04, -8.284322260e-08, 2.006475941e-11, -1.556993656e-15, 2.019640206e+04, -1.101282337e+01)
    chOH = (2.847234193e+08, -1.859532612e+05, 5.008240900e+01, -5.142374980e-03, 2.875536589e-07, -8.228817960e-12, 9.567229020e-17, 1.468393908e+06, -4.023555580e+02)
    clAIR = (1.009950160e+04, -1.968275610e+02, 5.009155110e+00, -5.761013730e-03, 1.066859930e-05, -7.940297970e-09, 2.185231910e-12, -1.767967310e+02, -3.921504225e+00)
    cmAIR = (2.415214430e+05, -1.257874600e+03, 5.144558670e+00, -2.138541790e-04, 7.065227840e-08, -1.071483490e-11, 6.577800150e-16, 6.462263190e+03, -8.147411905e+00)
    clNO2 = (-5.642038780e+04, 9.633085720e+02, -2.434510974e+00, 1.927760886e-02, -1.874559328e-05, 9.145497730e-09, -1.777647635e-12, -1.547925037e+03, 4.067851210e+01)
    cmNO2 = (7.213001570e+05, -3.832615200e+03, 1.113963285e+01, -2.238062246e-03, 6.547723430e-07, -7.611335900e-11, 3.328361050e-15, 2.502497403e+04, -4.305130040e+01)
    clN2O = (4.288225970e+04, -6.440118440e+02, 6.034351430e+00, 2.265394436e-04, 3.472782850e-06, -3.627748640e-09, 1.137969552e-12, 1.179405506e+04, -1.003128570e+01)
    cmN2O = (3.438448040e+05, -2.404557558e+03, 9.125636220e+00, -5.401667930e-04, 1.315124031e-07, -1.414215100e-11, 6.381066870e-16,  2.198632638e+04, -3.147805016e+01)
    clN2O3 =(-9.204444170e+04, 9.295520150e+02, 3.203664810e+00, 1.356473078e-02, -6.262966070e-06, -1.402915559e-09, 1.431620930e-12, 3.313622080e+03, 1.844430953e+01)
    cmN2O3 = (7.783881860e+05, -4.483024660e+03, 1.666668024e+01, -2.062143878e-03, 5.309541710e-07, -6.190451220e-11, 2.692956658e-15, 3.360912450e+04, -6.739212388e+01)
    clN2O4 = (-3.804751440e+04, 5.612828890e+02, -2.083648324e-01, 3.887087820e-02, -4.422412260e-05, 2.498812310e-08, -5.679102380e-12, -3.310794730e+03, 2.963924840e+01)
    cmN2O4 = (-4.582843760e+05, -1.604749805e+03, 1.674102133e+01, -5.091385080e-04, 1.143634670e-07, -1.316288176e-11, 5.976316620e-16, 4.306900520e+03,-6.569450380e+01)
    clN2O5 = (4.007828170e+04, -8.769675120e+02, 1.055932981e+01, 1.394613859e-02, -8.884346920e-06, 8.500431150e-10, 7.791550910e-13, 3.038962037e+03, -2.386831860e+01)
    cmN2O5 = (-5.325578960e+04, -3.109277389e+03, 2.036088958e+01, -9.959901140e-04, 2.401398635e-07, -3.057161911e-11, 1.495915511e-15, 1.336957281e+04, -8.298623341e+01)
    clHO2 = (-7.598882540e+04, 1.329383918e+03, -4.677388240e+00, 2.508308202e-02, -3.006551588e-05, 1.895600056e-08, -4.828567390e-12, -5.873350960e+03, 5.193602140e+01)
    cmHO2 = (-1.810669724e+06, 4.963192030e+03, -1.039498992e+00, 4.560148530e-03, -1.061859447e-06, 1.144567878e-10, -4.763064160e-15, -3.200817190e+04, 4.066850920e+01)
    clH2O2 = (-9.279533580e+04, 1.564748385e+03, -5.976460140e+00, 3.270744520e-02, -3.932193260e-05, 2.509255235e-08, -6.465045290e-12, -2.494004728e+04, 5.877174180e+01)
    cmH2O2 = (1.489428027e+06, -5.170821780e+03, 1.128204970e+01, -8.042397790e-05, -1.818383769e-08, 6.947265590e-12, -4.827831900e-16, 1.418251038e+04, -4.650855660e+01)
    clCH4 = (-1.766850998e+05, 2.786181020e+03, -1.202577850e+01, 3.917619290e-02, -3.619054430e-05, 2.026853043e-08, -4.976705490e-12, -2.331314360e+04, 8.904322750e+01)
    cmCH4 = (3.730042760e+06, -1.383501485e+04, 2.049107091e+01, -1.961974759e-03,  4.727313040e-07, -3.728814690e-11, 1.623737207e-15, 7.532066910e+04, -1.219124889e+02)
    clC2H2_acety = (1.598112089e+05, -2.216644118e+03, 1.265707813e+01, -7.979651080e-03, 8.054992750e-06, -2.433307673e-09, -7.529233180e-14, 3.712619060e+04, -5.244338900e+01)
    cmC2H2_acety = (1.713847410e+06, -5.929106660e+03, 1.236127943e+01, 1.314186993e-04, -1.362764431e-07, 2.712655786e-11, -1.302066204e-15, 6.266578970e+04, -5.818960590e+01)
    clC2H2_viny = (-1.466042239e+04, 2.789475593e+02, 1.276229776e+00, 1.395015463e-02, -1.475702649e-05, 9.476298110e-09, -2.567602217e-12, 4.736110180e+04, 1.658225704e+01)
    cmC2H2_viny = (1.940838725e+06, -6.892718150e+03, 1.339582494e+01, -9.368968670e-04, 1.470804368e-07, -1.220040365e-11, 4.122391660e-16, 9.107112930e+04, -6.337502930e+01)
    clC2H4 = (-1.163605836e+05, 2.554851510e+03, -1.609746428e+01, 6.625779320e-02, -7.885081860e-05, 5.125224820e-08, -1.370340031e-11, -6.176191070e+03, 1.093338343e+02)
    cmC2H4 = (3.408763670e+06, -1.374847903e+04, 2.365898074e+01, -2.423804419e-03, 4.431395660e-07, -4.352683390e-11, 1.775410633e-15, 8.820429380e+04, -1.371278108e+02)
    clC2H6 = (-1.862044161e+05, 3.406191860e+03, -1.951705092e+01, 7.565835590e-02, -8.204173220e-05, 5.061135800e-08, -1.319281992e-11, -2.702932890e+04, 1.298140496e+02)
    cmC2H6 = (5.025782130e+06, -2.033022397e+04, 3.322552930e+01, -3.836703410e-03, 7.238405860e-07, -7.319182500e-11, 3.065468699e-15, 1.115963950e+05, -2.039410584e+02)
    clC3H8 = (-2.433144337e+05, 4.656270810e+03, -2.939466091e+01, 1.188952745e-01, -1.376308269e-04, 8.814823910e-08, -2.342987994e-11, -3.540335270e+04, 1.841749277e+02)
    cmC3H8 = (6.420731680e+06, -2.659791134e+04, 4.534356840e+01, -5.020663920e-03, 9.471216940e-07, -9.575405230e-11, 4.009672880e-15, 1.455582459e+05, -2.818374734e+02)
    clC4H8_isobut = (-2.327205032e+05, 3.941994240e+03, -2.224581184e+01, 1.012790864e-01, -1.073194065e-04, 6.454696910e-08, -1.646330345e-11, -2.233766063e+04, 1.479597621e+02)
    cmC4H8_isobut = (6.484970990e+06, -2.732504764e+04, 4.836321080e+01, -4.768004050e-03, 8.233875840e-07, -7.449253000e-11, 2.782303056e-15, 1.595941773e+05, -2.982986237e+02)
    clC4H10_nbut = (-3.175872540e+05, 6.176331820e+03, -3.891562120e+01, 1.584654284e-01, -1.860050159e-04, 1.199676349e-07, -3.201670550e-11, -4.540363390e+04, 2.379488665e+02)
    cmC4H10_nbut = (7.682322450e+06, -3.256051510e+04, 5.736732750e+01, -6.197916810e-03, 1.180186048e-06, -1.221893698e-10, 5.250635250e-15, 1.774526560e+05, -3.587918760e+02)
    clC4H10_isobut = (-3.834469330e+05, 7.000039640e+03, -4.440026900e+01, 1.746183447e-01, -2.078195348e-04, 1.339792433e-07, -3.551681630e-11, -5.034018890e+04, 2.658966497e+02)
    cmC4H10_isobut = (7.528018920e+06, -3.202517060e+04, 5.700161000e+01, -6.060013090e-03, 1.143975809e-06, -1.157061835e-10, 4.846042910e-15, 1.728500802e+05, -3.576176890e+02)
    clC5H12_npent = (-2.768894625e+05, 5.834283470e+03, -3.617541480e+01, 1.533339707e-01, -1.528395882e-04, 8.191092000e-08, -1.792327902e-11,  -4.665375250e+04, 2.265544053e+02)
    cmC5H12_npent = (-2.530779286e+06, -8.972593260e+03, 4.536223260e+01, -2.626989916e-03, 3.135136419e-06, -5.318728940e-10, 2.886896868e-14,  1.484616529e+04, -2.516550384e+02)
    clC6H14_nhex = (-5.815926700e+05, 1.079097724e+04, -6.633947030e+01, 2.523715155e-01, -2.904344705e-04, 1.802201514e-07, -4.617223680e-11,  -7.271544570e+04, 3.938283540e+02)
    cmC6H14_nhex = (-3.106625684e+06, -7.346087920e+03, 4.694131760e+01, 1.693963977e-03, 2.068996667e-06, -4.212141680e-10, 2.452345845e-14, 5.237503120e+02, -2.549967718e+02)
    clC7H16_nhep = (-6.127432890e+05, 1.184085437e+04, -7.487188600e+01, 2.918466052e-01, -3.416795490e-04, 2.159285269e-07, -5.655852730e-11, -8.013408940e+04, 4.407213320e+02)
    cmC7H16_nhep = (9.135632470e+06, -3.923319690e+04, 7.889780850e+01, -4.654251930e-03, 2.071774142e-06, -3.442539300e-10, 1.976834775e-14,  2.050708295e+05, -4.851104020e+02)
    clC7H16_2meth = (-7.104777770e+05, 1.191251120e+04, -7.345339440e+01, 2.902952369e-01, -3.462767680e-04, 2.260184498e-07, -6.128813920e-11, -8.202147700e+04, 4.320042290e+02)
    cmC7H16_2meth = (1.289912969e+06, -1.784340963e+03, 1.083537673e+01, 5.270609240e-02, -1.886832314e-05, 2.432255843e-09, -1.135553789e-13, -1.637529884e+04, -2.981862410e+01)
    clC8H18_noct = (-6.986647150e+05, 1.338501096e+04, -8.415165920e01, 3.271936660e-01, -3.777209590e-04, 2.339836988e-07, -6.010892650e-11, -9.026223250e+04, 4.939222140e+02)
    cmC8H18_noct = (6.365406950e+06, -3.105364657e+04, 6.969162340e+01, 1.048059637e-02, -4.129621950e-06, 5.543226320e-10, -2.651436499e-14, 1.500968785e+05,-4.169895650e+02)
    clC8H18_isoct = (-1.688758565e+05, 3.126903227e+03, -2.123502828e+01, 1.489151508e-01, -1.151180135e-04, 4.473216170e-08, -5.554882070e-12, -4.468060620e+04, 1.417455793e+02)
    cmC8H18_isoct = (1.352765032e+07,-4.663370340e+04, 7.795313180e+01, 1.423729984e-02, -5.073593910e-06, 7.248232970e-10, -3.819190110e-14, 2.541178017e+05, -4.933887190e+02)
    clNH = (1.359651320e+04, -1.900296604e+02, 4.518496790e+00, -2.432776899e-03, 2.377587464e-06, -2.592797084e-10, -2.659680792e-13, 4.280972190e+04, -3.886561616e+00)
    cmNH = (1.958141991e+06, -5.782861300e+03, 9.335742020e+00, -2.292910311e-03, 6.076092480e-07, -6.647942750e-11, 2.384234783e-15, 7.898912340e+04, -4.116970400e+01)
    chNH = (9.524636790e+07,-8.585826910e+04, 2.980445181e+01, -2.979563697e-03, 1.656334158e-07, -4.744791840e-12, 5.570148290e-17, 6.961434270e+05, -2.229027419e+02)
    clNH2 = (-3.118240659e+04, 4.754243390e+02, 1.372395176e+00, 6.306429720e-03, -5.987893560e-06, 4.492752340e-09, -1.414073548e-12, 1.928939662e+04, 1.540126885e+01)
    cmNH2 = (2.111053740e+06,-6.880627230e+03, 1.132305924e+01, -1.829236741e-03, 5.643890090e-07, -7.886452480e-11, 4.078593450e-15, 6.503778560e+04,-5.359155744e+01)
    clNH3 = (-7.681226150e+04, 1.270951578e+03, -3.893229130e+00, 2.145988418e-02, -2.183766703e-05, 1.317385706e-08, -3.332322060e-12, -1.264886413e+04, 4.366014588e+01)
    cmNH3 = (2.452389535e+06,-8.040894240e+03, 1.271346201e+01, -3.980186580e-04, 3.552502750e-08, 2.530923570e-12, -3.322700530e-16, 4.386191960e+04,-6.462330602e+01)
    clN2H2 = (-1.504005163e+05, 2.346687716e+03, -9.405430290e+00, 3.284299800e-02, -3.121920401e-05, 1.721283190e-08, -4.014537220e-12, 1.319384041e+04, 7.832382630e+01)
    cmN2H2 = (6.217567870e+06, -1.753952096e+04, 2.022730509e+01, -9.757297660e-04, -4.208416740e-07, 1.117921171e-10, -7.627102210e-15, 1.374152574e+05, -1.199559168e+02)
    clN2H4 = (-1.660756354e+05, 3.035416736e+03, -1.736889823e+01, 7.159834020e-02, -8.866799300e-05, 5.798970280e-08, -1.530037218e-11, -3.731927230e+03, 1.190002218e+02)
    cmN2H4 = (3.293486700e+06,-1.199850628e+04, 2.104406814e+01, -1.399381724e-03, 1.933173351e-07, -1.318016127e-11, 3.166400170e-16, 8.348433700e+04, -1.155751024e+02)
    clN3H = (3.242576060e+03, 6.692664890e+01, 1.766142217e+00, 1.487411419e-02, -1.539086440e-05, 9.172303550e-09, -2.337205474e-12, 3.392069700e+04, 1.513752057e+01)
    cmN3H = (1.170469241e+06,-5.102451990e+03, 1.278288910e+01, -8.409487160e-04, 1.592142834e-07, -1.512289051e-11, 6.102906630e-16, 6.428344470e+04, -5.513119108e+01)

class Termodata:
    """
    The Termodata class is dedicated to assign NASA-Glenn Coefficients for each individual specie, and then,
    use it to calculate the specie heat capacity, enthalpy and entropy for ideal gases.

  
    """

    def __init__(self, T: float, specie:str):
        """
        The __init__ method, from Termodata class, is dedicated to assign NASA-Glenn Coefficients for each individual specie, and then,
        use it to calculate the specie heat capacity, enthalpy and entropy for ideal gases.

        Args:
            T (float): Temprature from 200 < T < 6000
            specie (str): Name of the specie, ex: CO2, CO, H2O, N2...
        """

        global a1, a2, a3, a4, a5, a6, a7, b1, b2
        self.T = T
        R = 8314.46261815324
        specie = specie.upper()

        if specie == 'CO2':
            if 200 <= T < 1000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.clCO2[0], Coefficients.clCO2[1], Coefficients.clCO2[
                    2], \
                                                     Coefficients.clCO2[3], Coefficients.clCO2[4], Coefficients.clCO2[
                                                         5], \
                                                     Coefficients.clCO2[6], Coefficients.clCO2[7], Coefficients.clCO2[8]
            elif 1000 <= T < 6000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.cmCO2[0], Coefficients.cmCO2[1], Coefficients.cmCO2[
                    2], \
                                                     Coefficients.cmCO2[3], Coefficients.cmCO2[4], Coefficients.cmCO2[
                                                         5], \
                                                     Coefficients.cmCO2[6], Coefficients.cmCO2[7], Coefficients.cmCO2[8]
            elif 6000 <= T <= 20000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.chCO2[0], Coefficients.chCO2[1], Coefficients.chCO2[
                    2], \
                                                     Coefficients.chCO2[3], Coefficients.chCO2[4], Coefficients.chCO2[
                                                         5], \
                                                     Coefficients.chCO2[6], Coefficients.chCO2[7], Coefficients.chCO2[8]
        elif specie == 'CO':
            if 200 <= T < 1000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.clCO[0], Coefficients.clCO[1], Coefficients.clCO[2], \
                                                     Coefficients.clCO[3], Coefficients.clCO[4], Coefficients.clCO[5], \
                                                     Coefficients.clCO[6], Coefficients.clCO[7], Coefficients.clCO[8]
            elif 1000 <= T < 6000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.cmCO[0], Coefficients.cmCO[1], Coefficients.cmCO[2], \
                                                     Coefficients.cmCO[3], Coefficients.cmCO[4], Coefficients.cmCO[5], \
                                                     Coefficients.cmCO[6], Coefficients.cmCO[7], Coefficients.cmCO[8]
            elif 6000 <= T <= 20000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.chCO[0], Coefficients.chCO[1], Coefficients.chCO[2], \
                                                     Coefficients.chCO[3], Coefficients.chCO[4], Coefficients.chCO[5], \
                                                     Coefficients.chCO[6], Coefficients.chCO[7], Coefficients.chCO[8]
        elif specie == 'H2O':
            if 200 <= T < 1000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.clH2O[0], Coefficients.clH2O[1], Coefficients.clH2O[
                    2], \
                                                     Coefficients.clH2O[3], Coefficients.clH2O[4], Coefficients.clH2O[
                                                         5], \
                                                     Coefficients.clH2O[6], Coefficients.clH2O[7], Coefficients.clH2O[8]
            elif 1000 <= T <= 6000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.cmH2O[0], Coefficients.cmH2O[1], Coefficients.cmH2O[
                    2], \
                                                     Coefficients.cmH2O[3], Coefficients.cmH2O[4], Coefficients.cmH2O[
                                                         5], \
                                                     Coefficients.cmH2O[6], Coefficients.cmH2O[7], Coefficients.cmH2O[8]
        elif specie == 'H2':
            if 200 <= T < 1000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.clH2[0], Coefficients.clH2[1], Coefficients.clH2[2], \
                                                     Coefficients.clH2[3], Coefficients.clH2[4], Coefficients.clH2[5], \
                                                     Coefficients.clH2[6], Coefficients.clH2[7], Coefficients.clH2[8]
            elif 1000 <= T < 6000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.cmH2[0], Coefficients.cmH2[1], Coefficients.cmH2[2], \
                                                     Coefficients.cmH2[3], Coefficients.cmH2[4], Coefficients.cmH2[5], \
                                                     Coefficients.cmH2[6], Coefficients.cmH2[7], Coefficients.cmH2[8]
            elif 6000 <= T <= 20000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.chH2[0], Coefficients.chH2[1], Coefficients.chH2[2], \
                                                     Coefficients.chH2[3], Coefficients.chH2[4], Coefficients.chH2[5], \
                                                     Coefficients.chH2[6], Coefficients.chH2[7], Coefficients.chH2[8]
        elif specie == 'N2':
            if 200 <= T < 1000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.clN2[0], Coefficients.clN2[1], Coefficients.clN2[2], \
                                                     Coefficients.clN2[3], Coefficients.clN2[4], Coefficients.clN2[5], \
                                                     Coefficients.clN2[6], Coefficients.clN2[7], Coefficients.clN2[8]
            elif 1000 <= T < 6000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.cmN2[0], Coefficients.cmN2[1], Coefficients.cmN2[2], \
                                                     Coefficients.cmN2[3], Coefficients.cmN2[4], Coefficients.cmN2[5], \
                                                     Coefficients.cmN2[6], Coefficients.cmN2[7], Coefficients.cmN2[8]
            elif 6000 <= T <= 20000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.chN2[0], Coefficients.chN2[1], Coefficients.chN2[2], \
                                                     Coefficients.chN2[3], Coefficients.chN2[4], Coefficients.chN2[5], \
                                                     Coefficients.chN2[6], Coefficients.chN2[7], Coefficients.chN2[8]
        elif specie == 'Ar':
            if 200 <= T < 1000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.clAr[0], Coefficients.clAr[1], Coefficients.clAr[2], \
                                                     Coefficients.clAr[3], Coefficients.clAr[4], Coefficients.clAr[5], \
                                                     Coefficients.clAr[6], Coefficients.clAr[7], Coefficients.clAr[8]
            elif 1000 <= T < 6000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.cmAr[0], Coefficients.cmAr[1], Coefficients.cmAr[2], \
                                                     Coefficients.cmAr[3], Coefficients.cmAr[4], Coefficients.cmAr[5], \
                                                     Coefficients.cmAr[6], Coefficients.cmAr[7], Coefficients.cmAr[8]
            elif 6000 <= T <= 20000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.chAr[0], Coefficients.chAr[1], Coefficients.chAr[2], \
                                                     Coefficients.chAr[3], Coefficients.chAr[4], Coefficients.chAr[5], \
                                                     Coefficients.chAr[6], Coefficients.chAr[7], Coefficients.chAr[8]
        elif specie == 'O2':
            if 200 <= T < 1000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.clO2[0], Coefficients.clO2[1], Coefficients.clO2[2], \
                                                     Coefficients.clO2[3], Coefficients.clO2[4], Coefficients.clO2[5], \
                                                     Coefficients.clO2[6], Coefficients.clO2[7], Coefficients.clO2[8]
            elif 1000 <= T < 6000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.cmO2[0], Coefficients.cmO2[1], Coefficients.cmO2[2], \
                                                     Coefficients.cmO2[3], Coefficients.cmO2[4], Coefficients.cmO2[5], \
                                                     Coefficients.cmO2[6], Coefficients.cmO2[7], Coefficients.cmO2[8]
            elif 6000 <= T <= 20000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.chO2[0], Coefficients.chO2[1], Coefficients.chO2[2], \
                                                     Coefficients.chO2[3], Coefficients.chO2[4], Coefficients.chO2[5], \
                                                     Coefficients.chO2[6], Coefficients.chO2[7], Coefficients.chO2[8]
        elif specie == 'H':
            if 200 <= T < 1000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.clH1[0], Coefficients.clH1[1], Coefficients.clH1[2], \
                                                     Coefficients.clH1[3], Coefficients.clH1[4], Coefficients.clH1[5], \
                                                     Coefficients.clH1[6], Coefficients.clH1[7], Coefficients.clH1[8]
            elif 1000 <= T < 6000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.cmH1[0], Coefficients.cmH1[1], Coefficients.cmH1[2], \
                                                     Coefficients.cmH1[3], Coefficients.cmH1[4], Coefficients.cmH1[5], \
                                                     Coefficients.cmH1[6], Coefficients.cmH1[7], Coefficients.cmH1[8]
            elif 6000 <= T <= 20000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.chH1[0], Coefficients.chH1[1], Coefficients.chH1[2], \
                                                     Coefficients.chH1[3], Coefficients.chH1[4], Coefficients.chH1[5], \
                                                     Coefficients.chH1[6], Coefficients.chH1[7], Coefficients.chH1[8]
        elif specie == 'N':
            if 200 <= T < 1000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.clN1[0], Coefficients.clN1[1], Coefficients.clN1[2], \
                                                     Coefficients.clN1[3], Coefficients.clN1[4], Coefficients.clN1[5], \
                                                     Coefficients.clN1[6], Coefficients.clN1[7], Coefficients.clN1[8]
            elif 1000 <= T < 6000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.cmN1[0], Coefficients.cmN1[1], Coefficients.cmN1[2], \
                                                     Coefficients.cmN1[3], Coefficients.cmN1[4], Coefficients.cmN1[5], \
                                                     Coefficients.cmN1[6], Coefficients.cmN1[7], Coefficients.cmN1[8]
            elif 6000 <= T <= 20000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.chN1[0], Coefficients.chN1[1], Coefficients.chN1[2], \
                                                     Coefficients.chN1[3], Coefficients.chN1[4], Coefficients.chN1[5], \
                                                     Coefficients.chN1[6], Coefficients.chN1[7], Coefficients.chN1[8]
        elif specie == 'O':
            if 200 <= T < 1000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.clO1[0], Coefficients.clO1[1], Coefficients.clO1[2], \
                                                     Coefficients.clO1[3], Coefficients.clO1[4], Coefficients.clO1[5], \
                                                     Coefficients.clO1[6], Coefficients.clO1[7], Coefficients.clO1[8]
            elif 1000 <= T < 6000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.cmO1[0], Coefficients.cmO1[1], Coefficients.cmO1[2], \
                                                     Coefficients.cmO1[3], Coefficients.cmO1[4], Coefficients.cmO1[5], \
                                                     Coefficients.cmO1[6], Coefficients.cmO1[7], Coefficients.cmO1[8]
            elif 6000 <= T <= 20000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.chO1[0], Coefficients.chO1[1], Coefficients.chO1[2], \
                                                     Coefficients.chO1[3], Coefficients.chO1[4], Coefficients.chO1[5], \
                                                     Coefficients.chO1[6], Coefficients.chO1[7], Coefficients.chO1[8]
        elif specie == 'NO':
            if 200 <= T < 1000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.clNO[0], Coefficients.clNO[1], Coefficients.clNO[2], \
                                                     Coefficients.clNO[3], Coefficients.clNO[4], Coefficients.clNO[5], \
                                                     Coefficients.clNO[6], Coefficients.clNO[7], Coefficients.clNO[8]
            elif 1000 <= T < 6000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.cmNO[0], Coefficients.cmNO[1], Coefficients.cmNO[2], \
                                                     Coefficients.cmNO[3], Coefficients.cmNO[4], Coefficients.cmNO[5], \
                                                     Coefficients.cmNO[6], Coefficients.cmNO[7], Coefficients.cmNO[8]
            elif 6000 <= T <= 20000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.chNO[0], Coefficients.chNO[1], Coefficients.chNO[2], \
                                                     Coefficients.chNO[3], Coefficients.chNO[4], Coefficients.chNO[5], \
                                                     Coefficients.chNO[6], Coefficients.chNO[7], Coefficients.chNO[8]
        elif specie == 'OH':
            if 200 <= T < 1000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.clOH[0], Coefficients.clOH[1], Coefficients.clOH[2], \
                                                     Coefficients.clOH[3], Coefficients.clOH[4], Coefficients.clOH[5], \
                                                     Coefficients.clOH[6], Coefficients.clOH[7], Coefficients.clOH[8]
            elif 1000 <= T < 6000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.cmOH[0], Coefficients.cmOH[1], Coefficients.cmOH[2], \
                                                     Coefficients.cmOH[3], Coefficients.cmOH[4], Coefficients.cmOH[5], \
                                                     Coefficients.cmOH[6], Coefficients.cmOH[7], Coefficients.cmOH[8]
            elif 6000 <= T <= 20000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.chOH[0], Coefficients.chOH[1], Coefficients.chOH[2], \
                                                     Coefficients.chOH[3], Coefficients.chOH[4], Coefficients.chOH[5], \
                                                     Coefficients.chOH[6], Coefficients.chOH[7], Coefficients.chOH[8]
        elif specie == 'AIR':
            if 200 <= T < 1000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.clAIR[0], Coefficients.clAIR[1],  Coefficients.clAIR[2], \
                                                 Coefficients.clAIR[3], Coefficients.clAIR[4], Coefficients.clAIR[5], \
                                                 Coefficients.clAIR[6], Coefficients.clAIR[7], Coefficients.clAIR[8]
            elif 1000 <= T < 6000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.cmAIR[0], Coefficients.cmAIR[1],  Coefficients.cmAIR[2], \
                                                 Coefficients.cmAIR[3], Coefficients.cmAIR[4], Coefficients.cmAIR[5], \
                                                 Coefficients.cmAIR[6], Coefficients.cmAIR[7], Coefficients.cmAIR[8]
        elif specie == 'NO2':
            if 200 <= T < 1000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.clNO2[0], Coefficients.clNO2[1],  Coefficients.clNO2[2], \
                                                 Coefficients.clNO2[3], Coefficients.clNO2[4], Coefficients.clNO2[5], \
                                                 Coefficients.clNO2[6], Coefficients.clNO2[7], Coefficients.clNO2[8]
            elif 1000 <= T < 6000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.cmNO2[0], Coefficients.cmNO2[1],  Coefficients.cmNO2[2], \
                                                 Coefficients.cmNO2[3], Coefficients.cmNO2[4], Coefficients.cmNO2[5], \
                                                 Coefficients.cmNO2[6], Coefficients.cmNO2[7], Coefficients.cmNO2[8]
        elif specie == 'N2O':
            if 200 <= T < 1000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.clN2O[0], Coefficients.clN2O[1],  Coefficients.clN2O[2], \
                                                 Coefficients.clN2O[3], Coefficients.clN2O[4], Coefficients.clN2O[5], \
                                                 Coefficients.clN2O[6], Coefficients.clN2O[7], Coefficients.clN2O[8]
            elif 1000 <= T < 6000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.cmN2O[0], Coefficients.cmN2O[1],  Coefficients.cmN2O[2], \
                                                 Coefficients.cmN2O[3], Coefficients.cmN2O[4], Coefficients.cmN2O[5], \
                                                 Coefficients.cmN2O[6], Coefficients.cmN2O[7], Coefficients.cmN2O[8]
        elif specie == 'N2O3':
            if 200 <= T < 1000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.clN2O3[0], Coefficients.clN2O3[1],  Coefficients.clN2O3[2], \
                                                 Coefficients.clN2O3[3], Coefficients.clN2O3[4], Coefficients.clN2O3[5], \
                                                 Coefficients.clN2O3[6], Coefficients.clN2O3[7], Coefficients.clN2O3[8]
            elif 1000 <= T < 6000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.cmN2O3[0], Coefficients.cmN2O3[1],  Coefficients.cmN2O3[2], \
                                                 Coefficients.cmN2O3[3], Coefficients.cmN2O3[4], Coefficients.cmN2O3[5], \
                                                 Coefficients.cmN2O3[6], Coefficients.cmN2O3[7], Coefficients.cmN2O3[8]
        elif specie == 'N2O4':
            if 200 <= T < 1000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.clN2O4[0], Coefficients.clN2O4[1],  Coefficients.clN2O4[2], \
                                                 Coefficients.clN2O4[3], Coefficients.clN2O4[4], Coefficients.clN2O4[5], \
                                                 Coefficients.clN2O4[6], Coefficients.clN2O4[7], Coefficients.clN2O4[8]
            elif 1000 <= T < 6000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.cmN2O4[0], Coefficients.cmN2O4[1],  Coefficients.cmN2O4[2], \
                                                 Coefficients.cmN2O4[3], Coefficients.cmN2O4[4], Coefficients.cmN2O4[5], \
                                                 Coefficients.cmN2O4[6], Coefficients.cmN2O4[7], Coefficients.cmN2O4[8]
        elif specie == 'N2O5':
            if 200 <= T < 1000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.clN2O5[0], Coefficients.clN2O5[1],  Coefficients.clN2O5[2], \
                                                 Coefficients.clN2O5[3], Coefficients.clN2O5[4], Coefficients.clN2O5[5], \
                                                 Coefficients.clN2O5[6], Coefficients.clN2O5[7], Coefficients.clN2O5[8]
            elif 1000 <= T < 6000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.cmN2O5[0], Coefficients.cmN2O5[1],  Coefficients.cmN2O5[2], \
                                                 Coefficients.cmN2O5[3], Coefficients.cmN2O5[4], Coefficients.cmN2O5[5], \
                                                 Coefficients.cmN2O5[6], Coefficients.cmN2O5[7], Coefficients.cmN2O5[8]
        elif specie == 'HO2':
            if 200 <= T < 1000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.clHO2[0], Coefficients.clHO2[1],  Coefficients.clHO2[2], \
                                                 Coefficients.clHO2[3], Coefficients.clHO2[4], Coefficients.clHO2[5], \
                                                 Coefficients.clHO2[6], Coefficients.clHO2[7], Coefficients.clHO2[8]
            elif 1000 <= T < 6000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.cmHO2[0], Coefficients.cmHO2[1],  Coefficients.cmHO2[2], \
                                                 Coefficients.cmHO2[3], Coefficients.cmHO2[4], Coefficients.cmHO2[5], \
                                                 Coefficients.cmHO2[6], Coefficients.cmHO2[7], Coefficients.cmHO2[8]
        elif specie == 'H2O2':
            if 200 <= T < 1000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.clH2O2[0], Coefficients.clH2O2[1],  Coefficients.clH2O2[2], \
                                                 Coefficients.clH2O2[3], Coefficients.clH2O2[4], Coefficients.clH2O2[5], \
                                                 Coefficients.clH2O2[6], Coefficients.clH2O2[7], Coefficients.clH2O2[8]
            elif 1000 <= T < 6000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.cmH2O2[0], Coefficients.cmH2O2[1],  Coefficients.cmH2O2[2], \
                                                 Coefficients.cmH2O2[3], Coefficients.cmH2O2[4], Coefficients.cmH2O2[5], \
                                                 Coefficients.cmH2O2[6], Coefficients.cmH2O2[7], Coefficients.cmH2O2[8]
        elif specie == 'CH4':
            if 200 <= T < 1000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.clCH4[0], Coefficients.clCH4[1],  Coefficients.clCH4[2], \
                                                 Coefficients.clCH4[3], Coefficients.clCH4[4], Coefficients.clCH4[5], \
                                                 Coefficients.clCH4[6], Coefficients.clCH4[7], Coefficients.clCH4[8]
            elif 1000 <= T < 6000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.cmCH4[0], Coefficients.cmCH4[1],  Coefficients.cmCH4[2], \
                                                 Coefficients.cmCH4[3], Coefficients.cmCH4[4], Coefficients.cmCH4[5], \
                                                 Coefficients.cmCH4[6], Coefficients.cmCH4[7], Coefficients.cmCH4[8]
        elif specie == 'C2H2(ACETY)':
            if 200 <= T < 1000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.clC2H2_acety[0], Coefficients.clC2H2_acety[1],  Coefficients.clC2H2_acety[2], \
                                                 Coefficients.clC2H2_acety[3], Coefficients.clC2H2_acety[4], Coefficients.clC2H2_acety[5], \
                                                 Coefficients.clC2H2_acety[6], Coefficients.clC2H2_acety[7], Coefficients.clC2H2_acety[8]
            elif 1000 <= T < 6000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.cmC2H2_acety[0], Coefficients.cmC2H2_acety[1],  Coefficients.cmC2H2_acety[2], \
                                                 Coefficients.cmC2H2_acety[3], Coefficients.cmC2H2_acety[4], Coefficients.cmC2H2_acety[5], \
                                                 Coefficients.cmC2H2_acety[6], Coefficients.cmC2H2_acety[7], Coefficients.cmC2H2_acety[8]
        elif specie == 'C2H2(VINY)':
            if 200 <= T < 1000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.clC2H2_viny[0], Coefficients.clC2H2_viny[1],  Coefficients.clC2H2_viny[2], \
                                                 Coefficients.clC2H2_viny[3], Coefficients.clC2H2_viny[4], Coefficients.clC2H2_viny[5], \
                                                 Coefficients.clC2H2_viny[6], Coefficients.clC2H2_viny[7], Coefficients.clC2H2_viny[8]
            elif 1000 <= T < 6000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.cmC2H2_viny[0], Coefficients.cmC2H2_viny[1],  Coefficients.cmC2H2_viny[2], \
                                                 Coefficients.cmC2H2_viny[3], Coefficients.cmC2H2_viny[4], Coefficients.cmC2H2_viny[5], \
                                                 Coefficients.cmC2H2_viny[6], Coefficients.cmC2H2_viny[7], Coefficients.cmC2H2_viny[8]
        elif specie == 'C2H4':
            if 200 <= T < 1000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.clC2H4[0], Coefficients.clC2H4[1],  Coefficients.clC2H4[2], \
                                                 Coefficients.clC2H4[3], Coefficients.clC2H4[4], Coefficients.clC2H4[5], \
                                                 Coefficients.clC2H4[6], Coefficients.clC2H4[7], Coefficients.clC2H4[8]
            elif 1000 <= T < 6000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.cmC2H4[0], Coefficients.cmC2H4[1],  Coefficients.cmC2H4[2], \
                                                 Coefficients.cmC2H4[3], Coefficients.cmC2H4[4], Coefficients.cmC2H4[5], \
                                                 Coefficients.cmC2H4[6], Coefficients.cmC2H4[7], Coefficients.cmC2H4[8]
        elif specie == 'C2H6':
            if 200 <= T < 1000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.clC2H6[0], Coefficients.clC2H6[1],  Coefficients.clC2H6[2], \
                                                 Coefficients.clC2H6[3], Coefficients.clC2H6[4], Coefficients.clC2H6[5], \
                                                 Coefficients.clC2H6[6], Coefficients.clC2H6[7], Coefficients.clC2H6[8]
            elif 1000 <= T < 6000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.cmC2H6[0], Coefficients.cmC2H6[1],  Coefficients.cmC2H6[2], \
                                                 Coefficients.cmC2H6[3], Coefficients.cmC2H6[4], Coefficients.cmC2H6[5], \
                                                 Coefficients.cmC2H6[6], Coefficients.cmC2H6[7], Coefficients.cmC2H6[8]
        elif specie == 'C3H8':
            if 200 <= T < 1000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.clC3H8[0], Coefficients.clC3H8[1],  Coefficients.clC3H8[2], \
                                                 Coefficients.clC3H8[3], Coefficients.clC3H8[4], Coefficients.clC3H8[5], \
                                                 Coefficients.clC3H8[6], Coefficients.clC3H8[7], Coefficients.clC3H8[8]
            elif 1000 <= T < 6000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.cmC3H8[0], Coefficients.cmC3H8[1],  Coefficients.cmC3H8[2], \
                                                 Coefficients.cmC3H8[3], Coefficients.cmC3H8[4], Coefficients.cmC3H8[5], \
                                                 Coefficients.cmC3H8[6], Coefficients.cmC3H8[7], Coefficients.cmC3H8[8]
        elif specie == 'C4H8(ISOBUT)':
            if 200 <= T < 1000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.clC4H8_isobut[0], Coefficients.clC4H8_isobut[1],  Coefficients.clC4H8_isobut[2], \
                                                 Coefficients.clC4H8_isobut[3], Coefficients.clC4H8_isobut[4], Coefficients.clC4H8_isobut[5], \
                                                 Coefficients.clC4H8_isobut[6], Coefficients.clC4H8_isobut[7], Coefficients.clC4H8_isobut[8]
            elif 1000 <= T < 6000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.cmC4H8_isobut[0], Coefficients.cmC4H8_isobut[1],  Coefficients.cmC4H8_isobut[2], \
                                                 Coefficients.cmC4H8_isobut[3], Coefficients.cmC4H8_isobut[4], Coefficients.cmC4H8_isobut[5], \
                                                 Coefficients.cmC4H8_isobut[6], Coefficients.cmC4H8_isobut[7], Coefficients.cmC4H8_isobut[8]
        elif specie == 'C4H10(NBUT)':
            if 200 <= T < 1000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.clC4H10_nbut[0], Coefficients.clC4H10_nbut[1],  Coefficients.clC4H10_nbut[2], \
                                                 Coefficients.clC4H10_nbut[3], Coefficients.clC4H10_nbut[4], Coefficients.clC4H10_nbut[5], \
                                                 Coefficients.clC4H10_nbut[6], Coefficients.clC4H10_nbut[7], Coefficients.clC4H10_nbut[8]
            elif 1000 <= T < 6000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.cmC4H10_nbut[0], Coefficients.cmC4H10_nbut[1],  Coefficients.cmC4H10_nbut[2], \
                                                 Coefficients.cmC4H10_nbut[3], Coefficients.cmC4H10_nbut[4], Coefficients.cmC4H10_nbut[5], \
                                                 Coefficients.cmC4H10_nbut[6], Coefficients.cmC4H10_nbut[7], Coefficients.cmC4H10_nbut[8]
        elif specie == 'C4H10(ISOBUT)':
            if 200 <= T < 1000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.clC4H10_isobut[0], Coefficients.clC4H10_isobut[1],  Coefficients.clC4H10_isobut[2], \
                                                 Coefficients.clC4H10_isobut[3], Coefficients.clC4H10_isobut[4], Coefficients.clC4H10_isobut[5], \
                                                 Coefficients.clC4H10_isobut[6], Coefficients.clC4H10_isobut[7], Coefficients.clC4H10_isobut[8]
            elif 1000 <= T < 6000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.cmC4H10_isobut[0], Coefficients.cmC4H10_isobut[1],  Coefficients.cmC4H10_isobut[2], \
                                                 Coefficients.cmC4H10_isobut[3], Coefficients.cmC4H10_isobut[4], Coefficients.cmC4H10_isobut[5], \
                                                 Coefficients.cmC4H10_isobut[6], Coefficients.cmC4H10_isobut[7], Coefficients.cmC4H10_isobut[8]
        elif specie == 'C5H12(NPENT)':
            if 200 <= T < 1000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.clC5H12_npent[0], Coefficients.clC5H12_npent[1],  Coefficients.clC5H12_npent[2], \
                                                 Coefficients.clC5H12_npent[3], Coefficients.clC5H12_npent[4], Coefficients.clC5H12_npent[5], \
                                                 Coefficients.clC5H12_npent[6], Coefficients.clC5H12_npent[7], Coefficients.clC5H12_npent[8]
            elif 1000 <= T < 6000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.cmC5H12_npent[0], Coefficients.cmC5H12_npent[1],  Coefficients.cmC5H12_npent[2], \
                                                 Coefficients.cmC5H12_npent[3], Coefficients.cmC5H12_npent[4], Coefficients.cmC5H12_npent[5], \
                                                 Coefficients.cmC5H12_npent[6], Coefficients.cmC5H12_npent[7], Coefficients.cmC5H12_npent[8]
        elif specie == 'C6H14(NHEX)':
            if 200 <= T < 1000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.clC6H14_nhex[0], Coefficients.clC6H14_nhex[1],  Coefficients.clC6H14_nhex[2], \
                                                 Coefficients.clC6H14_nhex[3], Coefficients.clC6H14_nhex[4], Coefficients.clC6H14_nhex[5], \
                                                 Coefficients.clC6H14_nhex[6], Coefficients.clC6H14_nhex[7], Coefficients.clC6H14_nhex[8]
            elif 1000 <= T < 6000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.cmC6H14_nhex[0], Coefficients.cmC6H14_nhex[1],  Coefficients.cmC6H14_nhex[2], \
                                                 Coefficients.cmC6H14_nhex[3], Coefficients.cmC6H14_nhex[4], Coefficients.cmC6H14_nhex[5], \
                                                 Coefficients.cmC6H14_nhex[6], Coefficients.cmC6H14_nhex[7], Coefficients.cmC6H14_nhex[8]
        elif specie == 'C7H16(NHEP)':
            if 200 <= T < 1000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.clC7H16_nhep[0], Coefficients.clC7H16_nhep[1],  Coefficients.clC7H16_nhep[2], \
                                                 Coefficients.clC7H16_nhep[3], Coefficients.clC7H16_nhep[4], Coefficients.clC7H16_nhep[5], \
                                                 Coefficients.clC7H16_nhep[6], Coefficients.clC7H16_nhep[7], Coefficients.clC7H16_nhep[8]
            elif 1000 <= T < 6000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.cmC7H16_nhep[0], Coefficients.cmC7H16_nhep[1],  Coefficients.cmC7H16_nhep[2], \
                                                 Coefficients.cmC7H16_nhep[3], Coefficients.cmC7H16_nhep[4], Coefficients.cmC7H16_nhep[5], \
                                                 Coefficients.cmC7H16_nhep[6], Coefficients.cmC7H16_nhep[7], Coefficients.cmC7H16_nhep[8]
        elif specie == 'C7H16(2METH)':
            if 200 <= T < 1000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.clC7H16_2meth[0], Coefficients.clC7H16_2meth[1],  Coefficients.clC7H16_2meth[2], \
                                                 Coefficients.clC7H16_2meth[3], Coefficients.clC7H16_2meth[4], Coefficients.clC7H16_2meth[5], \
                                                 Coefficients.clC7H16_2meth[6], Coefficients.clC7H16_2meth[7], Coefficients.clC7H16_2meth[8]
            elif 1000 <= T < 6000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.cmC7H16_2meth[0], Coefficients.cmC7H16_2meth[1],  Coefficients.cmC7H16_2meth[2], \
                                                 Coefficients.cmC7H16_2meth[3], Coefficients.cmC7H16_2meth[4], Coefficients.cmC7H16_2meth[5], \
                                                 Coefficients.cmC7H16_2meth[6], Coefficients.cmC7H16_2meth[7], Coefficients.cmC7H16_2meth[8]
        elif specie == 'C8H18(NOCT)':
            if 200 <= T < 1000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.clC8H18_noct[0], Coefficients.clC8H18_noct[1],  Coefficients.clC8H18_noct[2], \
                                                 Coefficients.clC8H18_noct[3], Coefficients.clC8H18_noct[4], Coefficients.clC8H18_noct[5], \
                                                 Coefficients.clC8H18_noct[6], Coefficients.clC8H18_noct[7], Coefficients.clC8H18_noct[8]
            elif 1000 <= T < 6000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.cmC8H18_noct[0], Coefficients.cmC8H18_noct[1],  Coefficients.cmC8H18_noct[2], \
                                                 Coefficients.cmC8H18_noct[3], Coefficients.cmC8H18_noct[4], Coefficients.cmC8H18_noct[5], \
                                                 Coefficients.cmC8H18_noct[6], Coefficients.cmC8H18_noct[7], Coefficients.cmC8H18_noct[8]
        elif specie == 'C8H18(ISOCT)':
            if 200 <= T < 1000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.clC8H18_isoct[0], Coefficients.clC8H18_isoct[1],  Coefficients.clC8H18_isoct[2], \
                                                 Coefficients.clC8H18_isoct[3], Coefficients.clC8H18_isoct[4], Coefficients.clC8H18_isoct[5], \
                                                 Coefficients.clC8H18_isoct[6], Coefficients.clC8H18_isoct[7], Coefficients.clC8H18_isoct[8]
            elif 1000 <= T < 6000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.cmC8H18_isoct[0], Coefficients.cmC8H18_isoct[1],  Coefficients.cmC8H18_isoct[2], \
                                                 Coefficients.cmC8H18_isoct[3], Coefficients.cmC8H18_isoct[4], Coefficients.cmC8H18_isoct[5], \
                                                 Coefficients.cmC8H18_isoct[6], Coefficients.cmC8H18_isoct[7], Coefficients.cmC8H18_isoct[8]
        elif specie == 'NH':
            if 200 <= T < 1000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.clNH[0], Coefficients.clNH[1],  Coefficients.clNH[2], \
                                                 Coefficients.clNH[3], Coefficients.clNH[4], Coefficients.clNH[5], \
                                                 Coefficients.clNH[6], Coefficients.clNH[7], Coefficients.clNH[8]
            elif 1000 <= T < 6000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.cmNH[0], Coefficients.cmNH[1],  Coefficients.cmNH[2], \
                                                 Coefficients.cmNH[3], Coefficients.cmNH[4], Coefficients.cmNH[5], \
                                                 Coefficients.cmNH[6], Coefficients.cmNH[7], Coefficients.cmNH[8]
            elif 6000 <= T <= 20000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.chNH[0], Coefficients.chNH[1],  Coefficients.chNH[2], \
                                                 Coefficients.chNH[3], Coefficients.chNH[4], Coefficients.chNH[5], \
                                                 Coefficients.chNH[6], Coefficients.chNH[7], Coefficients.chNH[8]
        elif specie == 'NH2':
            if 200 <= T < 1000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.clNH2[0], Coefficients.clNH2[1],  Coefficients.clNH2[2], \
                                                 Coefficients.clNH2[3], Coefficients.clNH2[4], Coefficients.clNH2[5], \
                                                 Coefficients.clNH2[6], Coefficients.clNH2[7], Coefficients.clNH2[8]
            elif 1000 <= T < 6000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.cmNH2[0], Coefficients.cmNH2[1],  Coefficients.cmNH2[2], \
                                                 Coefficients.cmNH2[3], Coefficients.cmNH2[4], Coefficients.cmNH2[5], \
                                                 Coefficients.cmNH2[6], Coefficients.cmNH2[7], Coefficients.cmNH2[8]
        elif specie == 'NH3':
            if 200 <= T < 1000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.clNH3[0], Coefficients.clNH3[1],  Coefficients.clNH3[2], \
                                                 Coefficients.clNH3[3], Coefficients.clNH3[4], Coefficients.clNH3[5], \
                                                 Coefficients.clNH3[6], Coefficients.clNH3[7], Coefficients.clNH3[8]
            elif 1000 <= T < 6000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.cmNH3[0], Coefficients.cmNH3[1],  Coefficients.cmNH3[2], \
                                                 Coefficients.cmNH3[3], Coefficients.cmNH3[4], Coefficients.cmNH3[5], \
                                                 Coefficients.cmNH3[6], Coefficients.cmNH3[7], Coefficients.cmNH3[8]
        elif specie == 'N2H2':
            if 200 <= T < 1000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.clN2H2[0], Coefficients.clN2H2[1],  Coefficients.clN2H2[2], \
                                                 Coefficients.clN2H2[3], Coefficients.clN2H2[4], Coefficients.clN2H2[5], \
                                                 Coefficients.clN2H2[6], Coefficients.clN2H2[7], Coefficients.clN2H2[8]
            elif 1000 <= T < 6000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.cmN2H2[0], Coefficients.cmN2H2[1],  Coefficients.cmN2H2[2], \
                                                 Coefficients.cmN2H2[3], Coefficients.cmN2H2[4], Coefficients.cmN2H2[5], \
                                                 Coefficients.cmN2H2[6], Coefficients.cmN2H2[7], Coefficients.cmN2H2[8]
        elif specie == 'N2H4':
            if 200 <= T < 1000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.clN2H4[0], Coefficients.clN2H4[1],  Coefficients.clN2H4[2], \
                                                 Coefficients.clN2H4[3], Coefficients.clN2H4[4], Coefficients.clN2H4[5], \
                                                 Coefficients.clN2H4[6], Coefficients.clN2H4[7], Coefficients.clN2H4[8]
            elif 1000 <= T < 6000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.cmN2H4[0], Coefficients.cmN2H4[1],  Coefficients.cmN2H4[2], \
                                                 Coefficients.cmN2H4[3], Coefficients.cmN2H4[4], Coefficients.cmN2H4[5], \
                                                 Coefficients.cmN2H4[6], Coefficients.cmN2H4[7], Coefficients.cmN2H4[8]
        elif specie == 'N3H':
            if 200 <= T < 1000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.clN3H[0], Coefficients.clN3H[1],  Coefficients.clN3H[2], \
                                                 Coefficients.clN3H[3], Coefficients.clN3H[4], Coefficients.clN3H[5], \
                                                 Coefficients.clN3H[6], Coefficients.clN3H[7], Coefficients.clN3H[8]
            elif 1000 <= T < 6000:
                a1, a2, a3, a4, a5, a6, a7, b1, b2 = Coefficients.cmN3H[0], Coefficients.cmN3H[1],  Coefficients.cmN3H[2], \
                                                 Coefficients.cmN3H[3], Coefficients.cmN3H[4], Coefficients.cmN3H[5], \
                                                 Coefficients.cmN3H[6], Coefficients.cmN3H[7], Coefficients.cmN3H[8]

        self.Cp = (a1*self.T**-2 + a2*self.T**-1 + a3 +a4*self.T + a5*self.T**2 + a6*self.T**3 + a7*self.T**4)*R
        self.H = (- a1*self.T**-1 + a2*cmath.log(self.T, cmath.e) + a3*self.T + a4*self.T**2/2 + a5*self.T**3/3 + a6*self.T**4/4 + a7*self.T**5/5 + b1)*R
        self.S = (- a1*self.T**-2/2 - a2*self.T**-1 + a3*cmath.log(self.T, cmath.e) + a4*self.T + a5*self.T**2/2 + a6*self.T**3/3 + a7*self.T**4/4 + b2)*R


    def idealgas_func(self, termo_func: str)->float:
        """
        The idealgas_func method, from Termodata class, simply return the value of the specie heat capacity,
        enthalpy and entropy chosen by the user.

        Args:
            termo_func (str): Name of the phisical function you want to return 
                'heatcap_pconst', 'enthalpy' and 'entropy'

        Returns:
            float: Value for heat capacity, enthalpy and entropy
        """

        termo_func = termo_func.lower()
        if termo_func == 'heatcap_pconst':
            return self.Cp
        elif termo_func == 'enthalpy':
            return self.H
        elif termo_func == 'entropy':
            return self.S
