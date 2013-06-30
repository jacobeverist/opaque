// testNelmin.cpp : Defines the entry point for the console application.
//

#include <stdio.h>
#include "stdlib.h"
#include "stdafx.h"
#include "runNelmin.h"

#define NUM_PAIRS 34
#define NUM_POSES 100

int _tmain(int argc, _TCHAR* argv[])
{

	int numPairs = NUM_PAIRS;
	int numPoses = NUM_POSES;
	double u1 = 1.0;
	double currU = 0.37554483820209283;
	double currAng = -0.20129568662575242;
	double uHigh = 0.77300001382827754;
	double uLow = 0.37300001382827758;
	double c_matchPairs[408] = {-1.4381000951056415, -0.45711270964811412, 0.021087196732377626, 0.024099318814465719, 0.024099318814465719, 0.029912803267622377, 2.8290796554329649, 1.2322455083401485, 0.011730418727443509, 0.020264467216743243, 0.020264467216743243, 0.039269581272556489, -1.4381000951056415, -0.45711270964811412, 0.021087196732377626, 0.024099318814465719, 0.024099318814465719, 0.029912803267622377, 2.812131303656169, 1.1960135697769143, 0.011730418727443509, 0.020264467216743243, 0.020264467216743243, 0.039269581272556489, -1.4381000951056415, -0.45711270964811412, 0.021087196732377626, 0.024099318814465719, 0.024099318814465719, 0.029912803267622377, 2.7951829518793732, 1.1597816312136799, 0.011730418727443509, 0.020264467216743243, 0.020264467216743243, 0.039269581272556489, -1.4381000951056415, -0.45711270964811412, 0.021087196732377626, 0.024099318814465719, 0.024099318814465719, 0.029912803267622377, 2.7782346001025768, 1.1235496926504458, 0.011730418727443509, 0.020264467216743243, 0.020264467216743243, 0.039269581272556489, -1.4381000951056415, -0.45711270964811412, 0.021087196732377626, 0.024099318814465719, 0.024099318814465719, 0.029912803267622377, 2.761286248325781, 1.0873177540872117, 0.011730418727443509, 0.020264467216743243, 0.020264467216743243, 0.039269581272556489, -1.4381000951056415, -0.45711270964811412, 0.021087196732377626, 0.024099318814465719, 0.024099318814465719, 0.029912803267622377, 2.7443378965489851, 1.0510858155239773, 0.011730418727443509, 0.020264467216743243, 0.020264467216743243, 0.039269581272556489, -1.4381000951056415, -0.45711270964811412, 0.021087196732377626, 0.024099318814465719, 0.024099318814465719, 0.029912803267622377, 2.7366247474765091, 1.0345967582396205, 0.0043777679466941079, 0.012413754995419281, 0.012413754995419283, 0.046622232053305908, -1.4381000951056415, -0.45711270964811412, 0.021087196732377626, 0.024099318814465719, 0.024099318814465719, 0.029912803267622377, 2.7412484315893439, 0.99486488759117608, 0.0043777679466941079, 0.012413754995419281, 0.012413754995419283, 0.046622232053305908, -1.4381000951056415, -0.45711270964811412, 0.021087196732377626, 0.024099318814465719, 0.024099318814465719, 0.029912803267622377, 2.7458721157021784, 0.95513301694273167, 0.0043777679466941079, 0.012413754995419281, 0.012413754995419283, 0.046622232053305908, -1.4381000951056415, -0.45711270964811412, 0.021087196732377626, 0.024099318814465719, 0.024099318814465719, 0.029912803267622377, 2.7504957998150132, 0.91540114629428726, 0.0043777679466941079, 0.012413754995419281, 0.012413754995419283, 0.046622232053305908, -1.4381000951056415, -0.45711270964811412, 0.021087196732377626, 0.024099318814465719, 0.024099318814465719, 0.029912803267622377, 2.7551194839278481, 0.87566927564584285, 0.0043777679466941079, 0.012413754995419281, 0.012413754995419283, 0.046622232053305908, -1.4097655006153476, -0.42887885133808856, 0.021087196732377626, 0.024099318814465719, 0.024099318814465719, 0.029912803267622377, 2.7597431680406825, 0.83593740499739833, 0.0043777679466941079, 0.012413754995419281, 0.012413754995419283, 0.046622232053305908, -1.4097655006153476, -0.42887885133808856, 0.021087196732377626, 0.024099318814465719, 0.024099318814465719, 0.029912803267622377, 2.7643668521535174, 0.79620553434895391, 0.0043777679466941079, 0.012413754995419281, 0.012413754995419283, 0.046622232053305908, -1.3814309061250534, -0.400644993028063, 0.021087196732377626, 0.024099318814465719, 0.024099318814465719, 0.029912803267622377, 2.7689905362663523, 0.7564736637005095, 0.0043777679466941079, 0.012413754995419281, 0.012413754995419283, 0.046622232053305908, -1.3814309061250534, -0.400644993028063, 0.021087196732377626, 0.024099318814465719, 0.024099318814465719, 0.029912803267622377, 2.7736142203791871, 0.71674179305206509, 0.0043777679466941079, 0.012413754995419281, 0.012413754995419283, 0.046622232053305908, -1.3530963116347594, -0.37241113471803744, 0.021087196732377626, 0.024099318814465719, 0.024099318814465719, 0.029912803267622377, 2.7782379044920216, 0.67700992240362068, 0.0043777679466941079, 0.012413754995419281, 0.012413754995419283, 0.046622232053305908, -1.3530963116347594, -0.37241113471803744, 0.021087196732377626, 0.024099318814465719, 0.024099318814465719, 0.029912803267622377, 2.7828615886048564, 0.63727805175517616, 0.0043777679466941079, 0.012413754995419281, 0.012413754995419283, 0.046622232053305908, -1.3247617171444654, -0.34417727640801188, 0.021087196732377626, 0.024099318814465719, 0.024099318814465719, 0.029912803267622377, 2.7874852727176913, 0.59754618110673174, 0.0043777679466941079, 0.012413754995419281, 0.012413754995419283, 0.046622232053305908, -1.3247617171444654, -0.34417727640801188, 0.021087196732377626, 0.024099318814465719, 0.024099318814465719, 0.029912803267622377, 2.7921089568305257, 0.55781431045828733, 0.0043777679466941079, 0.012413754995419281, 0.012413754995419283, 0.046622232053305908, -1.3247617171444654, -0.34417727640801188, 0.021087196732377626, 0.024099318814465719, 0.024099318814465719, 0.029912803267622377, 2.7929106895488438, 0.55092492530967618, 0.016994415081616022, -0.022976183869306034, -0.022976183869306034, 0.034005584918383971, -1.2964271226541713, -0.31594341809798632, 0.021087196732377626, 0.024099318814465719, 0.024099318814465719, 0.029912803267622377, 2.8119906641699695, 0.51576877978066238, 0.016994415081616022, -0.022976183869306034, -0.022976183869306034, 0.034005584918383971, -1.2793997908124697, -0.29897662244537621, 0.032459176957488109, 0.023490846218779903, 0.023490846218779903, 0.018540823042511874, 2.8310706387910956, 0.48061263425164868, 0.016994415081616022, -0.022976183869306034, -0.022976183869306034, 0.034005584918383971, -1.2793997908124697, -0.29897662244537621, 0.032459176957488109, 0.023490846218779903, 0.023490846218779903, 0.018540823042511874, 2.8501506134122212, 0.44545648872263494, 0.016994415081616022, -0.022976183869306034, -0.022976183869306034, 0.034005584918383971, -1.24529755316485, -0.27807118326258762, 0.032459176957488109, 0.023490846218779903, 0.023490846218779903, 0.018540823042511874, 2.8692305880333473, 0.41030034319362119, 0.016994415081616022, -0.022976183869306034, -0.022976183869306034, 0.034005584918383971, -1.2111953155172304, -0.25716574407979909, 0.032459176957488109, 0.023490846218779903, 0.023490846218779903, 0.018540823042511874, 2.8883105626544729, 0.37514419766460738, 0.016994415081616022, -0.022976183869306034, -0.022976183869306034, 0.034005584918383971, -1.1770930778696107, -0.2362603048970105, 0.032459176957488109, 0.023490846218779903, 0.023490846218779903, 0.018540823042511874, 2.907390537275599, 0.33998805213559369, 0.016994415081616022, -0.022976183869306034, -0.022976183869306034, 0.034005584918383971, -1.142990840221991, -0.21535486571422197, 0.032459176957488109, 0.023490846218779903, 0.023490846218779903, 0.018540823042511874, 2.9264705118967247, 0.30483190660657988, 0.016994415081616022, -0.022976183869306034, -0.022976183869306034, 0.034005584918383971, -1.1088886025743712, -0.19444942653143338, 0.032459176957488109, 0.023490846218779903, 0.023490846218779903, 0.018540823042511874, 2.9455504865178508, 0.26967576107756613, 0.016994415081616022, -0.022976183869306034, -0.022976183869306034, 0.034005584918383971, -1.0888381755437173, -0.18215806350682748, 0.041797704452277193, 0.018293026802205917, 0.018293026802205917, 0.0092022955477228038, 2.9646304611389764, 0.23451961554855238, 0.016994415081616022, -0.022976183869306034, -0.022976183869306034, 0.034005584918383971, -1.0512444931790028, -0.16849369208839204, 0.041797704452277193, 0.018293026802205917, 0.018293026802205917, 0.0092022955477228038, 2.9837104357601021, 0.19936347001953858, 0.016994415081616022, -0.022976183869306034, -0.022976183869306034, 0.034005584918383971, -1.0136508108142883, -0.15482932066995658, 0.041797704452277193, 0.018293026802205917, 0.018293026802205917, 0.0092022955477228038, 3.0027904103812282, 0.16420732449052489, 0.016994415081616022, -0.022976183869306034, -0.022976183869306034, 0.034005584918383971, -0.97605712844957382, -0.14116494925152115, 0.041797704452277193, 0.018293026802205917, 0.018293026802205917, 0.0092022955477228038, 3.0218703850023538, 0.12905117896151114, 0.016994415081616022, -0.022976183869306034, -0.022976183869306034, 0.034005584918383971, -0.93846344608485932, -0.12750057783308572, 0.041797704452277193, 0.018293026802205917, 0.018293026802205917, 0.0092022955477228038, 3.0409503596234799, 0.093895033432497332, 0.016994415081616022, -0.022976183869306034, -0.022976183869306034, 0.034005584918383971, -0.90086976372014482, -0.11383620641465027, 0.041797704452277193, 0.018293026802205917, 0.018293026802205917, 0.0092022955477228038, 3.0600303342446056, 0.058738887903483528, 0.016994415081616022, -0.022976183869306034, -0.022976183869306034, 0.034005584918383971};
	double c_poses_1[300] = {6.31394, 10.0204, -1.8756, 6.281425, 9.9170, -1.882345, 6.2482700257896493, 9.8140944647065886, -1.8887090283257753, 6.2145125571664854, 9.7115114000787646, -1.8947072984514133, 6.1801856088437566, 9.6092673859868647, -1.9003437041159068, 6.145321714504167, 9.5073319031314014, -1.9056230641650309, 6.1099529473676277, 9.4056761514832754, -1.9105510845119502, 6.0741108916185933, 9.304273002252657, -1.915134269371896, 6.0378266138334213, 9.2030969498579402, -1.9193798370541333, 6.0011306344077253, 9.1021240638946512, -1.923295640814896, 5.9640528989837351, 9.0013319411043771, -1.9268900951184769, 5.9266227498776383, 8.9006996573436918, -1.9301721075161284, 5.8888688975069474, 8.8002077195530894, -1.9331510162359977, 5.8508193918178488, 8.6998380177258845, -1.9358365334773628, 5.8125015937125557, 8.5995737768771647, -1.9382386943205558, 5.7739421464766592, 8.4993995090126901, -1.9403678110987912, 5.7351669472064906, 8.3993009650978365, -1.9422344330256278, 5.6962011182364742, 8.2992650870265141, -1.9438493108328141, 5.6570689785664729, 8.1992799595900845, -1.9452233661443679, 5.6177940152891459, 8.0993347624462935, -1.9463676652938748, 5.5783988550173111, 7.9994197220881942, -1.9472933972796804, 5.5389052353112911, 7.8995260638130693, -1.9480118555471411, 5.4993339761062661, 7.7996459636913595, -1.9485344232863782, 5.4597049511396323, 7.6997725005355795, -1.9488725619373377, 5.4200370593783571, 7.5998996078692516, -1.9490378025995039, 5.3803481964463309, 7.5000220258958317, -1.9490417400516187, 5.3406552260517159, 7.400135253467619, -1.9488960290962367, 5.3009739514143126, 7.300235500054697, -1.9486123829535813, 5.2613190866929038, 7.2003196377138474, -1.9482025734391051, 5.2217042284126052, 7.1003851530574789, -1.9476784326689607, 5.1821418268922361, 7.0004300992225517, -1.9470518560475634, 5.1426431576716611, 6.9004530478395028, -1.9463348062985355, 5.1032182929391405, 6.8004530410011661, -1.9455393183086946, 5.0638760729586965, 6.7004295432316976, -1.9446775045611322, 5.0246240774974602, 6.6003823934555079, -1.9437615609374299, 4.985468597253023, 6.5003117569661768, -1.9428037726736214, 4.9464146052807987, 6.4002180773953778, -1.9418165202565609, 4.9074657284213732, 6.3001020286818159, -1.9408122850478837, 4.868624218727855, 6.1999644670401333, -1.9398036544224853, 4.8298909248932356, 6.099806382929847, -1.9388033262068565, 4.7912652636777446, 5.9996288530242667, -1.9378241121992781, 4.7527451913361922, 5.8994329921794249, -1.9368789405490878, 4.7143271750453382, 5.799219905402996, -1.9359808567683836, 4.6760061643312376, 5.6989906398232257, -1.9351430231414783, 4.6377755624965964, 5.5987461366578488, -1.9343787162918618, 4.5996271980481218, 5.4984871831830171, -1.9337013226580482, 4.561551296123886, 5.3982143647022252, -1.9331243316227729, 4.5235364499206749, 5.2979280165152405, -1.932661326031526, 4.485569592121335, 5.1976281758870115, -1.9323259698298192, 4.4476359663221405, 5.0973145340166042, -1.9321319925426859, 4.4097190984601387, 4.9969863880061265, -1.9320931568582456, 4.371800770150001, 4.8966425937994575, -1.932222492050768, 4.3338610976268912, 4.796281572756885, -1.9325293753971886, 4.2958790218412135, 4.6959015270483953, -1.9330178027523788, 4.2578330428450206, 4.595500779119476, -1.9336862787014273, 4.2197019715448736, 4.4950781199771095, -1.9345278505397316, 4.1814656814578566, 4.3946331574773891, -1.9355301624656305, 4.1431058604675917, 4.2941666646131065, -1.9366755266360702, 4.1046067625802474, 4.1936809278013598, -1.9379410057390236, 4.0659559596805561, 4.0931800951711521, -1.9392984997645646, 4.027145093287829, 3.9926705248509964, -1.9407148277228581, 3.9881706263119647, 3.8921611332565127, -1.9421517931458971, 3.9490345948094681, 3.7916637433780296, -1.9435662202818129, 3.9097453597394596, 3.6911934330681917, -1.9449103212570586, 3.8703183096126859, 3.5907688690654593, -1.9461403578720176, 3.8307754198016695, 3.4904123089124082, -1.9472371433804234, 3.7911421054926162, 3.390148816759154, -1.9482134538025353, 3.7514430809337265, 3.2900051905260597, -1.9491142627119988, 3.7116981845878265, 3.1900088791626939, -1.9500168989206517, 3.6719182042849829, 3.0901868999067688, -1.951031340282634, 3.6321007023751282, 2.9905647555430965, -1.9523005413943733, 3.5922258408806766, 2.8911653516625324, -1.9540006350217412, 3.5522522066491526, 2.7920079139209264, -1.9563407843793583, 3.5121126365058029, 2.6931069052980678, -1.9595623950428762, 3.4717100424062282, 2.5944709433566411, -1.9639373214763418, 3.4309132365889967, 2.4961017175011673, -1.9697645396072503, 3.3895527670390178, 2.3979929069671329, -1.9773448266972085, 3.347419160124431, 2.3001292692794997, -1.9868326353574493, 3.3042781625439868, 2.2024869325554177, -1.9980943036514403, 3.2599029140483657, 2.1050348867975943, -2.0106979038797821, 3.2141084939565689, 2.0077366432932635, -2.0239320972786929, 3.1667864775380266, 1.9105518947118192, -2.0368310634301445, 3.117939492394691, 1.8134381752024546, -2.0481985871314965, 3.0677157748431338, 1.7163525204917913, -2.0566209572729499, 3.016443726296647, 1.6192531279815194, -2.0604554754311279, 2.9646664696473333, 1.5221010168460234, -2.0577792492673583, 2.9131764056482137, 1.4248616881300267, -2.0462813649316955, 2.863049769295317, 1.3275067848462208, -2.0230812755367467, 2.8156811862097806, 1.2300157520728996, -1.9844626277271271, 2.7728182829952859, 1.1323775058849093, -1.9253081032477146, 2.7366247474765091, 1.0345967582396205, -1.836463327769742, 2.7099986144482573, 0.93674211070184155, -1.7010766046893322, 2.6971980353577294, 0.8390444697615359, -1.4973009982349614, 2.7043443792299935, 0.74198425881804009, -1.2425965616963404, 2.7370120717737061, 0.64604802954395379, -1.0395109109482232, 2.7929106895488438, 0.55092492530967618, -0.96267227953713053, 2.8595270634631622, 0.45522964367845942, -1.0088785941540768, 2.9212203121612661, 0.35724541645337032, -1.135264792085076, 2.9683842896262709, 0.25589011756793684, -1.2189802437431847, 3.0066366415056267, 0.15168524275119885, -1.0519018421160453};
	double c_poses_2[300] = {-2.00106, -1.5116, 0.64040, -1.966762, -1.4860811, 0.760263, -1.9355664498879441, -1.456415545610011, 0.86331538355978366, -1.9071017381708562, -1.4231296365840809, 0.949505584683624, -1.881007807360894, -1.3866784409876545, 1.0200993257276185, -1.8569448867601943, -1.3474921410620604, 1.0769047676654462, -1.8345930737145224, -1.3059766774388499, 1.121794181190344, -1.8136518348413813, -1.2625143404583803, 1.1564691738048587, -1.7938395072581286, -1.2174643614884064, 1.1823768225730542, -1.7748927998100896, -1.1711635042426665, 1.2007029926953634, -1.7565662942986711, -1.123926656099471, 1.2123986316620388, -1.7386319467094786, -1.0760474194202934, 1.2182163223376146, -1.7208785884404281, -1.027798702868356, 1.2187468587410597, -1.703111427529862, -0.97943331272722023, 1.2144520161279797, -1.6851515498846605, -0.93118454421937291, 1.2056926802381791, -1.6668354205083611, -0.88326677282481647, 1.1927527492407013, -1.6480143847292676, -0.835876045599658, 1.1758596339195297, -1.6285541694285679, -0.78919067249469443, 1.1552022109343159, -1.6083343842684477, -0.74337181767400451, 1.1309469425549465, -1.5872480229202028, -0.69856409083353432, 1.1032526559864095, -1.5652009642923561, -0.65489613851968875, 1.0722842137943156, -1.5421114737587718, -0.61248123544791744, 1.0382250227881367, -1.5179097043867673, -0.57141787582130399, 1.0012880422961887, -1.4925371981652309, -0.5317903646491543, 0.9617246971038883, -1.4659463872327325, -0.49366940906558565, 0.9198309246043761, -1.4381000951056415, -0.45711270964811412, 0.87594954912468204, -1.4089710379062395, -0.4221655517362442, 0.83046833133870657, -1.3785413255908336, -0.38886139675005554, 0.78381340788412712, -1.3468019631778727, -0.35722247350879266, 0.73643838092026959, -1.3137523519760599, -0.32726036954945331, 0.68880993983547723, -1.2793997908124697, -0.29897662244537632, 0.64139145073922132, -1.2437589772606601, -0.27236331112483059, 0.59462628371192128, -1.2068515088687857, -0.2474036471896027, 0.54892266472399442, -1.1687053843877169, -0.22407256623358668, 0.50464153140323309, -1.129354504999148, -0.20233731916137082, 0.46208832599428246, -1.0888381755437173, -0.18215806350682748, 0.42150902006398688, -1.047200605749119, -0.1634884547517006, 0.383090083058227, -1.0044904114582145, -0.14627623764419442, 0.34696168709917891, -0.96076011585715182, -0.13046383751756202, 0.31320322541231405, -0.91606565070347845, -0.11598895160869369, 0.28185019501866526, -0.87046585755425454, -0.10278514037670534, 0.25290160359759611, -0.82402198899416668, -0.090782418821526795, 0.22632724325529108, -0.77679720986364542, -0.079907847802490301, 0.20207437627001576, -0.72885609848697563, -0.070086125356918827, 0.18007356282706788, -0.68026414790041456, -0.061240178018714862, 0.16024350917796396, -0.63108726708030316, -0.05329175213694845, 0.14249492131353378, -0.5813912821711823, -0.046162005194445847, 0.1267334175638242, -0.53124143771390719, -0.039772097126377637, 0.11286159110299544, -0.48070189787375989, -0.034043781638847553, 0.10078032870426826, -0.42983524766856562, -0.028899997527480646, 0.090389493030775009, -0.37870199419680667, -0.024265459996011651, 0.081588068505663344, -0.32736006786573579, -0.020067251974873693, 0.074273859981873205, -0.27586432361949192, -0.01623541543978638, 0.068342822200485639, -0.22426604216721327, -0.012703542730344526, 0.06368808846677175, -0.17261243121219666, -0.0094093678674048529, 0.060199078967586153, -0.12094614032882631, -0.0062953421639754909, 0.057765095411477119, -0.069304989649642754, -0.0033089637733270086, 0.056282869922715757, -0.017722086746104239, -0.00040266048551917955, 0.055657668206407947, 0.033773948304409662, 0.0024664519818943316, 0.055802143828743989, 0.085158755563482499, 0.0053368143264310814, 0.056635354401293407, 0.1364120078765409, 0.0082426680032447169, 0.058081968427286627, 0.18751718351773869, 0.011214299567689977, 0.060071640724373342, 0.23846133883483933, 0.014278285017887787, 0.062538536705216557, 0.28923488089409888, 0.017457734137290175, 0.06542098773306243, 0.33983134012514793, 0.020772534837245289, 0.068661261286858288, 0.39024714296587509, 0.024239597499562473, 0.072205430823340397, 0.44048138450730856, 0.02787309931907709, 0.076003331086489184, 0.49053560113849948, 0.031684728646215694, 0.080008585273723815, 0.54041354319140422, 0.035683929329560923, 0.084178691009633946, 0.59012094758576783, 0.03987814505841656, 0.088475152584275699, 0.63966531047400532, 0.044273063705372449, 0.092863647456917428, 0.68905565988608464, 0.048872861668869567, 0.097314215667788709, 0.73830232837441112, 0.053680448215765071, 0.10180146158420994, 0.78741672565870757, 0.058697709823897125, 0.10630475836198501, 0.83641111127089818, 0.063925754524650102, 0.11080844663904894, 0.88529836719999067, 0.06936515624551931, 0.11530202029184458, 0.93409177053695935, 0.075016199152676308, 0.1197802935564298, 0.9828047661196273, 0.080879121993533656, 0.12424354541609277, 1.03145073917755, 0.086954362439310109, 0.12869763884555258, 1.0800427879768968, 0.093242801427595418, 0.13315411423418461, 1.1285934964653339, 0.099746007504915468, 0.13763025803996629, 1.1771147069169074, 0.10646648116929723, 0.14214914940437895, 1.2256172925769251, 0.11340789921283376, 0.1467396890417885, 1.2741109303068394, 0.12057535906424913, 0.15143661616312087, 1.3226038732291312, 0.12797562313146355, 0.15628052046605706, 1.3711027233721906, 0.1356173631441584, 0.16131785728907891, 1.4196122043152006, 0.14351140449634098, 0.16660097485424832, 1.4681349338330203, 0.15167097058890977, 0.17218816308244472, 1.5166711965410655, 0.16011192717221928, 0.17814373372111961, 1.5652187165401923, 0.16885302668864499, 0.18453814143365468, 1.6137724300615823, 0.17791615261514865, 0.19144815499933387, 1.66232425811162, 0.18732656380584284, 0.19895708677306342, 1.7108628791167808, 0.19711313883455642, 0.20715508692060691, 1.7593735015685099, 0.20730862033739919, 0.21613950648197836, 1.8078376366681064, 0.21794985935532701, 0.22601532974133057, 1.8562328709716063, 0.22907805967670689, 0.23689567129697059, 1.9045326390346624, 0.24073902217988147, 0.24890232607949034, 1.9527059960574336, 0.25298338917573532, 0.26216635061512344, 2.0007173905294589, 0.26586688875025827, 0.27682864009837826, 2.0485264368745448, 0.27945057910711135, 0.2930404470851638};
	double *resultParam;
	double *resultSum;
	double *initGuess;
	double *resultOffset;
	int i, j;

	initGuess = (double*) malloc(3*sizeof(double));
	initGuess[0] = u1;
	initGuess[1] = currU;
	initGuess[2] = currAng;

	resultSum = (double*) malloc(sizeof(double));
	resultParam = (double*) malloc(2*sizeof(double));
	resultOffset = (double*) malloc(3*sizeof(double));

	doTest(c_matchPairs, numPairs, initGuess, uHigh, uLow, c_poses_1, c_poses_2, numPoses, resultParam, resultSum);
	//doCost(c_matchPairs, numPairs, initGuess, uHigh, uLow, c_poses_1, c_poses_2, numPoses, resultParam, resultSum, resultOffset);

	printf("%lf %lf %lf\n", resultSum[0], resultParam[0], resultParam[1]);
	//printf("%lf %lf %lf\n", resultOffset[0], resultOffset[1], resultOffset[2]);

	/*
	for ( i = 0  ; i < numPoses ; i++ ) {
		printf("%d:  ", i);
		for ( j = 0 ; j < 3 ; j++ ) {
			printf("%lf ", c_poses_1[i*3 + j]);
		}
		for ( j = 0 ; j < 3 ; j++ ) {
			printf("%lf ", c_poses_2[i*3 + j]);
		}
		printf("\n");
	}
	*/

	free(resultSum);
	free(resultParam);
	free(initGuess);
	return 0;
}
