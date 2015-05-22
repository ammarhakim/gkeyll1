import numpy

def makeMatrix(*ll):
    nrow = len(ll)
    ncol = len(ll[0])
    mat = numpy.zeros((nrow,ncol), numpy.float)
    for i in range(nrow):
        mat[i,:] = ll[i]
    return mat

## Some notes: The interpolation coefficients are generated in Maxima
## and cut-paste here. Perhaps it may be easier and more compact to
## actually write the Python code to generate these on the
## fly. However, in higher dimensions it may be pretty slow to
## generate these every time.

#################
class GkeDgLobatto1DPolyOrder1Basis:
    cMat_i2 = makeMatrix([0.75,0.25],[0.25,0.75])

#################
class GkeDgLobatto1DPolyOrder2Basis:
    cMat_i3 = makeMatrix([.5555555555555556,.5555555555555556,-.1111111111111111],[0.0,1.0,0.0],[-.1111111111111111,.5555555555555556,.5555555555555556])

#################
class GkeDgLobatto1DPolyOrder3Basis:
    cMat_i4 = makeMatrix([.3964843749999999,.7320061281981992,-.1851311281981993,.05664062500000011],[-0.107421875,.9134865201415714,.2583884798584289,-.06445312499999978],[-.06445312500000017,.2583884798584296,.9134865201415707,-.1074218749999999],[.05664062499999961,-.1851311281981994,.7320061281981994,.3964843749999998])

#################
class GkeDgLobatto1DPolyOrder4Basis:
    cMat_i5 = makeMatrix([.2663999999999975,.8553363583762964,-0.177600000000006,.08546364162371375,-.02960000000000157],[-.1316000000000021,.7234924181056768,.5263999999999952,-.1746924181056689,.05639999999999855],[2.775557561562891e-17,1.110223024625157e-16,1.0,-5.551115123125783e-16,-1.110223024625157e-16],[.05640000000000095,-.1746924181056717,.5264000000000018,.7234924181056706,-.1315999999999997],[-.02959999999999938,0.0854636416237109,-.1775999999999993,.8553363583762902,.2663999999999995])

#################
class GkeDgLobatto2DPolyOrder1Basis:
    cMat_i2 = makeMatrix([0.5625,0.1875,0.1875,0.0625],[0.1875,0.0625,0.5625,0.1875],[0.1875,0.5625,0.0625,0.1875],[0.0625,0.1875,0.1875,0.5625])

#################
class GkeDgLobatto2DPolyOrder2Basis:
    cMat_i3 = makeMatrix([0.308641975308642,0.308641975308642,-.06172839506172839,0.308641975308642,0.308641975308642,-.06172839506172839,-.06172839506172839,-.06172839506172839,.01234567901234568],[0.0,3.700743415417188e-17,0.0,.5555555555555556,.5555555555555556,-.1111111111111111,0.0,0.0,0.0],[-.06172839506172839,-.06172839506172836,.01234567901234568,0.308641975308642,0.308641975308642,-.06172839506172839,0.308641975308642,0.308641975308642,-.06172839506172839],[0.0,.5555555555555556,0.0,3.700743415417188e-17,.5555555555555556,0.0,0.0,-.1111111111111111,0.0],[0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0],[0.0,-.1111111111111111,0.0,-3.700743415417188e-17,.5555555555555556,0.0,0.0,.5555555555555556,0.0],[-.06172839506172839,0.308641975308642,0.308641975308642,-.06172839506172836,0.308641975308642,0.308641975308642,.01234567901234568,-.06172839506172839,-.06172839506172839],[0.0,-3.700743415417188e-17,0.0,-.1111111111111111,.5555555555555556,.5555555555555556,0.0,0.0,0.0],[.01234567901234568,-.06172839506172843,-.06172839506172839,-.06172839506172843,0.308641975308642,0.308641975308642,-.06172839506172839,0.308641975308642,0.308641975308642])

#################
class GkeDgLobatto2DPolyOrder3Basis:
    cMat_i4 = makeMatrix([.1571998596191407,.2902289922348328,-.07340159965670805,0.0224571228027344,.2902289922348327,.5358329717197183,-.1355171203613287,.04146128460497624,-.07340159965670807,-.1355171203613285,.03427353462793812,-.01048594280810112,.02245712280273439,.04146128460497629,-0.0104859428081009,.003208160400390668],[-.04259109497070331,-.07863347080254104,.01988713291191598,-.006084442138671904,.3621831320092554,.6686777307700775,-.1691147900676565,.05174044742989372,.1024469949438692,.1891419507121868,-.04783575078960892,.01463528499198142,-.02555465698242185,-0.0471800824815243,0.0119322797471494,-.003650665283203144],[-0.0255546569824223,-.04718008248152482,.01193227974714984,-.003650665283203138,.1024469949438695,.1891419507121866,-.04783575078960881,.01463528499198118,.3621831320092555,0.668677730770077,-.1691147900676563,.05174044742989382,-.04259109497070321,-.07863347080254075,.01988713291191576,-.006084442138671863],[.02245712280273402,.04146128460497608,-.01048594280810097,.003208160400390558,-.07340159965670794,-0.135517120361329,0.0342735346279384,-.01048594280810136,.2902289922348328,.5358329717197182,-.1355171203613287,.04146128460497649,.1571998596191403,.2902289922348334,-.07340159965670784,.02245712280273433],[-.04259109497070318,.3621831320092558,.1024469949438693,-.02555465698242184,-.07863347080254117,.6686777307700776,.1891419507121867,-.04718008248152458,.01988713291191599,-.1691147900676562,-.04783575078960876,.01193227974714961,-.006084442138671804,.05174044742989376,.01463528499198156,-.003650665283203205],[.01153945922851552,-0.0981284347808325,-.02775657498479214,.006923675537109358,-.09812843478083325,.8344576224803566,.2360343933105466,-.05887706086849967,-.02775657498479209,0.236034393310547,.06676460652354946,-.01665394499087517,.006923675537109367,-.05887706086849945,-.01665394499087542,.004154205322265587],[.006923675537109233,-.05887706086849953,-.01665394499087519,.004154205322265624,-.02775657498479295,.2360343933105468,.06676460652355028,-.01665394499087518,-.09812843478083277,.8344576224803563,.2360343933105467,-.05887706086849968,.01153945922851541,-.09812843478083273,-.02775657498479235,0.00692367553710938],[-.006084442138672251,0.0517404474298937,.01463528499198139,-0.00365066528320318,0.0198871329119148,-.1691147900676563,-.04783575078960812,.01193227974714997,-0.0786334708025413,.6686777307700777,0.189141950712187,-.04718008248152454,-.04259109497070367,.3621831320092557,.1024469949438693,-.02555465698242183],[-.02555465698242213,.1024469949438698,.3621831320092556,-.04259109497070312,-.04718008248152474,.1891419507121872,.6686777307700772,-0.0786334708025412,.01193227974714993,-.04783575078960897,-.1691147900676566,.01988713291191614,-.003650665283203055,.01463528499198124,.05174044742989373,-.006084442138671915],[.006923675537109431,-.02775657498479143,-.09812843478083291,.01153945922851558,-0.0588770608685005,.2360343933105468,.8344576224803563,-.09812843478083284,-.01665394499087489,.06676460652354996,.2360343933105462,-.02775657498479199,.004154205322265519,-.01665394499087555,-.05887706086849948,.006923675537109403],[.004154205322266276,-.01665394499087465,-.05887706086849976,.006923675537109343,-.01665394499087727,.06676460652354975,.2360343933105472,-.02775657498479176,-.05887706086849955,.2360343933105468,.8344576224803568,-.09812843478083295,.006923675537108965,-.02775657498479227,-.09812843478083273,.01153945922851563],[-.003650665283202641,.01463528499198117,.05174044742989377,-.006084442138671961,0.0119322797471468,-.04783575078960856,-0.169114790067656,.01988713291191711,-.04718008248152522,.1891419507121875,.6686777307700791,-.07863347080254138,-.02555465698242272,.1024469949438693,.3621831320092554,-.04259109497070305],[.02245712280273428,-.07340159965670776,.2902289922348332,.1571998596191407,.04146128460497531,-.1355171203613274,0.535832971719718,.2902289922348328,-.01048594280810084,.03427353462793805,-.1355171203613295,-.07340159965670796,.003208160400390632,-.01048594280810129,.04146128460497597,.02245712280273444],[-.006084442138671532,.01988713291191624,-.07863347080254093,-.04259109497070316,.05174044742989217,-.1691147900676554,.6686777307700774,.3621831320092558,.01463528499198182,-.04783575078960819,.1891419507121864,.1024469949438694,-0.00365066528320339,0.0119322797471493,-.04718008248152417,-.02555465698242179],[-.003650665283201892,.01193227974714944,-.04718008248152489,-0.0255546569824219,.01463528499197896,-.04783575078960874,.1891419507121876,.1024469949438699,.05174044742989418,-.1691147900676555,.6686777307700781,.3621831320092556,-.006084442138672542,.01988713291191575,-.07863347080254038,-0.0425910949707032],[.003208160400392545,-.01048594280810275,.04146128460497584,0.0224571228027343,-.01048594280810436,.03427353462793804,-.1355171203613276,-0.0734015996567064,.04146128460497557,-.1355171203613278,.5358329717197208,.2902289922348325,.02245712280273318,-.07340159965670842,.2902289922348334,.1571998596191406])

#################
class GkeDgLobatto2DPolyOrder4Basis:
    cMat_i5 = makeMatrix([.07096895999999885,.2278616058714436,-.04731264000000123,.02276751412855717,-.007885440000000328,.2278616058714438,.7316002859604237,-.1519077372476352,.07310016000000394,-.02531795620793976,-.04731264000000118,-.1519077372476355,.03154176000000217,-.01517834275237207,.005256960000000446,.02276751412855704,.07310016000000423,-.01517834275237238,.007304034039586597,-.002529723792061912,-.007885440000000321,-.02531795620793973,.005256960000000335,-.002529723792062067,8.761600000000652e-4],[-.03505823999999992,-0.112562264762322,0.0233721600000012,-.01124701523768097,.003895360000000297,.1927383801833501,.6188293702153704,-.1284922534555725,.06183229673845706,-.02141537557592904,.1402329599999976,.4502490590492774,-.09348864000000268,.04498806095072296,-.01558144000000069,-0.04653806018335,-.1494207767384517,.03102537345556795,-.01492985021536294,.005170895575928046,.01502495999999952,.04824097061242193,-.01001664000000023,.004820149387577149,-0.00166944],[4.874692427290324e-16,6.777357974399691e-16,1.732127615185474e-16,3.137551000242133e-17,2.557329902629354e-17,-1.062096714399712e-15,2.910415832948089e-16,2.679918385318459e-16,-2.89144421535523e-16,9.982345831032652e-17,.2663999999999981,.8553363583762955,-.1776000000000064,.08546364162371439,-.02960000000000168,-1.960845491077117e-16,7.188652539266563e-16,5.322708374108828e-17,-4.624906611956046e-16,1.184391621454331e-16,-9.718625852890279e-18,-1.695616808550743e-17,-2.141560729598039e-16,-4.918250383626351e-17,-8.039215960458653e-17],[.01502496000000069,.04824097061242465,-.01001664000000037,0.00482014938757766,-.001669440000000097,-.04653806018335165,-.1494207767384557,.03102537345556896,-.01492985021536331,.005170895575928228,.1402329599999998,.4502490590492836,-.09348864000000291,.04498806095072304,-.01558144000000105,.1927383801833494,.6188293702153669,-.1284922534555723,.06183229673845671,-.02141537557592861,-.03505823999999967,-.1125622647623205,.02337216000000037,-.01124701523768036,.003895359999999887],[-0.00788543999999879,-.02531795620793804,.005256960000000225,-.002529723792061741,8.761599999999625e-4,.02276751412855474,.07310016000000165,-.01517834275237124,.007304034039586257,-.002529723792061973,-.04731263999999943,-.1519077372476298,.03154176000000163,-.01517834275237242,.005256960000000597,.2278616058714427,.7316002859604207,-.1519077372476356,.07310016000000451,-.02531795620793999,.07096895999999889,.2278616058714449,-.04731264000000198,.02276751412855701,-.007885440000000321],[-.03505824000000029,.1927383801833506,.1402329599999979,-.04653806018334997,0.0150249599999996,-.1125622647623222,.6188293702153705,.4502490590492777,-.1494207767384521,.04824097061242204,.02337216000000109,-.1284922534555729,-.09348864000000219,.03102537345556793,-.01001664000000019,-.01124701523768085,.06183229673845734,.04498806095072218,-.01492985021536296,.004820149387577304,.003895360000000148,-.02141537557592901,-.01558144000000086,.005170895575927936,-.001669440000000057],[.01731856000000063,-.09521160222270748,-.06927424000000056,.02298952222270606,-0.00742223999999988,-.09521160222270801,0.523441279056398,0.380846408890825,-.1263886399999983,.04080497238115912,-.06927424000000018,.3808464088908251,0.277096959999995,-.09195808889082344,.02968895999999939,.02298952222270598,-.1263886399999973,-.09195808889082414,.03051744094360613,-0.0098526523811596,-.007422239999999703,.04080497238115893,.02968895999999911,-.009852652381159916,0.00318095999999977],[-1.067553715400607e-16,2.096003485637495e-15,-5.91378786276663e-16,-1.810562543624595e-16,1.194382703799655e-16,2.989354942576149e-16,-2.610789562633086e-15,1.802489414230277e-16,-8.43649072312012e-17,3.946914823923329e-17,-0.131600000000002,.7234924181056764,.5263999999999962,-0.174692418105669,.05639999999999901,-4.875826177476976e-17,2.164323156037239e-15,-2.11397918844838e-15,1.05141529997695e-15,-2.776496041605517e-16,3.872997922359546e-16,-5.137958835802737e-16,3.906999447564368e-16,-9.515058794014654e-16,-5.746814158100204e-18],[-.007422240000000881,.04080497238116312,.02968895999999987,-.009852652381160372,0.00318096000000036,.02298952222270649,-.1263886400000037,-.09195808889082681,0.030517440943607,-.009852652381159604,-.06927424000000212,.3808464088908266,.2770969600000008,-.09195808889082385,0.0296889599999992,-.09521160222270614,.5234412790563981,.3808464088908191,-.1263886399999956,0.0408049723811585,.01731856000000033,-.09521160222270707,-.06927423999999896,.02298952222270456,-.007422239999999619],[.003895359999999418,-.02141537557592463,-.01558144000000006,.005170895575925955,-0.00166943999999906,-.01124701523767967,.06183229673845228,.04498806095071929,-.01492985021535944,.004820149387576108,0.0233721599999987,-.1284922534555709,-.09348863999999578,.03102537345556647,-.01001663999999937,-0.112562264762318,.6188293702153697,.4502490590492727,-.1494207767384497,.04824097061242123,-.03505824000000018,.1927383801833518,.1402329599999979,-.04653806018335085,.01502495999999957],[-3.847884678950865e-16,5.742387388817286e-16,.2663999999999982,-2.676706939334655e-16,2.083708762299356e-17,4.291770424355782e-16,-7.374600022263116e-16,.8553363583762963,-1.114494814437701e-16,1.826802347907917e-16,3.822019147792809e-16,-2.854244873972616e-16,-.1776000000000059,4.25544261564614e-16,-2.634606044871052e-16,1.576304264778235e-16,-5.560422650415168e-16,.08546364162371334,-3.861768610833561e-16,1.670355789872684e-16,-1.525114257849683e-16,3.238735394008555e-16,-.02960000000000189,-4.627978071665251e-17,-9.929044257199922e-17],[-4.918164099131217e-16,2.006024163300577e-15,-.1316000000000022,-5.092817417930436e-16,1.595242049599236e-16,8.735924436589949e-16,-3.597368236001792e-15,.7234924181056782,-7.768230227048496e-18,-9.477184019782624e-17,4.619847701530193e-16,1.187273568668316e-15,0.526399999999995,6.853120128359735e-16,3.521466450986568e-17,-2.107942854301935e-16,7.263229825269456e-16,-.1746924181056703,3.622561476968433e-16,5.624761519569559e-17,3.042405438116214e-16,-4.805301319412178e-16,.05639999999999879,-7.002942409404135e-16,-3.031217503422521e-17],[-3.621991595913266e-16,3.365305951800332e-15,-8.531510203353255e-16,-7.108465402362325e-16,3.149702720802255e-16,5.152632901933026e-16,-5.370371653652188e-15,7.144039334471126e-16,1.551707583048934e-16,-1.381247568723915e-16,-2.321298160396437e-17,4.112295232311789e-17,1.000000000000001,1.129077841305945e-15,-6.265440840264258e-17,1.953945144974528e-16,2.847994812836151e-15,-3.200564689814534e-15,1.469051770354068e-15,5.916757055518943e-17,8.328226996457407e-16,-1.174954075917282e-15,5.257498140553746e-16,-1.321363287921881e-15,-8.649192477787416e-17],[2.663739895488464e-16,3.915028251117109e-15,.05639999999999985,-1.085492535674946e-15,4.761816829583767e-16,-7.308486702356176e-16,-5.484139716372857e-15,-.1746924181056752,5.111269773010277e-16,3.816879816209253e-18,-1.32093066782102e-15,-4.039448419784682e-15,.5264000000000041,1.822053268979825e-15,-5.772278478001558e-16,2.423818041676767e-15,4.254308611274434e-15,.7234924181056674,2.021392719935926e-15,3.193285033724458e-16,1.052351450390311e-15,-1.026521639319217e-15,-.1315999999999996,-1.397682535698542e-15,-3.578286334366862e-16],[2.349999039864425e-15,2.777409850398821e-15,-.02960000000000022,-2.072551138973959e-15,8.202938106305229e-16,-2.670337993951597e-16,-6.959247368267879e-15,0.0854636416237079,2.257592410212913e-15,-8.724355871591745e-16,-2.401860424890855e-15,-5.971377846905007e-15,-.1775999999999971,1.986425327552778e-15,-4.967681478731725e-16,6.447869634503143e-15,4.381108760045625e-15,.8553363583762875,1.889307410233463e-15,1.10843380496841e-16,1.640049340760161e-15,-1.258770221922788e-15,.2663999999999987,-5.537252738719212e-16,-5.160446293012055e-16],[.01502495999999925,-.04653806018334852,0.140232959999999,.1927383801833493,-.03505823999999979,.04824097061242492,-.1494207767384594,0.450249059049285,0.618829370215365,-.1125622647623201,-.01001663999999889,.03102537345556858,-.09348864000000319,-.1284922534555711,.02337216000000006,.004820149387577509,-.01492985021536361,.04498806095072241,0.0618322967384565,-.01124701523768009,-.001669440000000118,.005170895575928493,-.01558144000000134,-.02141537557592877,0.00389535999999987],[-.007422240000001753,.02298952222270992,-.06927424000000094,-.09521160222270833,.01731856000000042,0.0408049723811602,-.1263886400000075,.3808464088908323,.5234412790563947,-.09521160222270665,.02968896000000095,-.09195808889082391,.2770969599999956,.3808464088908254,-.06927424000000033,-.009852652381159668,.03051744094360627,-.09195808889082455,-0.126388639999998,.02298952222270656,0.00318096000000012,-.009852652381160142,.02968895999999922,.04080497238115807,-.007422239999999723],[-3.503097698274765e-16,4.151591125874708e-15,-1.767176786142608e-16,-1.264317347038826e-15,4.951005494407048e-16,-1.652990852891285e-15,-8.345262337175288e-15,2.910033769850537e-15,-6.789285793789714e-16,-7.391116998528605e-17,.05640000000000102,-.1746924181056717,.5263999999999984,.7234924181056767,-0.131600000000001,6.461345801641522e-16,2.213729537907645e-15,-1.381388187711739e-15,-6.612823256135264e-16,1.413040786869399e-15,1.055698461890639e-15,-1.470032283827525e-15,-1.219801954220413e-16,-5.89782289306416e-16,-6.155529087125914e-16],[.003180960000003812,-.009852652381154516,.02968895999999896,.04080497238115922,-0.00742223999999981,-.009852652381162385,.03051744094359589,-.09195808889082557,-.1263886400000013,.02298952222270669,0.029688960000001,-.09195808889082892,.2770969600000009,.3808464088908303,-0.069274240000001,.04080497238116294,-.1263886399999946,.3808464088908228,.5234412790563914,-0.0952116022227033,-0.0074222399999977,.02298952222270433,-.06927424000000086,-0.0952116022227051,.01731855999999762],[-.001669439999989926,.005170895575926957,-.01558143999999929,-.02141537557592742,0.00389535999999953,.004820149387571407,-.01492985021537551,.04498806095072227,.06183229673845154,-.01124701523767895,-.01001664000000012,.03102537345555745,-.09348864000000058,-.1284922534555627,.02337215999999888,.04824097061243237,-.1494207767384496,.4502490590492768,.6188293702153645,-.1125622647623187,.01502496000000311,-.04653806018335403,.1402329599999991,.1927383801833507,-.03505824000000115],[-.007885440000000637,.02276751412855924,-.04731263999999968,.2278616058714425,.07096895999999916,-.02531795620793675,.07310015999999506,-.1519077372476272,.7316002859604187,.2278616058714449,.005256960000002872,-.01517834275237174,0.0315417600000014,-.1519077372476343,-.04731264000000266,-0.00252972379206149,.007304034039584581,-.01517834275237174,.07310016000000256,.02276751412855796,8.761599999997514e-4,-.002529723792061185,.005256959999999673,-.02531795620793946,-.007885440000000746],[.003895359999997417,-.01124701523767721,.02337216000000099,-.1125622647623219,-.03505824000000041,-.02141537557592876,.06183229673844507,-.1284922534555638,.6188293702153655,.1927383801833524,-.01558143999999894,.04498806095072414,-.09348864000000429,.4502490590492791,.1402329599999974,.005170895575928298,-0.0149298502153633,.03102537345556814,-.1494207767384532,-.04653806018334884,-.001669439999999692,0.0048201493875766,-.01001664000000001,.04824097061242102,.01502495999999932],[-2.113114283615757e-15,5.074237812982654e-15,-5.517312762637234e-16,-8.82853413372274e-16,5.088589327649176e-16,-2.996737653120149e-15,-1.640521453460988e-14,6.877431893526599e-15,-1.571219076740223e-15,4.91315530770302e-16,-0.0296000000000001,.08546364162371432,-.1776000000000073,.8553363583762965,.2663999999999992,1.430600821554417e-15,1.618399192990128e-15,6.494984605077245e-16,-2.006794241206969e-15,2.106085330896492e-15,1.174880459442219e-15,-1.779553517627956e-15,9.694532485844193e-17,-2.750220458336811e-16,-1.399298255129971e-15],[-.001669439999997961,.004820149387586325,-.01001664000000303,.04824097061242302,.01502496000000095,.005170895575922909,-.01492985021538886,.03102537345557668,-.1494207767384573,-.04653806018334992,-.01558143999999854,.04498806095071781,-.09348864000000333,.4502490590492845,.1402329600000006,-.02141537557592282,.06183229673846021,-.1284922534555699,.6188293702153617,.1927383801833537,0.00389536000000206,-.01124701523768254,.02337215999999993,-.1125622647623187,-.03505824000000323],[8.761600000080947e-4,-.002529723792058178,0.00525696000000159,-0.0253179562079375,-.007885439999999208,-.002529723792072462,.007304034039555922,-.01517834275235924,.07310015999999653,.02276751412855941,.005256960000002532,-.01517834275239237,.03154176000000135,-.1519077372476262,-.04731263999999853,-.02531795620792556,.07310016000001204,-.1519077372476334,.7316002859604187,.2278616058714447,-.007885439999996385,.02276751412855125,-.04731263999999923,.2278616058714434,.07096895999999683])
