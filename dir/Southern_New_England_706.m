food_web_filename = 'Southern_New_England_706';
crc = 'crc_40eab51b';
n = 33;
l = 31;
b0 = [7.5316 5.08305359 11.825181 5.652736 1.19584119 4.232795 0.228357792 35.4362373 6.39222431 17.8053951 18.9326153 3.70247149 3.37289071 0.269999981 0.422238052 14.8511448 1.945768 3.06976056 0.659100056 0.349553078 2.3374362 3.63476229 2.33421755 0.0214726068 0.03875167 0.00547841657 0.0178088564 0.14004381 0.07523035 0.01070364 26.5280533 0.974253953 40].';
p = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5821.08195383212 35.4385769897232 0].';
q = [0 0 0 0 6.35856762827351e-07 0.00961714731419753 0.099927710053883 0.000354340703646022 6.14034442968885e-05 0.000431866732710461 1.66859287290282 0.672310690749345 0.0523976205078634 0.0403582921748046 0.440481314311901 0.195238402514672 0.0351548692511008 0.248356204891113 0.104237892643025 0.102663080167094 0.0399639921007455 0.301180205185157 0.0994959731324433 0.00476592150450134 0.00309969463282205 0.00227873565232052 0.032051565568209 1.70705723414356e-08 1.30602316644622e-08 6.32691479874417e-05 0 78.6295642390996 765.161155064204].';
r = [1603.58162843629 743.041687467187 595.75547903885 113.838102607544 65.6514590145689 381.858864291522 0.137213273896076 221.46317576458 47.9378110711117 64.0959567606758 128.350504113484 2.59171531557894 33.2380490606725 0.404812234728254 7.28353513442902 19.8542160903984 0.680284916024986 4.25801707026818 0.817954334171028 0.219777966411245 0.420637963311627 0.254354119058245 0.96103714076408 0.0137909536403125 0.00935954894173757 0.030744353915144 0.123168524477973 0.497356189100192 0.938227255484526 0.0875098210993969 970.180297955093 0 0].';
F = [
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 429.536387743401 0 2434.00291024341
263.530104510836 159.211408375872 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 242.852593070327 0 566.647465328088
0 124.682996928022 139.472997103815 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1054.66896990557 0 191.837576745884
0 36.6757672834099 116.864962411732 75.552747136479 36.6491510551131 0 0 0 0.743456097703079 0 0.396104838106831 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 337.420790371872 0 14.6702464522517
3.51300973615239 8.78254116906305 64.8156927666583 72.0174947571635 3.51304263985323 0 0.151688060723427 0 0 0 0 0 0 0 1.75650990993654 2.40828533100626 0.227609881369521 0.526943658300077 0.0167608605787457 0 0 0 0 0 0 0 0 0 0 0 15.1061395689541 0 1.75650023414616
0 0 147.470374001762 141.863294195692 0 33.7473261758479 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 132.662596214541 0 133.824260762809
0.00847783576618199 0.0254335537933198 0.21194798774503 0.156841661907922 0 0.00423895242228719 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.00133221818480693 0.00847794475808622 0 0
182.900943388134 0 0 0 0 0 0 29.5239555704688 2.71055323870293 1.4598123510232 5.067793518763 2.1104120425991 0.299255979143055 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 75.9752878320844 3.49013085445994 316.558776222712
19.8077247846993 0 2.05553706216216 3.92421054086688 0 0 0 18.685636968028 2.50242303631914 5.82376023858841 6.76681739839354 2.05553135354805 0.331970279938444 0 0 0 0 0 0 0 0.0858942887764519 0.250319378253598 0 0 0 0 0 0 0 0 22.4241093601396 0.927101099913281 48.5848396558249
48.7978312577637 0 0 0 0 0 0 0 0 4.16430564350806 2.07393400250421 3.06731462087694 0.438557684785387 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 105.683496984376 1.383442772745 83.6532452774052
74.5303697034025 0 0 0 0 0 0 23.2896302618243 3.18414801617097 11.4176648608979 8.21258618884762 1.86327244286939 2.75977823118388 0 0 0 0 0 0 0 0.276034351929138 0.450546189810617 0.225701472941654 0 0 0 0 0 0 0 69.407318019498 3.2354876130339 135.085941537099
14.0691687238986 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 42.2082087956106 0 14.0691856113205
10.9073594626563 0 0 0 0 0 0 6.08965300364498 4.11617479485301 8.27235331797534 8.95151650575644 2.99954500485336 1.52952805011186 0 0 0 0 0 0 0 0.0750614768711661 0.139306624034093 0.151154127131208 0 0 0 0 0 0 0 0 2.97632862178914 10.9073337921645
0.459732776511054 0 0 0 0 0.155161046781011 0 0 0.0775746491634193 0 0.0409889982164518 0 0 0.00162412091182926 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0775809077862688 0.0769800000591577 0.459731694529063
0 0 8.24589218005945 4.77393970916797 0 1.30198067962132 0 0.867938735711909 0.504740163809564 0.245447776234325 0.198096714696382 0 0 0 0.605731106504612 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.12838796780203 0 1.12837080569084
0 0 3.78937405362056 12.9921519937699 1.8469205018698 6.7826480740857 0 0.159208331364655 1.05991045754204 0.14323228320509 0.142770256749436 0 0 0 1.58018362551319 0 0.766172517587458 0 0 0 0.025533049833959 0.0318336450921783 0.0318410368166481 0 0 0 0 0 0 0 0.350279290004961 0 0
0 0 0.448943503135171 1.75673721313031 0.866655550794254 0.16396160405929 0 0 0.00545344186770544 0.00382038817243556 0.00571999536650765 0 0 0 0.0273268041222989 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.61290718251101 0 0
0 0 0 0.874160915635367 0 4.25423934807782 0 0 1.31113818164338 0 0.145679470591661 0 0 0.131400770813383 0.759868206818894 0.135979845886163 0.145693404013763 0.557692492239808 0.00920495985815608 0 0 0 0.116545691350169 0 0 0 0 0 0 0 0 0 0
0 0 0.104501898124964 1.15392251257379 0 0.0251333279919296 0 0.00396820050580004 0.00553691013121146 0 0 0 0 0 0.0105824070378665 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0145509395847965 0 0
0 0 0 0 0.00539006003583135 0 0 0 0.0118572230361396 0.000538978961922713 0.0118570137135597 0.000539004946047457 0.0215503734784272 0.0212338041033223 0.00646803306890293 0.134750955412504 0.0970211763020602 0.0431194155145121 0.0192064333677034 0.0231910040147214 0.0484991329484823 0.0215535933759268 0.0323378971317946 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0.00229090087897394 0 0.199297421956266 0.373388664827791 0.0893407566696469 0.419195666833381 0.142035795238223 0.171739134227405 0.119384642568523 0 0.114544712051016 0.0183272501875951 0.0687254243051518 0 0 0.105679151924349 0.0847373144765312 0.139734539227947 0 0 0 0 0 0 0 0 0.0272780163654474 0.0274905310945408
0 0 0 0 0.036936560043905 0.107731477024705 0 0.606340958091006 0.544772232931365 0.369347101023559 0.363175177153186 0.13235575008447 0.307663798534812 0.0817612881312053 0.0369363306258621 0.107731164869537 0.0109631502027925 0.0246237664561135 0 0 0.0417076355224496 0.0153854883692632 0.00923343551398109 0 0 0 0 0 0 0 0 0.0733011969899577 0.0369361271609171
0 0 0 0 0.0364454817716949 0.0026032453627184 0 0.00260310155008758 0.00260304791516192 0.00780935226266035 0.101517114514871 0.494616097550597 0.145715145535014 0.0768740143098347 0.0286355613886704 0.627380336556526 0.0520650184902943 0.408699988707482 0.099396436918594 0.0580731209008357 0.122323875291048 0.0546513439555867 0.0937097952816901 0 0 0 0 0 0 0 0 0 0
0 0 0 0 1.93388019591434e-05 0.000657518408120587 0 0 0.000212710432725992 0 0.000212706677622759 0 0.000212629461398044 0 0 0.00328758351469921 0.00328760006449337 0.00233993709419841 0.00218528555644344 2.69231073985891e-05 0.00110205270577801 0.00110197208741698 0.00110222796367683 0 0.000212364978811962 0.000212592854738716 0.00058178250638623 0.000212724763035123 0.000290072598093684 0.000448191395111123 0 0.00109377791192718 0
0 0 0 0.000512649794918401 1.70882998772187e-05 0 0.00017088252921432 0 0.000170869909559542 0 0.000170866893091842 0 0 0 0 0.00340055863891764 0.000854415977746842 0.00256318280634237 0.000341766292264224 0.00137134096335386 0.000854214217503189 0.00136664248688493 0.00119609046676156 0.00024311773844908 0.00017505795412266 0.000170882699055685 0.000512625705324557 0.000170881421071196 0.000341755044778437 0.000349367465765855 0 0 0.000854405280424819
0 0 0 0 0.00373768684947416 0 5.19121443287943e-05 0 0 0 0 0 0 0 0 0.00742342945368013 0.0254490308310077 0.00327039488726881 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0.0126558285631214 0 0 0 0 0 0 0 0 0 0.0761057700291498 0.0253117112042345 0 0.00860092063923263 1.95384030341083e-05 0.00831390870287268 0.0126519414772349 0.0126548792390004 0 0 0 0 0 0 0 0 0 0
0 0 0.0370748329017872 0.283175289579778 0.0127844150198326 0.184734603218146 0.00319609410225993 0 0.0370719480250338 0.00767026375620478 0.00725884600119155 0 0.00447249772973585 0 0 0.0338780173780896 0.00229804675404198 0.00703125596024869 0.000471596109446625 0 0 0 0 0 0 0 0 0 0 0 0 0 0.00767056011468061
0 0 0 0 0.0609712354812776 0.0304855820097682 0.0243884029051366 0 0 0 0 0 0 0 0 0.371922916858822 0.231690870277629 0.285338363265337 0.0162555520687038 0.00014706884703361 0 0.0768000416064387 0.0768178744604891 0 0 0 0 0 0.00180271935508707 0 0 0 0
0 0 0 0.00341456363947536 0 0.0139996720395225 0.000113818220185101 0 0 0 0 0 0 0.00110212995367046 0 0.0299341679507545 0.0265197579208334 0.00705658224691788 0.00326908936072785 5.57356142581473e-05 0 0.00420999011139046 0.00284524874643142 0 0 0 0 0 0 0 0 0.0138913037502163 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
6.41432935370859 1.59821701508174 4.31389753842754 1.8978249089979 0.659086881095444 1.60909518875435 0 3.19205115342245 0.690968133862283 1.53160548309019 1.7101397862413 0.542169128119492 0.178144277949164 0.00470822884268656 0.0646426710554409 0.0559797960454147 0.0154282180669966 0.0199782572998039 0.00220311150388166 0.000937883181830703 0.00851584845693112 0.0120690433929006 0.00460314634095648 0 2.96199619821715e-05 6.52588877464832e-05 0 0.00131048452296296 0.00235958838594994 0.000167091958577957 23.207639337782 0.703281646373177 7.72818423659429
635.018618016996 158.223519492698 427.075866304298 187.88469499073 65.2496004284518 159.300426686666 0 316.013107186239 68.4058491520513 151.629030821188 169.303846837102 53.6747458838172 17.6362800185819 0.466114667420408 6.39962473448573 5.54200003849428 1.52739384863195 1.97784729268566 0.218107956884526 0.0928504378009721 0.843068980240242 1.19483510595644 0.455711501753574 0 0.0029323761162468 0.00646063015690026 0 0.129737969773306 0.233599223210012 0.0165421003029722 2297.55593444043 0.0703281646373177 0
];

