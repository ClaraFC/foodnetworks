food_web_filename = 'East_Kleinemonde_November_2006_C_Ortega_2016';
crc = 'crc_db26f6b1';
n = 44;
l = 42;
b0 = [0.182 0.375 0.0429 0.39 0.00826 0.00976 0.0455 0.272 0.00228 0.0177 0.213 0.0805 0.000368 0.000634 0.000751 0.0794 0.021 0.0104 0.00272 0.002 0.00336 0.00446 0.862 0.0198 0.0202 1.11 0.174 0.00106 1.35 0.00228 0.163 0.925 0.00568 0.0805 0.00143 0.615 0.109 0.0104 0.00958 0.0279 0.00219 0.00189 0.622 71.7].';
p = [21.15666 43.4826 4.96944 9.64656 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 37.6362 0].';
q = [0.700070126936428 1.55006587729102 1.0406592 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 12.7811428029837 0].';
r = [17.4469460131599 38.6968355504761 2.6031096 3.43603620144191 1.39000371865208 0.332621712476908 1.60702309931902 4.51103938477474 0.250408894645117 0.697875118193878 4.44840982893279 1.89120758500171 0.0269649153221293 0.0449549367634356 0.0933542334360551 1.60015965226287 0.590618373480182 0.345639071242573 0.141423303780767 0.111344255254177 0.399829298072175 0.177844286861839 3.76724586308042 0.211508559650078 0.102427628155701 2.9604126139924 1.64570791991914 0.00274690316936628 1.00988289306549 0.0197940276693609 0.00531336974161492 0.123703746125534 0.00338546890574272 0.155287786196266 0.0072188922797846 0.066587496617889 6.5028264150363 0.537738114243925 0.276276383794995 2.32190455320162 0.183489039527994 0.0724152848654044 0 0].';
F = [
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0.46951467281488 0.504630121767776 0.16457616 0.211541141871181 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.588999491092331 0
0 0 0.15792588 0.0573708054478761 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.153270570868089
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.94019063132772 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5.67046307478602
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.924056893903539 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.00420386238986305 0 0 0 0 0 0.0201888212298966 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2.33188260059011 2.23639848458539
0 0 0 0 0.0461751563830809 0 0 0 0 0 0 0 0 0 0.00419655767676163 0 0 0 0 0 0.0201778569895473 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5.3160591326119
0 0 0 0.365419142980103 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.0179606589226 0.976254982978625
0 0 0 0.00420597433570098 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0117178744909002 0.0111472827584872
0 0 0 0.00755321368539873 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0210351741328694 0.0201754491617065
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.42384578763468 0
0 0 0 0.542282008182473 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.44867178424785
0 0 0 0.187935065234667 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.502011357265206
0 0 0 0.107360144207554 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.286701828272203
0 0 0 0.0244136529425007 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0679980181713913 0.0651080231911643
0 0 0 0.019082187646421 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.053126877174209 0.0508953290153153
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.629191764057688 0.603582796362533
0 0 0 0.0544109103897373 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.14519600682474
0 0 0 0.736684992247887 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 24.9155553955248 0
0 0 0 0.119711311240282 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.333340096602688 0
0 0 0 0.0795882893410664 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.221605577758994 0
0 0 0 0.736684992247887 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 25.7851700287753 0
0 0 0.16457616 0.736684992247887 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7.20903223278641 0
0 0 0 0.000861292920004102 0.000660845394754694 0.000150215078544718 0.000588291766690227 0.00107863152252954 0.000287340097050108 0.000747566749646427 0.000946703488655377 0.000602146464467315 7.2527583894887e-07 1.61691577432064e-05 0.00013334753766629 0.000486735020966921 0.000199918471716373 0.000118196870557053 5.16226112741471e-05 4.03061088104225e-05 0.000289053923076127 6.21194344904717e-05 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0.736684992247887 0.0461751563830809 0.0015657300208322 0.0337080220355115 0.20615632621594 0.0191304359045438 0.422378867346661 0.162145817593639 0.0610550207600051 7.2527583894887e-07 2.19913702850075e-05 0.00419655767676163 0.0495505021344703 0.0064142032515948 0.00200002989004423 0.000274553378432566 0.00015958003370364 0.0201778569895473 0.000511571813450943 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0.00867870464577744 0.00665959731695858 0 0 0.0108264945234594 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0.0461751563830809 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0.736684992247887 0.0461386254049297 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0.00208409230352434 0.000473555684247031 0.0018551287574449 0.00340172548910179 0.000895335085011205 0.00232856446505214 0.00306709998969392 0.00195039916706707 7.2527583894887e-07 2.19913702850075e-05 0 0.00157677298458879 0.000647633513494358 0.00038290670894526 0.000167216342393327 0.00013056548212116 0.000936711600510954 0.00020123042690531 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0.0461386254049297 0.0015657300208322 0.0337080220355115 0.129559973261506 0.0191194770907737 0.0897737782933172 0.118246555436299 0.0610550207600051 7.2527583894887e-07 2.19913702850075e-05 0 0.0495505021344703 0.0064142032515948 0.00199528009308212 0.000274553378432566 0.00015958003370364 0.0201778569895473 0.000511571813450943 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0.00086432294305672 0.000196437913465435 0.000769164080051509 0.00141029336555589 0.000375887312311926 0.000977390545114463 0.00128737053553832 0.000818816885237411 7.2527583894887e-07 2.19913702850075e-05 0 0.000661769611840161 0.000271918360067608 0.000160762358717476 7.01819081795022e-05 5.48133846016625e-05 0.000381555564156627 8.19611126836047e-05 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0.0461386254049297 0.0015657300208322 0.0337080220355115 0.206119799581246 0.0191194770907737 0.422378867346661 0.162145817593639 0.0610550207600051 3.62601381648795e-07 2.19913702850075e-05 0 0.0495505021344703 0.0064142032515948 0.00199528009308212 0.000274553378432566 0.00015958003370364 0.0201778569895473 0.000509379362821868 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3.34740439880501 0.00822797160704742 0.00398674785293863 3.48621287954844 1.07619434674402 2.44209778841212e-05 0.119585626230223 6.83182571829178e-05 0.016484601722438 0.187620393680499 0.000159847037153296 0.0712496901371103 1.58625659305793e-05 0.265399780702036 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.275290303903471 0.00822797160704742 0.00398674785293863 0.286687205250362 0.088477349212967 2.44209778841212e-05 0.00982845507361136 6.83182571829178e-05 0.00137475219862141 0.0154182010311263 0.000159847037153296 0.0197233757620575 1.58625659305793e-05 0.02181453168563 0 0 0 0 0 0 0 0
0 0 0 0 0.0417914390049403 0.0015657300208322 0.0336860981187404 0.0681952269746918 0.0181697132307047 0.0472801258085277 0.0622119435402105 0.0395706687510984 3.62601381648795e-07 2.19913702850075e-05 0.00419655767676163 0.0319849222111371 0.00640689361825965 0.00199528009308212 0.000274553378432566 0.000157972174421991 0.018997373778603 0.000509379362821868 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.2109265557952 0.00822797160704742 0.00398674785293863 1.26097061600688 0.389417234768881 2.42126576597281e-05 0.0432598171269734 6.83182571829178e-05 0.00605154077862057 0.0678634731329737 0.000159262453703135 0.0712496901371103 1.58625659305793e-05 0.0960073290972526 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0932859232632813 0.00337741427777737 0.00234125513233527 0.0971579503900986 0.0299953556724885 2.42126576597281e-05 0.00333217510302363 6.83182571829178e-05 0.00046629021941545 0.00522954389086081 0.000159262453703135 0.00667554789130772 1.58625659305793e-05 0.00738236790165399 0 0 0 0 0 0 0 0
0 0 0 0 0.0101665712194711 0.00155988365436001 0.00904726965419257 0.016553870843562 0.00440909607347419 0.0114729207912501 0.0151070483344119 0.00960949758221504 3.62601381648795e-07 2.19913702850075e-05 0.00212822816209661 0.00776875866798554 0.00319065495079331 0.00188640013195074 0.000274553378432566 0.000157972174421991 0.00461229044028413 0.000509379362821868 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
2.54012918708875 2.73106845046514 0.838593000000001 0 0.210089655347388 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 16.9580873429255 0.213481519094041 0.184464740253293 18.4304137558352 4.88050117871693 0 0 0.00609749578814384 0.016484601722438 0.382988259791903 0 0 0 0 2.07980847083155 0.193359228172085 0 0.836324079036744 0.0660224401487626 0 0 0
0 0 0 0.736684992247887 0 0.0273025314251066 0.186097513525051 0.516121348233741 0.592141237373791 2.89746056925521 0.413040518216438 0.232710608749535 0.000101502079627202 0.00363879146680296 0.311436443078733 0.199663675267512 0.0693684203505808 0.0378887649977323 0.014434602770287 0.0107397691562857 0.706828027853369 0.0188660376631944 0 0 0 0 0 0.00447705745406293 0.586418973722908 0 0 0 0.0160979667587982 0.254087356519726 0.00112301998320131 0.574143563044476 0 0 0.100739847919958 0 0 0.0260614645279985 8.94826149928733 0
];

