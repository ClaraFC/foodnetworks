food_web_filename = 'Richards_Bay_Summer_C_Scharler_2016';
crc = 'crc_cd9d7b71';
n = 33;
l = 29;
b0 = [0.5011 0.3424 0.52 0.0818 0.0375 0.0484 0.0251 0.0263 1482.3 435.9 1.006 1.218 0.585 8.078 385.78 39.644 34.767 5.5544 44.056 6.5719 22.552 3.569 345.6 48.384 9.1443 17.13 7.5971 18 6.6431 0.3321 378.55 271.2 1364.2].';
p = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 26.9 0 0 5.51 0 45800 51200 7700 3320].';
q = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 28704.8621088078 2679.24210018978].';
r = [74.3109614784351 53.0978812336593 27.4820836469452 5.8785987734423 3.43487614942203 4.43999151647495 2.30397001388113 0.191492287296119 73865.4855855142 65.3636133944975 2.01097475899592 4.94964384317596 2.48954090541664 51.7349746767881 654.287429073973 46.2463902191682 148.89223855899 12.4017577291591 124.818054931244 24.4126940069557 63.7524489097862 11.3811061798753 918.55123347456 335.498003105129 57.6760001617661 19.8268601511529 15.0343137632756 23.3091008374991 49.0439717075911 0 0 0 0].';
F = [
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 372.833689791247
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 267.924210018978
1.66940591835571 1.66632633760333 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4.47698676876691 1.11858470708761 102.520212378888 0
0.262804509508592 0.262319709676373 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.704785036943628 0.176092046911824 16.1464849362044 0
2.87691300181754 2.04688864368685 0.185474149431488 0.00119727062595566 0.123456418124154 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.929034821425692 0 0 0
1.79184892846767 1.22217137462856 0.788763721507026 0.00478908250382265 0.874151089621025 0.0796701675723923 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3.19355719865081 0 0 0
0.889951634472276 0.633939298384569 0.408840867026398 0.0024943138040743 0.206092568965321 0.330047878178496 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.66185108142958 0 0 0
0.030560978946643 0.0217606122799719 0.0204420433513199 0.000498862760814859 0.0576461823337783 0.0350987471658099 0.0223492674133628 0.00628334067690392 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.330367986067326 0 0 0
0 0 0 7.74235004784662 0 0 0 0 5219.18641466741 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3483.88058034634 95349.8413107773 0 0
0 0 58.6337633686638 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 159.177302377893 0 0 0
0.0467871664655447 0.0332867813415095 0 0.000798180417303775 0 0 0 0 1.20981141029647 0.415870291999739 0.00869556680082728 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.496553094210283 1.8109466188091 0 0
0 0 0 0.000698407865140803 0 0 0 0 1.89970386740768 0.998688513720528 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3.00334532788478 4.09213904471227 0 0
0.177193949592914 0.176867076978767 0 0 0.11847833674818 0.118657696384414 0.11873048313349 0 0.71388871648899 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2.14111920676877 0 0
0 0 0 0 0 0 0 0 46.7927057866733 15.1952606692212 0.403794136498186 0.646593582069457 0.00809168224945821 0.391272076802726 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 17.4091000924189 0 0
0 0 0 0 0 0 0 0 1099.82855481497 55.9825393076572 0 0 0.134861370824303 1.35093939561044 122.712564081758 22.3340262792266 0 0 44.6493844235924 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 118.981452748165 0 0 0 0.0154840833168645 1.5510785653305 0 21.2372660601574 0 0 0 0 0 0 0 12.3538982583626 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 227.964464088921 0 0.235879743102901 1.18108425207734 0.0235757655663227 0.236164220269676 0 0 2.38419920453741 0 0 0 0 0 2.38701963822619 2.37115143991152 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0.0187502864644117 8.86861752837163 0 0 0.0118108425207734 0.00189804892271242 0.83658172942987 0.150647131515004 0.0751779277434681 8.50497951534565 0.0342435099577943 0 0.759955489757994 0 0 0 0.197263859286757 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 196.969295725954 8.02749626143728 0 0 0.0272719661000258 1.37095331258244 10.6749954119903 8.00634959920489 7.88388560491993 0.273347316329762 0 0 0 0 40.3185669986104 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 30.7951995348192 1.40956036471066 0 0 0.0047950709626419 1.30090460318042 1.87560667051793 1.40584718989774 2.82497552806534 1.34170477612411 0 0 0 0 7.07079346617421 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0.0155587483428097 22.696461994818 0.30990334259596 0.155920508152765 0 0.0155839806285862 1.40097418804045 0 0 0 0 0 0 0 0 131.386375045223 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 7.82877962200111 6.82787041913034 0.0230882290918517 0 0.00169825429926901 1.08075151648835 1.16726585346063 0.0160525813881941 0 0 0 0 0 0 0.181533846436529 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 1839.71321896322 0 0 0 0 1.07074455800235 251.411106899212 0 10.7188787767018 0 42.9475020576707 0 0 0 0 0 0 0 0 0 0 0 0 0 0
249.863378358547 186.803429618024 0 0 0.0539624021155576 1.07689337895098 0.540772452591189 0.151598060776095 54.6914744985263 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 46.6519640931436 0 0 0
0 0 0 0 0 0 0 0.00967434993110603 0.537916147718595 0 0 0 0 0.534371583152573 0.871957569166314 2.03399167900099 0 0 0 0.00500629439893276 3.79932356387335 0.0488832755551956 7.66253362859162 37.9583486809366 1.64825675506026 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 10.198410235557 14.295541287491 0.222886367423504 0.672617472369466 0.0444543037161593 0 1.33686858430533 8.03626124154314 5.6098804812645 0.672854932504029 0 0 0 0 2.46725559245228 0 0 0.893011918639315 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 9.60850219251988 0 0 0 0.0470516338209237 1.18082110134838 0.368137692245274 3.89848405141857 0 0 0 3.03381440575326 0 0 1.995869361374 3.4072848422258 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0.000500459428846329 0.0063934279501892 0.0946658272775906 0.343196114179876 1.08678967162308 0.0846490894047947 0 0.0236261316680892 0.974224890032316 0.122300125275079 0.178938779639673 2.88849435213925 6.37620555270326 2.20766510829283 5.82702749145654 2.68862545028663 0.0397370947199719 2.99298022992276 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 11.4731259100831 5.06503810261032 0.292515196523078 0 0.14616166201445 2.92367592897673 3.99981557252102 1.42951092114376 7.94335946838294 28.8921814106866 3.69609090528663 0.44201595525946 0 0 0 0 0 0 0
40.9138838666785 21.959339332757 23.9321483137403 3.92106130000479 1.29430115775323 1.87459217817394 1.14739542523961 0.131650947516082 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 56.3895678567195 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 21296.6801977808 48.9847218942001 0.961509800275384 2.53232470996242 0.744234972326712 16.7116206716254 290.319968681233 34.6975051123699 49.5873363968916 4.7360176052739 60.9674329909591 15.9200161886062 84.3068896363453 4.08860116606851 1023.00841638265 56.3895678567195 16.7822505969771 17.4611269007687 5.81702606716289 8.09718186379327 14.2665390959652 21123.5288061229 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 21123.5288061229 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
];
