food_web_filename = 'Durban_Eddy_Summer_C_Scharler_2016';
crc = 'crc_058fcd49';
n = 32;
l = 28;
b0 = [1.1495 1.2276 1.84 0.1964 0.3662 0.1714 0.0727 0.0346 3515.2 2122.8 154.39 138.3 335.05 88.01 720.6 233.29 130.33 336.31 73.615 52.748 83.209 13.953 252.27 83.127 5.637 15.634 238.45 18 7.28 378.6 271.2 30000].';
p = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7.62 8800 173000 711 1820].';
q = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 301.756807737139 0].';
r = [170.586914714947 191.164305827543 98.9446267083655 14.345685448767 33.5999348427858 15.7524248973407 6.68412034729969 0.252466124372632 175655.727582863 319.238299973808 309.444092440263 562.484309681951 799.527569104001 564.71725940756 1228.21337276143 273.523475337213 558.518000823069 751.894079355541 208.867723909056 196.205087614759 235.214620133619 44.6191665980445 670.225047703319 576.772993436904 35.6654841269635 18.1806735704719 473.106054601095 23.3878199091082 0 0 0 0].';
F = [
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 858.565086638974
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 961.434913361012
5.90757369079949 5.91186674369733 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 15.8174555412004 3.9668111946103 369.779091501454 0
0.630479679538814 0.630937851160013 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.6881015494173 0.423353813540058 39.4641007613669 0
28.094325843337 20.0489393916692 1.8393838572559 0.0118710318799472 1.20522589336353 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 9.0606511632746 0 0 0
6.08069924941595 4.33945365276752 2.84059958362809 0.0170455842378729 3.38412592397691 0.28252065962506 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 11.3073216052869 0 0 0
2.58015819589569 1.84130638956882 1.20535496851454 0.00720378857672009 1.43597988287563 0.119881468016291 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4.79791486679777 0 0 0
0.040220083772885 0.0286925787358414 0.0273944311026032 0.000710232676578037 0.0827682849157482 0.039660785668723 0.0294709896192238 0.00829560107796058 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.434420258986711 0 0 0
0 0 0 18.6841924845494 0 0 0 0 12415.4939305612 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8252.18566518796 226725.921534617 0 0
0 0 215.381104980023 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 848.758826756111 0 0 0
7.17629831357585 5.1212267963097 0 0.121145402262025 0 0 0 0 186.437946893955 63.9679721583145 1.33857268693717 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 76.1195056191832 278.609325182627 0 0
0 0 0 0.0855323066221836 0 0 0 0 214.792155757419 113.694994557425 1.22923804581652 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 340.049304888838 466.200588652722 0 0
0 0 0 0 0 0 0 0 1847.63564687379 0.321634517231235 0.579643963558028 17.8552454777758 178.711085647585 4.09647099068521 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 610.397537343587 211.077664939813 22.8380122502142 21.7420471190331 10.8806816732921 5.44003808107343 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 2044.85180258253 50.5511607194479 0 0 123.75622363937 12.6376115023831 229.69094005381 41.9299111744263 0 0 0 25.2750651955445 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 621.326139769746 0 0 0 73.4105114332317 18.2470868181969 0 125.773699906846 0 0 0 0 0 0 0 72.9908536028585 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 679.167522702321 0 51.0044597601761 71.5562185661468 44.7619934226368 15.1119966819056 0 0 8.950823615811 0 0 0 0 0 8.95167299174021 8.90132117861376 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 723.53363805622 0 0 0.714350064531076 52.4492571886827 11.8024563749348 303.343637493038 4.57328845721355 23.6028113243815 2.07418019251004 0 46.0776168711635 0 0 0 11.9330280894631 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 179.379473400363 13.4579180947471 0 0 89.7255180852598 4.57570766645983 17.9449493434156 13.4512020862847 26.9129189053162 44.8591970196055 0 0 0 0 67.2881806478545 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 151.336078000396 11.3544604352984 0 0 75.6985427250168 3.86046196667402 15.1396097377124 11.3487941315037 22.7056137820022 37.8463445634173 0 0 0 0 56.7694209989236 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 206.099405020282 1.14908142206583 14.356810895457 115.041314557213 2.29825825023224 5.17074316242685 0 0 0 0 0 0 0 0 231.408702074553 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0.669341120165422 0.0669644906293204 1.6793953269555 26.7885058164008 0.0904942080036134 7.65128924117502 17.2554857719644 7.57935862242885 4.58549702883626 0.0631323647050538 0 0 0 0 0 0 0.709103955364922 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 880.995749159719 0 47.0058784762955 17.2111182985674 86.134050912966 15.6754506814565 155.035619453143 103.306594880063 0 51.675549249716 40.1934769428985 0 78.3631064519311 7.83684951089022 84.6443345786865 0 0 0 0 0 0 0 0 0
391.867063729777 349.511474725968 0 0 6.49571483268304 5.57778500312799 1.39442736645805 0 93.3001859971171 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 79.9949021771479 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 16.0389627135222 0.448039369200469 0 0 3.56768490150422 1.94711242031999 0 0.594568013419762 23.1889926305537 4.44403112925718 0.494919718203733 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 5.03808571845892 13.107007617718 0.900733047105733 0.613874238598575 3.68651772727244 0 1.2288149466583 7.36905011262037 5.11913851387074 0.614363972144692 0 0 0 0 2.25263067624691 0 0 0.819112724659937 0 0 0 0 0 0
0 0 0 0 0 0 0 0 298.671685938373 0 37.0384876368278 0 39.072997125709 14.8222489846276 11.6053630936438 122.957796020795 37.0497474511708 0 22.2332388297759 35.5725739073033 0 0 62.7426167660157 58.9523488780616 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0.0566461761114355 0.38881726200562 1.86320643962074 0.394466182040523 0.10967920782087 0.233100109484258 0.216958512824904 0.178227826378681 0.205707901456306 5.92657716971545 0.231781873918142 1.61823416263861 9.35007501428919 3.12297693431985 0.0399714732448136 0 0 0 0
245.601353137889 382.836709403534 81.1443341430295 9.56358737555609 14.0566475201492 6.4794934450572 3.2104397369606 0.363907030419933 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 96.9575978216372 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 50596.421360885 239.431231483763 133.065569642522 321.182046963895 451.831034763395 198.582444290124 545.477130980292 205.140101250342 205.15248832171 291.030870551375 182.499840556643 80.7649118190194 261.768360970352 13.8822759228408 353.970498403355 96.9575978216372 12.9456728881839 12.3994679859488 264.490073096876 8.12863486355243 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
];
