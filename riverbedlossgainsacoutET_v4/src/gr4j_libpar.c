void gr4j_libpar(int * ierr,int * npar,double * config, int * nlib,double * parlib)
{
    //int version;
	*ierr=0;
	if(*npar<4) 	*ierr=20501;
	if(*nlib<1) 	*ierr=20502;
	if(*ierr>0) return;
    //version = (int)config[6];

	// Return parameter library
       if(*nlib>=1){parlib[0]=100.000000;parlib[1]=-2.000000;parlib[2]=20.000000;parlib[3]=0.500000;}
       if(*nlib>=2){parlib[4]=733.333333;parlib[5]=-2.000000;parlib[6]=20.000000;parlib[7]=0.500000;}
       if(*nlib>=3){parlib[8]=1366.666667;parlib[9]=-2.000000;parlib[10]=20.000000;parlib[11]=0.500000;}
       if(*nlib>=4){parlib[12]=2000.000000;parlib[13]=-2.000000;parlib[14]=20.000000;parlib[15]=0.500000;}
       if(*nlib>=5){parlib[16]=100.000000;parlib[17]=-1.333333;parlib[18]=20.000000;parlib[19]=0.500000;}
       if(*nlib>=6){parlib[20]=733.333333;parlib[21]=-1.333333;parlib[22]=20.000000;parlib[23]=0.500000;}
       if(*nlib>=7){parlib[24]=1366.666667;parlib[25]=-1.333333;parlib[26]=20.000000;parlib[27]=0.500000;}
       if(*nlib>=8){parlib[28]=2000.000000;parlib[29]=-1.333333;parlib[30]=20.000000;parlib[31]=0.500000;}
       if(*nlib>=9){parlib[32]=100.000000;parlib[33]=-0.666667;parlib[34]=20.000000;parlib[35]=0.500000;}
       if(*nlib>=10){parlib[36]=733.333333;parlib[37]=-0.666667;parlib[38]=20.000000;parlib[39]=0.500000;}
       if(*nlib>=11){parlib[40]=1366.666667;parlib[41]=-0.666667;parlib[42]=20.000000;parlib[43]=0.500000;}
       if(*nlib>=12){parlib[44]=2000.000000;parlib[45]=-0.666667;parlib[46]=20.000000;parlib[47]=0.500000;}
       if(*nlib>=13){parlib[48]=100.000000;parlib[49]=0.000000;parlib[50]=20.000000;parlib[51]=0.500000;}
       if(*nlib>=14){parlib[52]=733.333333;parlib[53]=0.000000;parlib[54]=20.000000;parlib[55]=0.500000;}
       if(*nlib>=15){parlib[56]=1366.666667;parlib[57]=0.000000;parlib[58]=20.000000;parlib[59]=0.500000;}
       if(*nlib>=16){parlib[60]=2000.000000;parlib[61]=0.000000;parlib[62]=20.000000;parlib[63]=0.500000;}
       if(*nlib>=17){parlib[64]=100.000000;parlib[65]=-2.000000;parlib[66]=80.000000;parlib[67]=0.500000;}
       if(*nlib>=18){parlib[68]=733.333333;parlib[69]=-2.000000;parlib[70]=80.000000;parlib[71]=0.500000;}
       if(*nlib>=19){parlib[72]=1366.666667;parlib[73]=-2.000000;parlib[74]=80.000000;parlib[75]=0.500000;}
       if(*nlib>=20){parlib[76]=2000.000000;parlib[77]=-2.000000;parlib[78]=80.000000;parlib[79]=0.500000;}
       if(*nlib>=21){parlib[80]=100.000000;parlib[81]=-1.333333;parlib[82]=80.000000;parlib[83]=0.500000;}
       if(*nlib>=22){parlib[84]=733.333333;parlib[85]=-1.333333;parlib[86]=80.000000;parlib[87]=0.500000;}
       if(*nlib>=23){parlib[88]=1366.666667;parlib[89]=-1.333333;parlib[90]=80.000000;parlib[91]=0.500000;}
       if(*nlib>=24){parlib[92]=2000.000000;parlib[93]=-1.333333;parlib[94]=80.000000;parlib[95]=0.500000;}
       if(*nlib>=25){parlib[96]=100.000000;parlib[97]=-0.666667;parlib[98]=80.000000;parlib[99]=0.500000;}
       if(*nlib>=26){parlib[100]=733.333333;parlib[101]=-0.666667;parlib[102]=80.000000;parlib[103]=0.500000;}
       if(*nlib>=27){parlib[104]=1366.666667;parlib[105]=-0.666667;parlib[106]=80.000000;parlib[107]=0.500000;}
       if(*nlib>=28){parlib[108]=2000.000000;parlib[109]=-0.666667;parlib[110]=80.000000;parlib[111]=0.500000;}
       if(*nlib>=29){parlib[112]=100.000000;parlib[113]=0.000000;parlib[114]=80.000000;parlib[115]=0.500000;}
       if(*nlib>=30){parlib[116]=733.333333;parlib[117]=0.000000;parlib[118]=80.000000;parlib[119]=0.500000;}
       if(*nlib>=31){parlib[120]=1366.666667;parlib[121]=0.000000;parlib[122]=80.000000;parlib[123]=0.500000;}
       if(*nlib>=32){parlib[124]=2000.000000;parlib[125]=0.000000;parlib[126]=80.000000;parlib[127]=0.500000;}
       if(*nlib>=33){parlib[128]=100.000000;parlib[129]=-2.000000;parlib[130]=140.000000;parlib[131]=0.500000;}
       if(*nlib>=34){parlib[132]=733.333333;parlib[133]=-2.000000;parlib[134]=140.000000;parlib[135]=0.500000;}
       if(*nlib>=35){parlib[136]=1366.666667;parlib[137]=-2.000000;parlib[138]=140.000000;parlib[139]=0.500000;}
       if(*nlib>=36){parlib[140]=2000.000000;parlib[141]=-2.000000;parlib[142]=140.000000;parlib[143]=0.500000;}
       if(*nlib>=37){parlib[144]=100.000000;parlib[145]=-1.333333;parlib[146]=140.000000;parlib[147]=0.500000;}
       if(*nlib>=38){parlib[148]=733.333333;parlib[149]=-1.333333;parlib[150]=140.000000;parlib[151]=0.500000;}
       if(*nlib>=39){parlib[152]=1366.666667;parlib[153]=-1.333333;parlib[154]=140.000000;parlib[155]=0.500000;}
       if(*nlib>=40){parlib[156]=2000.000000;parlib[157]=-1.333333;parlib[158]=140.000000;parlib[159]=0.500000;}
       if(*nlib>=41){parlib[160]=100.000000;parlib[161]=-0.666667;parlib[162]=140.000000;parlib[163]=0.500000;}
       if(*nlib>=42){parlib[164]=733.333333;parlib[165]=-0.666667;parlib[166]=140.000000;parlib[167]=0.500000;}
       if(*nlib>=43){parlib[168]=1366.666667;parlib[169]=-0.666667;parlib[170]=140.000000;parlib[171]=0.500000;}
       if(*nlib>=44){parlib[172]=2000.000000;parlib[173]=-0.666667;parlib[174]=140.000000;parlib[175]=0.500000;}
       if(*nlib>=45){parlib[176]=100.000000;parlib[177]=0.000000;parlib[178]=140.000000;parlib[179]=0.500000;}
       if(*nlib>=46){parlib[180]=733.333333;parlib[181]=0.000000;parlib[182]=140.000000;parlib[183]=0.500000;}
       if(*nlib>=47){parlib[184]=1366.666667;parlib[185]=0.000000;parlib[186]=140.000000;parlib[187]=0.500000;}
       if(*nlib>=48){parlib[188]=2000.000000;parlib[189]=0.000000;parlib[190]=140.000000;parlib[191]=0.500000;}
       if(*nlib>=49){parlib[192]=100.000000;parlib[193]=-2.000000;parlib[194]=200.000000;parlib[195]=0.500000;}
       if(*nlib>=50){parlib[196]=733.333333;parlib[197]=-2.000000;parlib[198]=200.000000;parlib[199]=0.500000;}
       if(*nlib>=51){parlib[200]=1366.666667;parlib[201]=-2.000000;parlib[202]=200.000000;parlib[203]=0.500000;}
       if(*nlib>=52){parlib[204]=2000.000000;parlib[205]=-2.000000;parlib[206]=200.000000;parlib[207]=0.500000;}
       if(*nlib>=53){parlib[208]=100.000000;parlib[209]=-1.333333;parlib[210]=200.000000;parlib[211]=0.500000;}
       if(*nlib>=54){parlib[212]=733.333333;parlib[213]=-1.333333;parlib[214]=200.000000;parlib[215]=0.500000;}
       if(*nlib>=55){parlib[216]=1366.666667;parlib[217]=-1.333333;parlib[218]=200.000000;parlib[219]=0.500000;}
       if(*nlib>=56){parlib[220]=2000.000000;parlib[221]=-1.333333;parlib[222]=200.000000;parlib[223]=0.500000;}
       if(*nlib>=57){parlib[224]=100.000000;parlib[225]=-0.666667;parlib[226]=200.000000;parlib[227]=0.500000;}
       if(*nlib>=58){parlib[228]=733.333333;parlib[229]=-0.666667;parlib[230]=200.000000;parlib[231]=0.500000;}
       if(*nlib>=59){parlib[232]=1366.666667;parlib[233]=-0.666667;parlib[234]=200.000000;parlib[235]=0.500000;}
       if(*nlib>=60){parlib[236]=2000.000000;parlib[237]=-0.666667;parlib[238]=200.000000;parlib[239]=0.500000;}
       if(*nlib>=61){parlib[240]=100.000000;parlib[241]=0.000000;parlib[242]=200.000000;parlib[243]=0.500000;}
       if(*nlib>=62){parlib[244]=733.333333;parlib[245]=0.000000;parlib[246]=200.000000;parlib[247]=0.500000;}
       if(*nlib>=63){parlib[248]=1366.666667;parlib[249]=0.000000;parlib[250]=200.000000;parlib[251]=0.500000;}
       if(*nlib>=64){parlib[252]=2000.000000;parlib[253]=0.000000;parlib[254]=200.000000;parlib[255]=0.500000;}
       if(*nlib>=65){parlib[256]=100.000000;parlib[257]=-2.000000;parlib[258]=20.000000;parlib[259]=1.000000;}
       if(*nlib>=66){parlib[260]=733.333333;parlib[261]=-2.000000;parlib[262]=20.000000;parlib[263]=1.000000;}
       if(*nlib>=67){parlib[264]=1366.666667;parlib[265]=-2.000000;parlib[266]=20.000000;parlib[267]=1.000000;}
       if(*nlib>=68){parlib[268]=2000.000000;parlib[269]=-2.000000;parlib[270]=20.000000;parlib[271]=1.000000;}
       if(*nlib>=69){parlib[272]=100.000000;parlib[273]=-1.333333;parlib[274]=20.000000;parlib[275]=1.000000;}
       if(*nlib>=70){parlib[276]=733.333333;parlib[277]=-1.333333;parlib[278]=20.000000;parlib[279]=1.000000;}
       if(*nlib>=71){parlib[280]=1366.666667;parlib[281]=-1.333333;parlib[282]=20.000000;parlib[283]=1.000000;}
       if(*nlib>=72){parlib[284]=2000.000000;parlib[285]=-1.333333;parlib[286]=20.000000;parlib[287]=1.000000;}
       if(*nlib>=73){parlib[288]=100.000000;parlib[289]=-0.666667;parlib[290]=20.000000;parlib[291]=1.000000;}
       if(*nlib>=74){parlib[292]=733.333333;parlib[293]=-0.666667;parlib[294]=20.000000;parlib[295]=1.000000;}
       if(*nlib>=75){parlib[296]=1366.666667;parlib[297]=-0.666667;parlib[298]=20.000000;parlib[299]=1.000000;}
       if(*nlib>=76){parlib[300]=2000.000000;parlib[301]=-0.666667;parlib[302]=20.000000;parlib[303]=1.000000;}
       if(*nlib>=77){parlib[304]=100.000000;parlib[305]=0.000000;parlib[306]=20.000000;parlib[307]=1.000000;}
       if(*nlib>=78){parlib[308]=733.333333;parlib[309]=0.000000;parlib[310]=20.000000;parlib[311]=1.000000;}
       if(*nlib>=79){parlib[312]=1366.666667;parlib[313]=0.000000;parlib[314]=20.000000;parlib[315]=1.000000;}
       if(*nlib>=80){parlib[316]=2000.000000;parlib[317]=0.000000;parlib[318]=20.000000;parlib[319]=1.000000;}
       if(*nlib>=81){parlib[320]=100.000000;parlib[321]=-2.000000;parlib[322]=80.000000;parlib[323]=1.000000;}
       if(*nlib>=82){parlib[324]=733.333333;parlib[325]=-2.000000;parlib[326]=80.000000;parlib[327]=1.000000;}
       if(*nlib>=83){parlib[328]=1366.666667;parlib[329]=-2.000000;parlib[330]=80.000000;parlib[331]=1.000000;}
       if(*nlib>=84){parlib[332]=2000.000000;parlib[333]=-2.000000;parlib[334]=80.000000;parlib[335]=1.000000;}
       if(*nlib>=85){parlib[336]=100.000000;parlib[337]=-1.333333;parlib[338]=80.000000;parlib[339]=1.000000;}
       if(*nlib>=86){parlib[340]=733.333333;parlib[341]=-1.333333;parlib[342]=80.000000;parlib[343]=1.000000;}
       if(*nlib>=87){parlib[344]=1366.666667;parlib[345]=-1.333333;parlib[346]=80.000000;parlib[347]=1.000000;}
       if(*nlib>=88){parlib[348]=2000.000000;parlib[349]=-1.333333;parlib[350]=80.000000;parlib[351]=1.000000;}
       if(*nlib>=89){parlib[352]=100.000000;parlib[353]=-0.666667;parlib[354]=80.000000;parlib[355]=1.000000;}
       if(*nlib>=90){parlib[356]=733.333333;parlib[357]=-0.666667;parlib[358]=80.000000;parlib[359]=1.000000;}
       if(*nlib>=91){parlib[360]=1366.666667;parlib[361]=-0.666667;parlib[362]=80.000000;parlib[363]=1.000000;}
       if(*nlib>=92){parlib[364]=2000.000000;parlib[365]=-0.666667;parlib[366]=80.000000;parlib[367]=1.000000;}
       if(*nlib>=93){parlib[368]=100.000000;parlib[369]=0.000000;parlib[370]=80.000000;parlib[371]=1.000000;}
       if(*nlib>=94){parlib[372]=733.333333;parlib[373]=0.000000;parlib[374]=80.000000;parlib[375]=1.000000;}
       if(*nlib>=95){parlib[376]=1366.666667;parlib[377]=0.000000;parlib[378]=80.000000;parlib[379]=1.000000;}
       if(*nlib>=96){parlib[380]=2000.000000;parlib[381]=0.000000;parlib[382]=80.000000;parlib[383]=1.000000;}
       if(*nlib>=97){parlib[384]=100.000000;parlib[385]=-2.000000;parlib[386]=140.000000;parlib[387]=1.000000;}
       if(*nlib>=98){parlib[388]=733.333333;parlib[389]=-2.000000;parlib[390]=140.000000;parlib[391]=1.000000;}
       if(*nlib>=99){parlib[392]=1366.666667;parlib[393]=-2.000000;parlib[394]=140.000000;parlib[395]=1.000000;}
       if(*nlib>=100){parlib[396]=2000.000000;parlib[397]=-2.000000;parlib[398]=140.000000;parlib[399]=1.000000;}
       if(*nlib>=101){parlib[400]=100.000000;parlib[401]=-1.333333;parlib[402]=140.000000;parlib[403]=1.000000;}
       if(*nlib>=102){parlib[404]=733.333333;parlib[405]=-1.333333;parlib[406]=140.000000;parlib[407]=1.000000;}
       if(*nlib>=103){parlib[408]=1366.666667;parlib[409]=-1.333333;parlib[410]=140.000000;parlib[411]=1.000000;}
       if(*nlib>=104){parlib[412]=2000.000000;parlib[413]=-1.333333;parlib[414]=140.000000;parlib[415]=1.000000;}
       if(*nlib>=105){parlib[416]=100.000000;parlib[417]=-0.666667;parlib[418]=140.000000;parlib[419]=1.000000;}
       if(*nlib>=106){parlib[420]=733.333333;parlib[421]=-0.666667;parlib[422]=140.000000;parlib[423]=1.000000;}
       if(*nlib>=107){parlib[424]=1366.666667;parlib[425]=-0.666667;parlib[426]=140.000000;parlib[427]=1.000000;}
       if(*nlib>=108){parlib[428]=2000.000000;parlib[429]=-0.666667;parlib[430]=140.000000;parlib[431]=1.000000;}
       if(*nlib>=109){parlib[432]=100.000000;parlib[433]=0.000000;parlib[434]=140.000000;parlib[435]=1.000000;}
       if(*nlib>=110){parlib[436]=733.333333;parlib[437]=0.000000;parlib[438]=140.000000;parlib[439]=1.000000;}
       if(*nlib>=111){parlib[440]=1366.666667;parlib[441]=0.000000;parlib[442]=140.000000;parlib[443]=1.000000;}
       if(*nlib>=112){parlib[444]=2000.000000;parlib[445]=0.000000;parlib[446]=140.000000;parlib[447]=1.000000;}
       if(*nlib>=113){parlib[448]=100.000000;parlib[449]=-2.000000;parlib[450]=200.000000;parlib[451]=1.000000;}
       if(*nlib>=114){parlib[452]=733.333333;parlib[453]=-2.000000;parlib[454]=200.000000;parlib[455]=1.000000;}
       if(*nlib>=115){parlib[456]=1366.666667;parlib[457]=-2.000000;parlib[458]=200.000000;parlib[459]=1.000000;}
       if(*nlib>=116){parlib[460]=2000.000000;parlib[461]=-2.000000;parlib[462]=200.000000;parlib[463]=1.000000;}
       if(*nlib>=117){parlib[464]=100.000000;parlib[465]=-1.333333;parlib[466]=200.000000;parlib[467]=1.000000;}
       if(*nlib>=118){parlib[468]=733.333333;parlib[469]=-1.333333;parlib[470]=200.000000;parlib[471]=1.000000;}
       if(*nlib>=119){parlib[472]=1366.666667;parlib[473]=-1.333333;parlib[474]=200.000000;parlib[475]=1.000000;}
       if(*nlib>=120){parlib[476]=2000.000000;parlib[477]=-1.333333;parlib[478]=200.000000;parlib[479]=1.000000;}
       if(*nlib>=121){parlib[480]=100.000000;parlib[481]=-0.666667;parlib[482]=200.000000;parlib[483]=1.000000;}
       if(*nlib>=122){parlib[484]=733.333333;parlib[485]=-0.666667;parlib[486]=200.000000;parlib[487]=1.000000;}
       if(*nlib>=123){parlib[488]=1366.666667;parlib[489]=-0.666667;parlib[490]=200.000000;parlib[491]=1.000000;}
       if(*nlib>=124){parlib[492]=2000.000000;parlib[493]=-0.666667;parlib[494]=200.000000;parlib[495]=1.000000;}
       if(*nlib>=125){parlib[496]=100.000000;parlib[497]=0.000000;parlib[498]=200.000000;parlib[499]=1.000000;}
       if(*nlib>=126){parlib[500]=733.333333;parlib[501]=0.000000;parlib[502]=200.000000;parlib[503]=1.000000;}
       if(*nlib>=127){parlib[504]=1366.666667;parlib[505]=0.000000;parlib[506]=200.000000;parlib[507]=1.000000;}
       if(*nlib>=128){parlib[508]=2000.000000;parlib[509]=0.000000;parlib[510]=200.000000;parlib[511]=1.000000;}
       if(*nlib>=129){parlib[512]=100.000000;parlib[513]=-2.000000;parlib[514]=20.000000;parlib[515]=1.500000;}
       if(*nlib>=130){parlib[516]=733.333333;parlib[517]=-2.000000;parlib[518]=20.000000;parlib[519]=1.500000;}
       if(*nlib>=131){parlib[520]=1366.666667;parlib[521]=-2.000000;parlib[522]=20.000000;parlib[523]=1.500000;}
       if(*nlib>=132){parlib[524]=2000.000000;parlib[525]=-2.000000;parlib[526]=20.000000;parlib[527]=1.500000;}
       if(*nlib>=133){parlib[528]=100.000000;parlib[529]=-1.333333;parlib[530]=20.000000;parlib[531]=1.500000;}
       if(*nlib>=134){parlib[532]=733.333333;parlib[533]=-1.333333;parlib[534]=20.000000;parlib[535]=1.500000;}
       if(*nlib>=135){parlib[536]=1366.666667;parlib[537]=-1.333333;parlib[538]=20.000000;parlib[539]=1.500000;}
       if(*nlib>=136){parlib[540]=2000.000000;parlib[541]=-1.333333;parlib[542]=20.000000;parlib[543]=1.500000;}
       if(*nlib>=137){parlib[544]=100.000000;parlib[545]=-0.666667;parlib[546]=20.000000;parlib[547]=1.500000;}
       if(*nlib>=138){parlib[548]=733.333333;parlib[549]=-0.666667;parlib[550]=20.000000;parlib[551]=1.500000;}
       if(*nlib>=139){parlib[552]=1366.666667;parlib[553]=-0.666667;parlib[554]=20.000000;parlib[555]=1.500000;}
       if(*nlib>=140){parlib[556]=2000.000000;parlib[557]=-0.666667;parlib[558]=20.000000;parlib[559]=1.500000;}
       if(*nlib>=141){parlib[560]=100.000000;parlib[561]=0.000000;parlib[562]=20.000000;parlib[563]=1.500000;}
       if(*nlib>=142){parlib[564]=733.333333;parlib[565]=0.000000;parlib[566]=20.000000;parlib[567]=1.500000;}
       if(*nlib>=143){parlib[568]=1366.666667;parlib[569]=0.000000;parlib[570]=20.000000;parlib[571]=1.500000;}
       if(*nlib>=144){parlib[572]=2000.000000;parlib[573]=0.000000;parlib[574]=20.000000;parlib[575]=1.500000;}
       if(*nlib>=145){parlib[576]=100.000000;parlib[577]=-2.000000;parlib[578]=80.000000;parlib[579]=1.500000;}
       if(*nlib>=146){parlib[580]=733.333333;parlib[581]=-2.000000;parlib[582]=80.000000;parlib[583]=1.500000;}
       if(*nlib>=147){parlib[584]=1366.666667;parlib[585]=-2.000000;parlib[586]=80.000000;parlib[587]=1.500000;}
       if(*nlib>=148){parlib[588]=2000.000000;parlib[589]=-2.000000;parlib[590]=80.000000;parlib[591]=1.500000;}
       if(*nlib>=149){parlib[592]=100.000000;parlib[593]=-1.333333;parlib[594]=80.000000;parlib[595]=1.500000;}
       if(*nlib>=150){parlib[596]=733.333333;parlib[597]=-1.333333;parlib[598]=80.000000;parlib[599]=1.500000;}
       if(*nlib>=151){parlib[600]=1366.666667;parlib[601]=-1.333333;parlib[602]=80.000000;parlib[603]=1.500000;}
       if(*nlib>=152){parlib[604]=2000.000000;parlib[605]=-1.333333;parlib[606]=80.000000;parlib[607]=1.500000;}
       if(*nlib>=153){parlib[608]=100.000000;parlib[609]=-0.666667;parlib[610]=80.000000;parlib[611]=1.500000;}
       if(*nlib>=154){parlib[612]=733.333333;parlib[613]=-0.666667;parlib[614]=80.000000;parlib[615]=1.500000;}
       if(*nlib>=155){parlib[616]=1366.666667;parlib[617]=-0.666667;parlib[618]=80.000000;parlib[619]=1.500000;}
       if(*nlib>=156){parlib[620]=2000.000000;parlib[621]=-0.666667;parlib[622]=80.000000;parlib[623]=1.500000;}
       if(*nlib>=157){parlib[624]=100.000000;parlib[625]=0.000000;parlib[626]=80.000000;parlib[627]=1.500000;}
       if(*nlib>=158){parlib[628]=733.333333;parlib[629]=0.000000;parlib[630]=80.000000;parlib[631]=1.500000;}
       if(*nlib>=159){parlib[632]=1366.666667;parlib[633]=0.000000;parlib[634]=80.000000;parlib[635]=1.500000;}
       if(*nlib>=160){parlib[636]=2000.000000;parlib[637]=0.000000;parlib[638]=80.000000;parlib[639]=1.500000;}
       if(*nlib>=161){parlib[640]=100.000000;parlib[641]=-2.000000;parlib[642]=140.000000;parlib[643]=1.500000;}
       if(*nlib>=162){parlib[644]=733.333333;parlib[645]=-2.000000;parlib[646]=140.000000;parlib[647]=1.500000;}
       if(*nlib>=163){parlib[648]=1366.666667;parlib[649]=-2.000000;parlib[650]=140.000000;parlib[651]=1.500000;}
       if(*nlib>=164){parlib[652]=2000.000000;parlib[653]=-2.000000;parlib[654]=140.000000;parlib[655]=1.500000;}
       if(*nlib>=165){parlib[656]=100.000000;parlib[657]=-1.333333;parlib[658]=140.000000;parlib[659]=1.500000;}
       if(*nlib>=166){parlib[660]=733.333333;parlib[661]=-1.333333;parlib[662]=140.000000;parlib[663]=1.500000;}
       if(*nlib>=167){parlib[664]=1366.666667;parlib[665]=-1.333333;parlib[666]=140.000000;parlib[667]=1.500000;}
       if(*nlib>=168){parlib[668]=2000.000000;parlib[669]=-1.333333;parlib[670]=140.000000;parlib[671]=1.500000;}
       if(*nlib>=169){parlib[672]=100.000000;parlib[673]=-0.666667;parlib[674]=140.000000;parlib[675]=1.500000;}
       if(*nlib>=170){parlib[676]=733.333333;parlib[677]=-0.666667;parlib[678]=140.000000;parlib[679]=1.500000;}
       if(*nlib>=171){parlib[680]=1366.666667;parlib[681]=-0.666667;parlib[682]=140.000000;parlib[683]=1.500000;}
       if(*nlib>=172){parlib[684]=2000.000000;parlib[685]=-0.666667;parlib[686]=140.000000;parlib[687]=1.500000;}
       if(*nlib>=173){parlib[688]=100.000000;parlib[689]=0.000000;parlib[690]=140.000000;parlib[691]=1.500000;}
       if(*nlib>=174){parlib[692]=733.333333;parlib[693]=0.000000;parlib[694]=140.000000;parlib[695]=1.500000;}
       if(*nlib>=175){parlib[696]=1366.666667;parlib[697]=0.000000;parlib[698]=140.000000;parlib[699]=1.500000;}
       if(*nlib>=176){parlib[700]=2000.000000;parlib[701]=0.000000;parlib[702]=140.000000;parlib[703]=1.500000;}
       if(*nlib>=177){parlib[704]=100.000000;parlib[705]=-2.000000;parlib[706]=200.000000;parlib[707]=1.500000;}
       if(*nlib>=178){parlib[708]=733.333333;parlib[709]=-2.000000;parlib[710]=200.000000;parlib[711]=1.500000;}
       if(*nlib>=179){parlib[712]=1366.666667;parlib[713]=-2.000000;parlib[714]=200.000000;parlib[715]=1.500000;}
       if(*nlib>=180){parlib[716]=2000.000000;parlib[717]=-2.000000;parlib[718]=200.000000;parlib[719]=1.500000;}
       if(*nlib>=181){parlib[720]=100.000000;parlib[721]=-1.333333;parlib[722]=200.000000;parlib[723]=1.500000;}
       if(*nlib>=182){parlib[724]=733.333333;parlib[725]=-1.333333;parlib[726]=200.000000;parlib[727]=1.500000;}
       if(*nlib>=183){parlib[728]=1366.666667;parlib[729]=-1.333333;parlib[730]=200.000000;parlib[731]=1.500000;}
       if(*nlib>=184){parlib[732]=2000.000000;parlib[733]=-1.333333;parlib[734]=200.000000;parlib[735]=1.500000;}
       if(*nlib>=185){parlib[736]=100.000000;parlib[737]=-0.666667;parlib[738]=200.000000;parlib[739]=1.500000;}
       if(*nlib>=186){parlib[740]=733.333333;parlib[741]=-0.666667;parlib[742]=200.000000;parlib[743]=1.500000;}
       if(*nlib>=187){parlib[744]=1366.666667;parlib[745]=-0.666667;parlib[746]=200.000000;parlib[747]=1.500000;}
       if(*nlib>=188){parlib[748]=2000.000000;parlib[749]=-0.666667;parlib[750]=200.000000;parlib[751]=1.500000;}
       if(*nlib>=189){parlib[752]=100.000000;parlib[753]=0.000000;parlib[754]=200.000000;parlib[755]=1.500000;}
       if(*nlib>=190){parlib[756]=733.333333;parlib[757]=0.000000;parlib[758]=200.000000;parlib[759]=1.500000;}
       if(*nlib>=191){parlib[760]=1366.666667;parlib[761]=0.000000;parlib[762]=200.000000;parlib[763]=1.500000;}
       if(*nlib>=192){parlib[764]=2000.000000;parlib[765]=0.000000;parlib[766]=200.000000;parlib[767]=1.500000;}
       if(*nlib>=193){parlib[768]=100.000000;parlib[769]=-2.000000;parlib[770]=20.000000;parlib[771]=2.000000;}
       if(*nlib>=194){parlib[772]=733.333333;parlib[773]=-2.000000;parlib[774]=20.000000;parlib[775]=2.000000;}
       if(*nlib>=195){parlib[776]=1366.666667;parlib[777]=-2.000000;parlib[778]=20.000000;parlib[779]=2.000000;}
       if(*nlib>=196){parlib[780]=2000.000000;parlib[781]=-2.000000;parlib[782]=20.000000;parlib[783]=2.000000;}
       if(*nlib>=197){parlib[784]=100.000000;parlib[785]=-1.333333;parlib[786]=20.000000;parlib[787]=2.000000;}
       if(*nlib>=198){parlib[788]=733.333333;parlib[789]=-1.333333;parlib[790]=20.000000;parlib[791]=2.000000;}
       if(*nlib>=199){parlib[792]=1366.666667;parlib[793]=-1.333333;parlib[794]=20.000000;parlib[795]=2.000000;}
       if(*nlib>=200){parlib[796]=2000.000000;parlib[797]=-1.333333;parlib[798]=20.000000;parlib[799]=2.000000;}
       if(*nlib>=201){parlib[800]=100.000000;parlib[801]=-0.666667;parlib[802]=20.000000;parlib[803]=2.000000;}
       if(*nlib>=202){parlib[804]=733.333333;parlib[805]=-0.666667;parlib[806]=20.000000;parlib[807]=2.000000;}
       if(*nlib>=203){parlib[808]=1366.666667;parlib[809]=-0.666667;parlib[810]=20.000000;parlib[811]=2.000000;}
       if(*nlib>=204){parlib[812]=2000.000000;parlib[813]=-0.666667;parlib[814]=20.000000;parlib[815]=2.000000;}
       if(*nlib>=205){parlib[816]=100.000000;parlib[817]=0.000000;parlib[818]=20.000000;parlib[819]=2.000000;}
       if(*nlib>=206){parlib[820]=733.333333;parlib[821]=0.000000;parlib[822]=20.000000;parlib[823]=2.000000;}
       if(*nlib>=207){parlib[824]=1366.666667;parlib[825]=0.000000;parlib[826]=20.000000;parlib[827]=2.000000;}
       if(*nlib>=208){parlib[828]=2000.000000;parlib[829]=0.000000;parlib[830]=20.000000;parlib[831]=2.000000;}
       if(*nlib>=209){parlib[832]=100.000000;parlib[833]=-2.000000;parlib[834]=80.000000;parlib[835]=2.000000;}
       if(*nlib>=210){parlib[836]=733.333333;parlib[837]=-2.000000;parlib[838]=80.000000;parlib[839]=2.000000;}
       if(*nlib>=211){parlib[840]=1366.666667;parlib[841]=-2.000000;parlib[842]=80.000000;parlib[843]=2.000000;}
       if(*nlib>=212){parlib[844]=2000.000000;parlib[845]=-2.000000;parlib[846]=80.000000;parlib[847]=2.000000;}
       if(*nlib>=213){parlib[848]=100.000000;parlib[849]=-1.333333;parlib[850]=80.000000;parlib[851]=2.000000;}
       if(*nlib>=214){parlib[852]=733.333333;parlib[853]=-1.333333;parlib[854]=80.000000;parlib[855]=2.000000;}
       if(*nlib>=215){parlib[856]=1366.666667;parlib[857]=-1.333333;parlib[858]=80.000000;parlib[859]=2.000000;}
       if(*nlib>=216){parlib[860]=2000.000000;parlib[861]=-1.333333;parlib[862]=80.000000;parlib[863]=2.000000;}
       if(*nlib>=217){parlib[864]=100.000000;parlib[865]=-0.666667;parlib[866]=80.000000;parlib[867]=2.000000;}
       if(*nlib>=218){parlib[868]=733.333333;parlib[869]=-0.666667;parlib[870]=80.000000;parlib[871]=2.000000;}
       if(*nlib>=219){parlib[872]=1366.666667;parlib[873]=-0.666667;parlib[874]=80.000000;parlib[875]=2.000000;}
       if(*nlib>=220){parlib[876]=2000.000000;parlib[877]=-0.666667;parlib[878]=80.000000;parlib[879]=2.000000;}
       if(*nlib>=221){parlib[880]=100.000000;parlib[881]=0.000000;parlib[882]=80.000000;parlib[883]=2.000000;}
       if(*nlib>=222){parlib[884]=733.333333;parlib[885]=0.000000;parlib[886]=80.000000;parlib[887]=2.000000;}
       if(*nlib>=223){parlib[888]=1366.666667;parlib[889]=0.000000;parlib[890]=80.000000;parlib[891]=2.000000;}
       if(*nlib>=224){parlib[892]=2000.000000;parlib[893]=0.000000;parlib[894]=80.000000;parlib[895]=2.000000;}
       if(*nlib>=225){parlib[896]=100.000000;parlib[897]=-2.000000;parlib[898]=140.000000;parlib[899]=2.000000;}
       if(*nlib>=226){parlib[900]=733.333333;parlib[901]=-2.000000;parlib[902]=140.000000;parlib[903]=2.000000;}
       if(*nlib>=227){parlib[904]=1366.666667;parlib[905]=-2.000000;parlib[906]=140.000000;parlib[907]=2.000000;}
       if(*nlib>=228){parlib[908]=2000.000000;parlib[909]=-2.000000;parlib[910]=140.000000;parlib[911]=2.000000;}
       if(*nlib>=229){parlib[912]=100.000000;parlib[913]=-1.333333;parlib[914]=140.000000;parlib[915]=2.000000;}
       if(*nlib>=230){parlib[916]=733.333333;parlib[917]=-1.333333;parlib[918]=140.000000;parlib[919]=2.000000;}
       if(*nlib>=231){parlib[920]=1366.666667;parlib[921]=-1.333333;parlib[922]=140.000000;parlib[923]=2.000000;}
       if(*nlib>=232){parlib[924]=2000.000000;parlib[925]=-1.333333;parlib[926]=140.000000;parlib[927]=2.000000;}
       if(*nlib>=233){parlib[928]=100.000000;parlib[929]=-0.666667;parlib[930]=140.000000;parlib[931]=2.000000;}
       if(*nlib>=234){parlib[932]=733.333333;parlib[933]=-0.666667;parlib[934]=140.000000;parlib[935]=2.000000;}
       if(*nlib>=235){parlib[936]=1366.666667;parlib[937]=-0.666667;parlib[938]=140.000000;parlib[939]=2.000000;}
       if(*nlib>=236){parlib[940]=2000.000000;parlib[941]=-0.666667;parlib[942]=140.000000;parlib[943]=2.000000;}
       if(*nlib>=237){parlib[944]=100.000000;parlib[945]=0.000000;parlib[946]=140.000000;parlib[947]=2.000000;}
       if(*nlib>=238){parlib[948]=733.333333;parlib[949]=0.000000;parlib[950]=140.000000;parlib[951]=2.000000;}
       if(*nlib>=239){parlib[952]=1366.666667;parlib[953]=0.000000;parlib[954]=140.000000;parlib[955]=2.000000;}
       if(*nlib>=240){parlib[956]=2000.000000;parlib[957]=0.000000;parlib[958]=140.000000;parlib[959]=2.000000;}
       if(*nlib>=241){parlib[960]=100.000000;parlib[961]=-2.000000;parlib[962]=200.000000;parlib[963]=2.000000;}
       if(*nlib>=242){parlib[964]=733.333333;parlib[965]=-2.000000;parlib[966]=200.000000;parlib[967]=2.000000;}
       if(*nlib>=243){parlib[968]=1366.666667;parlib[969]=-2.000000;parlib[970]=200.000000;parlib[971]=2.000000;}
       if(*nlib>=244){parlib[972]=2000.000000;parlib[973]=-2.000000;parlib[974]=200.000000;parlib[975]=2.000000;}
       if(*nlib>=245){parlib[976]=100.000000;parlib[977]=-1.333333;parlib[978]=200.000000;parlib[979]=2.000000;}
       if(*nlib>=246){parlib[980]=733.333333;parlib[981]=-1.333333;parlib[982]=200.000000;parlib[983]=2.000000;}
       if(*nlib>=247){parlib[984]=1366.666667;parlib[985]=-1.333333;parlib[986]=200.000000;parlib[987]=2.000000;}
       if(*nlib>=248){parlib[988]=2000.000000;parlib[989]=-1.333333;parlib[990]=200.000000;parlib[991]=2.000000;}
       if(*nlib>=249){parlib[992]=100.000000;parlib[993]=-0.666667;parlib[994]=200.000000;parlib[995]=2.000000;}
       if(*nlib>=250){parlib[996]=733.333333;parlib[997]=-0.666667;parlib[998]=200.000000;parlib[999]=2.000000;}
       if(*nlib>=251){parlib[1000]=1366.666667;parlib[1001]=-0.666667;parlib[1002]=200.000000;parlib[1003]=2.000000;}
       if(*nlib>=252){parlib[1004]=2000.000000;parlib[1005]=-0.666667;parlib[1006]=200.000000;parlib[1007]=2.000000;}
       if(*nlib>=253){parlib[1008]=100.000000;parlib[1009]=0.000000;parlib[1010]=200.000000;parlib[1011]=2.000000;}
       if(*nlib>=254){parlib[1012]=733.333333;parlib[1013]=0.000000;parlib[1014]=200.000000;parlib[1015]=2.000000;}
       if(*nlib>=255){parlib[1016]=1366.666667;parlib[1017]=0.000000;parlib[1018]=200.000000;parlib[1019]=2.000000;}
       if(*nlib>=256){parlib[1020]=2000.000000;parlib[1021]=0.000000;parlib[1022]=200.000000;parlib[1023]=2.000000;}
	return;
}
