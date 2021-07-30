#ifndef SPECPOLY
#define SPECPOLY
#include <complex.h>
#include <stdlib.h>
void spec_poly(double complex** poly,double complex** exact_roots,unsigned int* deg,const unsigned int poly_num)
{
	switch(poly_num)
	{
		// Mandelbrot Polynomial, theoretical limiting accuracy 2.53E-5
		case 1:
		// set degree
		*deg = 63;
		// allocate space for poly and exact_roots
		*poly = (double complex*)malloc((*deg+1)*sizeof(double complex));
		*exact_roots = (double complex*)malloc((*deg)*sizeof(double complex));
		// store Mandelbrot exact_roots, rounded to double precision
	    (*exact_roots)[0] = CMPLX(-0.199909568232701847321062999922E+1,0E+0);
	    (*exact_roots)[1] = CMPLX(-0.199181417254912221573256094986E+1,0E+0);
	    (*exact_roots)[2] = CMPLX(-0.197717958700625738734608852066E+1,0E+0);
	    (*exact_roots)[3] = CMPLX(-0.195370589428439624542762219901E+1,0E+0);
	    (*exact_roots)[4] = CMPLX(-0.192714770936395026246006818895E+1,0E+0);
	    (*exact_roots)[5] = CMPLX(-0.188480357158668179232947809292E+1,0E+0);
	    (*exact_roots)[6] = CMPLX(-0.183231520275122919208489752604E+1,0E+0);
	    (*exact_roots)[7] = CMPLX(-0.176926167027683114607548022864E+1,-0.5691950039560031530490018730E-1);
	    (*exact_roots)[8] = CMPLX(-0.176926167027683114607548022864E+1,0.5691950039560031530490018730E-1);
	    (*exact_roots)[9] = CMPLX(-0.167406609147478797156529602917E+1,0E+0);
	    (*exact_roots)[10] = CMPLX(-0.157488913975230096981996555250E+1,0E+0);
	    (*exact_roots)[11] = CMPLX(-0.140844648574007265491757700881E+1,-0.13617199730465991568470779361E+0);
	    (*exact_roots)[12] = CMPLX(-0.140844648574007265491757700881E+1,0.13617199730465991568470779361E+0);
	    (*exact_roots)[13] = CMPLX(-0.1292558061033522087164184706361E+1,-0.438198816086631837129737124327E+0);
	    (*exact_roots)[14] = CMPLX(-0.1292558061033522087164184706361E+1,0.438198816086631837129737124327E+0);
	    (*exact_roots)[15] = CMPLX(-0.1262287281438472543010111941208E+1,-0.408104324112690383290160657426E+0);
	    (*exact_roots)[16] = CMPLX(-0.1262287281438472543010111941208E+1,0.408104324112690383290160657426E+0);
	    (*exact_roots)[17] = CMPLX(-0.1252735884012037946295811002564E+1,-0.342470647889750893861875786871E+0);
	    (*exact_roots)[18] = CMPLX(-0.1252735884012037946295811002564E+1,0.342470647889750893861875786871E+0);
	    (*exact_roots)[19] = CMPLX(-0.1028193852454817599302497455967E+1,-0.361376517118561592479460832998E+0);
	    (*exact_roots)[20] = CMPLX(-0.1028193852454817599302497455967E+1,0.361376517118561592479460832998E+0);
	    (*exact_roots)[21] = CMPLX(-0.623532485956252757990016587001E+0,-0.681064414225239608090835812687E+0);
	    (*exact_roots)[22] = CMPLX(-0.623532485956252757990016587001E+0,0.681064414225239608090835812687E+0);
	    (*exact_roots)[23] = CMPLX(-0.622436295041293587960163506947E+0,-0.424878436475629184311574438805E+0);
	    (*exact_roots)[24] = CMPLX(-0.622436295041293587960163506947E+0,0.424878436475629184311574438805E+0);
	    (*exact_roots)[25] = CMPLX(-0.530827804859942728921477297120E+0,-0.668288725559205771444092465565E+0);
	    (*exact_roots)[26] = CMPLX(-0.530827804859942728921477297120E+0,0.668288725559205771444092465565E+0);
	    (*exact_roots)[27] = CMPLX(-0.272102461488938894219383324518E+0,-0.842364690294128145503155708243E+0);
	    (*exact_roots)[28] = CMPLX(-0.272102461488938894219383324518E+0,0.842364690294128145503155708243E+0);
	    (*exact_roots)[29] = CMPLX(-0.224915951286740054685326255204E+0,-0.1116260157454991835001268254245E+1);
	    (*exact_roots)[30] = CMPLX(-0.224915951286740054685326255204E+0,0.1116260157454991835001268254245E+1);
	    (*exact_roots)[31] = CMPLX(-0.207283835455666412824133850187E+0,-0.1117480772494962911373775673122E+1);
	    (*exact_roots)[32] = CMPLX(-0.207283835455666412824133850187E+0,0.1117480772494962911373775673122E+1);
	    (*exact_roots)[33] = CMPLX(-0.174578221135716969451566432662E+0,-0.1071427671454031189229646310220E+1);
	    (*exact_roots)[34] = CMPLX(-0.174578221135716969451566432662E+0,0.1071427671454031189229646310220E+1);
	    (*exact_roots)[35] = CMPLX(-0.157516053475965356164335109645E+0,-0.1109006514113607177971751986155E+1);
	    (*exact_roots)[36] = CMPLX(-0.157516053475965356164335109645E+0,0.1109006514113607177971751986155E+1);
	    (*exact_roots)[37] = CMPLX(-0.127499973546363000199539565346E+0,-0.987460909489456792207407680793E+0);
	    (*exact_roots)[38] = CMPLX(-0.127499973546363000199539565346E+0,0.987460909489456792207407680793E+0);
	    (*exact_roots)[39] = CMPLX(-0.14233481920354066767761845314E-1,-0.1032914775213644109395013402655E+1);
	    (*exact_roots)[40] = CMPLX(-0.14233481920354066767761845314E-1,0.1032914775213644109395013402655E+1);
	    (*exact_roots)[41] = CMPLX(-0.6983568496261391817966491075E-2,-0.1003603862288289548530704966951E+1);
	    (*exact_roots)[42] = CMPLX(-0.6983568496261391817966491075E-2,0.1003603862288289548530704966951E+1);
	    (*exact_roots)[43] = CMPLX(0.14895466603687646529815779209E-1,-0.848148761908416527719331111783E+0);
	    (*exact_roots)[44] = CMPLX(0.14895466603687646529815779209E-1,0.848148761908416527719331111783E+0);
	    (*exact_roots)[45] = CMPLX(0.121192786105906486314704443411E+0,-0.610611692210754211675387244150E+0);
	    (*exact_roots)[46] = CMPLX(0.121192786105906486314704443411E+0,0.610611692210754211675387244150E+0);
	    (*exact_roots)[47] = CMPLX(0.352482539722363278193253964052E+0,-0.698337239583330331258141954760E+0);
	    (*exact_roots)[48] = CMPLX(0.352482539722363278193253964052E+0,0.698337239583330331258141954760E+0);
	    (*exact_roots)[49] = CMPLX(0.376008681846767559704804317729E+0,-0.144749371321632864747110182013E+0);
	    (*exact_roots)[50] = CMPLX(0.376008681846767559704804317729E+0,0.144749371321632864747110182013E+0);
	    (*exact_roots)[51] = CMPLX(0.376893240379311323690004017969E+0,-0.678568693190448141957540792997E+0);
	    (*exact_roots)[52] = CMPLX(0.376893240379311323690004017969E+0,0.678568693190448141957540792997E+0);
	    (*exact_roots)[53] = CMPLX(0.386539176596158026508293086904E+0,-0.569324711303102903213792357135E+0);
	    (*exact_roots)[54] = CMPLX(0.386539176596158026508293086904E+0,0.569324711303102903213792357135E+0);
	    (*exact_roots)[55] = CMPLX(0.412916024722700479197334566383E+0,-0.614806760143385694954549720401E+0);
	    (*exact_roots)[56] = CMPLX(0.412916024722700479197334566383E+0,0.614806760143385694954549720401E+0);
	    (*exact_roots)[57] = CMPLX(0.432376192641994507824669648087E+0,-0.226759904435348618697876559972E+0);
	    (*exact_roots)[58] = CMPLX(0.432376192641994507824669648087E+0,0.226759904435348618697876559972E+0);
	    (*exact_roots)[59] = CMPLX(0.452774498724915493508803077733E+0,-0.396170128033165002412596877271E+0);
	    (*exact_roots)[60] = CMPLX(0.452774498724915493508803077733E+0,0.396170128033165002412596877271E+0);
	    (*exact_roots)[61] = CMPLX(0.456823285823316651283953236253E+0,-0.347758700883481983632188723200E+0);
	    (*exact_roots)[62] = CMPLX(0.456823285823316651283953236253E+0,0.347758700883481983632188723200E+0);
		// store Mandelbrot poly, rounded to double precision
	    (*poly)[0] = CMPLX(1E+0,0E+0);
	    (*poly)[1] = CMPLX(1E+0,0E+0);
	    (*poly)[2] = CMPLX(2E+0,0E+0);
	    (*poly)[3] = CMPLX(5E+0,0E+0);
	    (*poly)[4] = CMPLX(14E+0,0E+0);
	    (*poly)[5] = CMPLX(42E+0,0E+0);
	    (*poly)[6] = CMPLX(132E+0,0E+0);
	    (*poly)[7] = CMPLX(365E+0,0E+0);
	    (*poly)[8] = CMPLX(950E+0,0E+0);
	    (*poly)[9] = CMPLX(2398E+0,0E+0);
	    (*poly)[10] = CMPLX(5916E+0,0E+0);
	    (*poly)[11] = CMPLX(14290E+0,0E+0);
	    (*poly)[12] = CMPLX(33708E+0,0E+0);
	    (*poly)[13] = CMPLX(77684E+0,0E+0);
	    (*poly)[14] = CMPLX(175048E+0,0E+0);
	    (*poly)[15] = CMPLX(385741E+0,0E+0);
	    (*poly)[16] = CMPLX(831014E+0,0E+0);
	    (*poly)[17] = CMPLX(1749654E+0,0E+0);
	    (*poly)[18] = CMPLX(3598964E+0,0E+0);
	    (*poly)[19] = CMPLX(7228014E+0,0E+0);
	    (*poly)[20] = CMPLX(14162220E+0,0E+0);
	    (*poly)[21] = CMPLX(27049196E+0,0E+0);
	    (*poly)[22] = CMPLX(50323496E+0,0E+0);
	    (*poly)[23] = CMPLX(91143114E+0,0E+0);
	    (*poly)[24] = CMPLX(160617860E+0,0E+0);
	    (*poly)[25] = CMPLX(275276716E+0,0E+0);
	    (*poly)[26] = CMPLX(458591432E+0,0E+0);
	    (*poly)[27] = CMPLX(742179284E+0,0E+0);
	    (*poly)[28] = CMPLX(1166067016E+0,0E+0);
	    (*poly)[29] = CMPLX(1777171560E+0,0E+0);
	    (*poly)[30] = CMPLX(2625062128E+0,0E+0);
	    (*poly)[31] = CMPLX(3754272037E+0,0E+0);
	    (*poly)[32] = CMPLX(5193067630E+0,0E+0);
	    (*poly)[33] = CMPLX(6939692682E+0,0E+0);
	    (*poly)[34] = CMPLX(8948546308E+0,0E+0);
	    (*poly)[35] = CMPLX(11120136162E+0,0E+0);
	    (*poly)[36] = CMPLX(13299362332E+0,0E+0);
	    (*poly)[37] = CMPLX(15286065700E+0,0E+0);
	    (*poly)[38] = CMPLX(16859410792E+0,0E+0);
	    (*poly)[39] = CMPLX(17813777994E+0,0E+0);
	    (*poly)[40] = CMPLX(17999433372E+0,0E+0);
	    (*poly)[41] = CMPLX(17357937708E+0,0E+0);
	    (*poly)[42] = CMPLX(15941684776E+0,0E+0);
	    (*poly)[43] = CMPLX(13910043524E+0,0E+0);
	    (*poly)[44] = CMPLX(11500901864E+0,0E+0);
	    (*poly)[45] = CMPLX(8984070856E+0,0E+0);
	    (*poly)[46] = CMPLX(6609143792E+0,0E+0);
	    (*poly)[47] = CMPLX(4562339774E+0,0E+0);
	    (*poly)[48] = CMPLX(2943492972E+0,0E+0);
	    (*poly)[49] = CMPLX(1766948340E+0,0E+0);
	    (*poly)[50] = CMPLX(981900168E+0,0E+0);
	    (*poly)[51] = CMPLX(502196500E+0,0E+0);
	    (*poly)[52] = CMPLX(234813592E+0,0E+0);
	    (*poly)[53] = CMPLX(99582920E+0,0E+0);
	    (*poly)[54] = CMPLX(37945904E+0,0E+0);
	    (*poly)[55] = CMPLX(12843980E+0,0E+0);
	    (*poly)[56] = CMPLX(3807704E+0,0E+0);
	    (*poly)[57] = CMPLX(971272E+0,0E+0);
	    (*poly)[58] = CMPLX(208336E+0,0E+0);
	    (*poly)[59] = CMPLX(36440E+0,0E+0);
	    (*poly)[60] = CMPLX(4976E+0,0E+0);
	    (*poly)[61] = CMPLX(496E+0,0E+0);
	    (*poly)[62] = CMPLX(32E+0,0E+0);
	    (*poly)[63] = CMPLX(1E+0,0E+0);
		break;
		// Kameny Polynomial c=10, theoretical limiting accuracy 1.11E-16
		case 2:
		// set degree
		*deg = 9;
		// allocate space for poly and exact_roots
		*poly = (double complex*)malloc((*deg+1)*sizeof(double complex));
		*exact_roots = (double complex*)malloc((*deg)*sizeof(double complex));
        // store Kameny exact_roots, rounded to double precision
        (*exact_roots)[0] = CMPLX(-0.250708394126194522431084083966E+1,0E+0);
        (*exact_roots)[1] = CMPLX(-0.77475933255040504312683246509E+0,-0.239347474286129783412456829950E+1);
        (*exact_roots)[2] = CMPLX(-0.77475933255040504312683246509E+0,0.239347474286129783412456829950E+1);
        (*exact_roots)[3] = CMPLX(-0.173313478195944654943157024747E+0,0E+0);
        (*exact_roots)[4] = CMPLX(-0.173097223325097844852013748504E+0,0E+0);
        (*exact_roots)[5] = CMPLX(0.173204810760521109164135307321E+0,-0.108125591879847889809932772E-3);
        (*exact_roots)[6] = CMPLX(0.173204810760521109164135307321E+0,0.108125591879847889809932772E-3);
        (*exact_roots)[7] = CMPLX(0.202830184318137779601570296422E+1,-0.147928157900441456648083759833E+1);
        (*exact_roots)[8] = CMPLX(0.202830184318137779601570296422E+1,0.147928157900441456648083759833E+1);
		// store Kameny poly, rounded to double precision
        (*poly)[0] = CMPLX(9E+0,0E+0);
        (*poly)[1] = CMPLX(0E+0,0E+0);
        (*poly)[2] = CMPLX(-600E+0,0E+0);
        (*poly)[3] = CMPLX(0E+0,0E+0);
        (*poly)[4] = CMPLX(10000E+0,0E+0);
        (*poly)[5] = CMPLX(0E+0,0E+0);
        (*poly)[6] = CMPLX(0E+0,0E+0);
        (*poly)[7] = CMPLX(0E+0,0E+0);
        (*poly)[8] = CMPLX(0E+0,0E+0);
        (*poly)[9] = CMPLX(100E+0,0E+0);
		break;
		// Kameny Polynomial c=10^3, theoretical limiting accuracy 1.11E-16
		case 3:
		// set degree
		*deg = 9;
		// allocate space for poly and exact_roots
		*poly = (double complex*)malloc((*deg+1)*sizeof(double complex));
		*exact_roots = (double complex*)malloc((*deg)*sizeof(double complex));
        // store Kameny exact_roots, rounded to double precision
        (*exact_roots)[0] = CMPLX(-0.158489318488962525196907753240E+2,0E+0);
        (*exact_roots)[1] = CMPLX(-0.48975892839992994551324676957E+1,-0.150732300552288401309352587876E+2);
        (*exact_roots)[2] = CMPLX(-0.48975892839992994551324676957E+1,0.150732300552288401309352587876E+2);
        (*exact_roots)[3] = CMPLX(-0.173205080767700380719051028416E-2,0E+0);
        (*exact_roots)[4] = CMPLX(-0.173205080746075077991838239885E-2,0E+0);
        (*exact_roots)[5] = CMPLX(0.173205080756887729350044634151E-2,-0.10812651363606394264E-12);
        (*exact_roots)[6] = CMPLX(0.173205080756887729350044634151E-2,0.10812651363606394264E-12);
        (*exact_roots)[7] = CMPLX(0.128220552084474257149779093577E+2,-0.93157684943778791903260859371E+1);
        (*exact_roots)[8] = CMPLX(0.128220552084474257149779093577E+2,0.93157684943778791903260859371E+1);
		// store Kameny poly, rounded to double precision
        (*poly)[0] = CMPLX(9E+0,0E+0);
        (*poly)[1] = CMPLX(0E+0,0E+0);
        (*poly)[2] = CMPLX(-6000000E+0,0E+0);
        (*poly)[3] = CMPLX(0E+0,0E+0);
        (*poly)[4] = CMPLX(1000000000000E+0,0E+0);
        (*poly)[5] = CMPLX(0E+0,0E+0);
        (*poly)[6] = CMPLX(0E+0,0E+0);
        (*poly)[7] = CMPLX(0E+0,0E+0);
        (*poly)[8] = CMPLX(0E+0,0E+0);
        (*poly)[9] = CMPLX(1000000E+0,0E+0);
		break;
		// Kameny Polynomial c=10^6, theoretical limiting accuracy 8.09E-9
		case 4:
		// set degree
		*deg = 9;
		// allocate space for poly and exact_roots
		*poly = (double complex*)malloc((*deg+1)*sizeof(double complex));
		*exact_roots = (double complex*)malloc((*deg)*sizeof(double complex));
        // store Kameny exact_roots, rounded to double precision
        (*exact_roots)[0] = CMPLX(-0.251188643150958006331217160138E3,0E0);
        (*exact_roots)[1] = CMPLX(-0.77621559527630266049907587459E2,-0.238894595888056617951698213041E3);
        (*exact_roots)[2] = CMPLX(-0.77621559527630266049907587459E2,0.238894595888056617951698213041E3);
        (*exact_roots)[3] = CMPLX(-0.173205080756887729353086560209E-5,0E0);
        (*exact_roots)[4] = CMPLX(-0.173205080756887729352402708092E-5,0E0);
        (*exact_roots)[5] = CMPLX(0.173205080756887729352744634151E-5,-0.341926059E-26);
        (*exact_roots)[6] = CMPLX(0.173205080756887729352744634151E-5,0.341926059E-26);
        (*exact_roots)[7] = CMPLX(0.203215881103109269215516167528E3,-0.147644979987489859882850769633E3);
        (*exact_roots)[8] = CMPLX(0.203215881103109269215516167528E3,0.147644979987489859882850769633E3);
		// store Kameny poly, rounded to double precision
        (*poly)[0] = CMPLX(9E0,0E0);
        (*poly)[1] = CMPLX(0E0,0E0);
        (*poly)[2] = CMPLX(-6000000000000E0,0E0);
        (*poly)[3] = CMPLX(0E0,0E0);
        (*poly)[4] = CMPLX(1000000000000000000000000E0,0E0);
        (*poly)[5] = CMPLX(0E0,0E0);
        (*poly)[6] = CMPLX(0E0,0E0);
        (*poly)[7] = CMPLX(0E0,0E0);
        (*poly)[8] = CMPLX(0E0,0E0);
        (*poly)[9] = CMPLX(1000000000000E0,0E0);
		break;
		// Multiplicity I
		case 5:
		// set degree
		*deg = 55;
		// allocate space for poly and exact_roots
		*poly = (double complex*)malloc((*deg+1)*sizeof(double complex));
		*exact_roots = (double complex*)malloc((*deg)*sizeof(double complex));
		// store Multiplicity I exact_roots, rounded to double precision
        (*exact_roots)[0] = CMPLX(-0.100000000000000000000000000000E1,0E0);
        (*exact_roots)[1] = CMPLX(-0.100000000000000000000000000000E1,0E0);
        (*exact_roots)[2] = CMPLX(-0.100000000000000000000000000000E1,0E0);
        (*exact_roots)[3] = CMPLX(-0.100000000000000000000000000000E1,0E0);
        (*exact_roots)[4] = CMPLX(-0.100000000000000000000000000000E1,0E0);
        (*exact_roots)[5] = CMPLX(-0.951677511788434760780566076462E0,-0.156577036126765257092467586019E0);
        (*exact_roots)[6] = CMPLX(-0.951677511788434760780566076462E0,0.156577036126765257092467586019E0);
        (*exact_roots)[7] = CMPLX(-0.947111794804360643443628545664E0,-0.45987543603302557306898551016E-1);
        (*exact_roots)[8] = CMPLX(-0.947111794804360643443628545664E0,0.45987543603302557306898551016E-1);
        (*exact_roots)[9] = CMPLX(-0.935147084823051581325738268702E0,-0.276366395363089095701385358741E0);
        (*exact_roots)[10] = CMPLX(-0.935147084823051581325738268702E0,0.276366395363089095701385358741E0);
        (*exact_roots)[11] = CMPLX(-0.899451950794801943937381293773E0,-0.394568367012592113985865567971E0);
        (*exact_roots)[12] = CMPLX(-0.899451950794801943937381293773E0,0.394568367012592113985865567971E0);
        (*exact_roots)[13] = CMPLX(-0.846928410223891427281972003960E0,-0.507591758019284469647521636947E0);
        (*exact_roots)[14] = CMPLX(-0.846928410223891427281972003960E0,0.507591758019284469647521636947E0);
        (*exact_roots)[15] = CMPLX(-0.779195724062158658036386916327E0,-0.613081457400403535591975580892E0);
        (*exact_roots)[16] = CMPLX(-0.779195724062158658036386916327E0,0.613081457400403535591975580892E0);
        (*exact_roots)[17] = CMPLX(-0.697743392981121831801589116937E0,-0.709081007070975090262670227767E0);
        (*exact_roots)[18] = CMPLX(-0.697743392981121831801589116937E0,0.709081007070975090262670227767E0);
        (*exact_roots)[19] = CMPLX(-0.604121886825094827201080735988E0,-0.793887314714735149770017023538E0);
        (*exact_roots)[20] = CMPLX(-0.604121886825094827201080735988E0,0.793887314714735149770017023538E0);
        (*exact_roots)[21] = CMPLX(-0.500000000000000000000000000000E0,-0.866025403784438646763723170753E0);
        (*exact_roots)[22] = CMPLX(-0.500000000000000000000000000000E0,0.866025403784438646763723170753E0);
        (*exact_roots)[23] = CMPLX(-0.387174132037879582862879778205E0,-0.924251993602053902039550620488E0);
        (*exact_roots)[24] = CMPLX(-0.387174132037879582862879778205E0,0.924251993602053902039550620488E0);
        (*exact_roots)[25] = CMPLX(-0.267556959428167474972897159502E0,-0.967566373494130224905327435791E0);
        (*exact_roots)[26] = CMPLX(-0.267556959428167474972897159502E0,0.967566373494130224905327435791E0);
        (*exact_roots)[27] = CMPLX(-0.143155839532593341810009310219E0,-0.995221832724034632237407792227E0);
        (*exact_roots)[28] = CMPLX(-0.143155839532593341810009310219E0,0.995221832724034632237407792227E0);
        (*exact_roots)[29] = CMPLX(-0.16045517947534302681825863292E-1,-0.1006735041626386010701028433184E1);
        (*exact_roots)[30] = CMPLX(-0.16045517947534302681825863292E-1,0.1006735041626386010701028433184E1);
        (*exact_roots)[31] = CMPLX(0.111662481142841735790134563805E0,-0.1001892206008032435299692737060E1);
        (*exact_roots)[32] = CMPLX(0.111662481142841735790134563805E0,0.1001892206008032435299692737060E1);
        (*exact_roots)[33] = CMPLX(0.237852321567861477202899799722E0,-0.980751410679651646035316169672E0);
        (*exact_roots)[34] = CMPLX(0.237852321567861477202899799722E0,0.980751410679651646035316169672E0);
        (*exact_roots)[35] = CMPLX(0.360437139274099643543556447604E0,-0.943640865401314173256925688139E0);
        (*exact_roots)[36] = CMPLX(0.360437139274099643543556447604E0,0.943640865401314173256925688139E0);
        (*exact_roots)[37] = CMPLX(0.477392314844102847110474634171E0,-0.891152945345826129924640822366E0);
        (*exact_roots)[38] = CMPLX(0.477392314844102847110474634171E0,0.891152945345826129924640822366E0);
        (*exact_roots)[39] = CMPLX(0.586788023132224235580014756671E0,-0.824134044230047196220784957068E0);
        (*exact_roots)[40] = CMPLX(0.586788023132224235580014756671E0,0.824134044230047196220784957068E0);
        (*exact_roots)[41] = CMPLX(0.686820425487486753027867016412E0,-0.743670357958749198637378689329E0);
        (*exact_roots)[42] = CMPLX(0.686820425487486753027867016412E0,0.743670357958749198637378689329E0);
        (*exact_roots)[43] = CMPLX(0.775840929682842739762490529077E0,-0.651069801113127816898897251309E0);
        (*exact_roots)[44] = CMPLX(0.775840929682842739762490529077E0,0.651069801113127816898897251309E0);
        (*exact_roots)[45] = CMPLX(0.852382996720766159190958138062E0,-0.547840332639938504859225568105E0);
        (*exact_roots)[46] = CMPLX(0.852382996720766159190958138062E0,0.547840332639938504859225568105E0);
        (*exact_roots)[47] = CMPLX(0.915186027756796108167454007547E0,-0.435665032553729539118325217757E0);
        (*exact_roots)[48] = CMPLX(0.915186027756796108167454007547E0,0.435665032553729539118325217757E0);
        (*exact_roots)[49] = CMPLX(0.963215921119370289466626861950E0,-0.316374328936911136376493374596E0);
        (*exact_roots)[50] = CMPLX(0.963215921119370289466626861950E0,0.316374328936911136376493374596E0);
        (*exact_roots)[51] = CMPLX(0.995681949864471509895340077408E0,-0.191915823987060561907739223463E0);
        (*exact_roots)[52] = CMPLX(0.995681949864471509895340077408E0,0.191915823987060561907739223463E0);
        (*exact_roots)[53] = CMPLX(0.1012049674656226877398138236602E1,-0.64322209034046437066803422584E-1);
        (*exact_roots)[54] = CMPLX(0.1012049674656226877398138236602E1,0.64322209034046437066803422584E-1);
		// store Multiplicty I poly, rounded to double precision
		(*poly)[0] = CMPLX(1E0,0E0);
        (*poly)[1] = CMPLX(6E0,0E0);
        (*poly)[2] = CMPLX(15E0,0E0);
        (*poly)[3] = CMPLX(20E0,0E0);
        (*poly)[4] = CMPLX(15E0,0E0);
        (*poly)[5] = CMPLX(6E0,0E0);
        (*poly)[6] = CMPLX(1E0,0E0);
        (*poly)[7] = CMPLX(0E0,0E0);
        (*poly)[8] = CMPLX(0E0,0E0);
        (*poly)[9] = CMPLX(0E0,0E0);
        (*poly)[10] = CMPLX(0E0,0E0);
        (*poly)[11] = CMPLX(0E0,0E0);
        (*poly)[12] = CMPLX(0E0,0E0);
        (*poly)[13] = CMPLX(0E0,0E0);
        (*poly)[14] = CMPLX(0E0,0E0);
        (*poly)[15] = CMPLX(0E0,0E0);
        (*poly)[16] = CMPLX(0E0,0E0);
        (*poly)[17] = CMPLX(0E0,0E0);
        (*poly)[18] = CMPLX(0E0,0E0);
        (*poly)[19] = CMPLX(0E0,0E0);
        (*poly)[20] = CMPLX(0E0,0E0);
        (*poly)[21] = CMPLX(0E0,0E0);
        (*poly)[22] = CMPLX(0E0,0E0);
        (*poly)[23] = CMPLX(0E0,0E0);
        (*poly)[24] = CMPLX(0E0,0E0);
        (*poly)[25] = CMPLX(0E0,0E0);
        (*poly)[26] = CMPLX(0E0,0E0);
        (*poly)[27] = CMPLX(0E0,0E0);
        (*poly)[28] = CMPLX(0E0,0E0);
        (*poly)[29] = CMPLX(0E0,0E0);
        (*poly)[30] = CMPLX(0E0,0E0);
        (*poly)[31] = CMPLX(0E0,0E0);
        (*poly)[32] = CMPLX(0E0,0E0);
        (*poly)[33] = CMPLX(0E0,0E0);
        (*poly)[34] = CMPLX(0E0,0E0);
        (*poly)[35] = CMPLX(0E0,0E0);
        (*poly)[36] = CMPLX(0E0,0E0);
        (*poly)[37] = CMPLX(0E0,0E0);
        (*poly)[38] = CMPLX(0E0,0E0);
        (*poly)[39] = CMPLX(0E0,0E0);
        (*poly)[40] = CMPLX(0E0,0E0);
        (*poly)[41] = CMPLX(0E0,0E0);
        (*poly)[42] = CMPLX(0E0,0E0);
        (*poly)[43] = CMPLX(0E0,0E0);
        (*poly)[44] = CMPLX(0E0,0E0);
        (*poly)[45] = CMPLX(0E0,0E0);
        (*poly)[46] = CMPLX(0E0,0E0);
        (*poly)[47] = CMPLX(0E0,0E0);
        (*poly)[48] = CMPLX(0E0,0E0);
        (*poly)[49] = CMPLX(0E0,0E0);
        (*poly)[50] = CMPLX(1E0,0E0);
        (*poly)[51] = CMPLX(5E0,0E0);
        (*poly)[52] = CMPLX(10E0,0E0);
        (*poly)[53] = CMPLX(10E0,0E0);
        (*poly)[54] = CMPLX(5E0,0E0);
        (*poly)[55] = CMPLX(1E0,0E0);
		break;
		// Multiplicity II
		case 6:
		// set degree
		*deg = 62;
		// allocate space for poly and exact_roots
		*poly = (double complex*)malloc((*deg+1)*sizeof(double complex));
		*exact_roots = (double complex*)malloc((*deg)*sizeof(double complex));
		// store Multiplicity II exact_roots, rounded to double precision
        (*exact_roots)[0] = CMPLX(-0.998026728428271561952336806863E0,-0.62790519529313376076178224566E-1);
        (*exact_roots)[1] = CMPLX(-0.998026728428271561952336806863E0,0.62790519529313376076178224566E-1);
        (*exact_roots)[2] = CMPLX(-0.982287250728688681085641742865E0,-0.187381314585724630542550734447E0);
        (*exact_roots)[3] = CMPLX(-0.982287250728688681085641742865E0,0.187381314585724630542550734447E0);
        (*exact_roots)[4] = CMPLX(-0.951056516295153572116439333379E0,-0.309016994374947424102293417183E0);
        (*exact_roots)[5] = CMPLX(-0.951056516295153572116439333379E0,0.309016994374947424102293417183E0);
        (*exact_roots)[6] = CMPLX(-0.904827052466019527713668647933E0,-0.425779291565072648862502445744E0);
        (*exact_roots)[7] = CMPLX(-0.904827052466019527713668647933E0,0.425779291565072648862502445744E0);
        (*exact_roots)[8] = CMPLX(-0.844327925502015078548558063967E0,-0.535826794978996618271308767868E0);
        (*exact_roots)[9] = CMPLX(-0.844327925502015078548558063967E0,0.535826794978996618271308767868E0);
        (*exact_roots)[10] = CMPLX(-0.770513242775789230803009636396E0,-0.637423989748689710176712811676E0);
        (*exact_roots)[11] = CMPLX(-0.770513242775789230803009636396E0,0.637423989748689710176712811676E0);
        (*exact_roots)[12] = CMPLX(-0.684547105928688673732283357621E0,-0.728968627421411523146730319055E0);
        (*exact_roots)[13] = CMPLX(-0.684547105928688673732283357621E0,0.728968627421411523146730319055E0);
        (*exact_roots)[14] = CMPLX(-0.587785252292473129168705954639E0,-0.809016994374947424102293417183E0);
        (*exact_roots)[15] = CMPLX(-0.587785252292473129168705954639E0,0.809016994374947424102293417183E0);
        (*exact_roots)[16] = CMPLX(-0.50000000000000000000000000000E0,-0.217944947177033677611849099193E1);
        (*exact_roots)[17] = CMPLX(-0.50000000000000000000000000000E0,-0.217944947177033677611849099193E1);
        (*exact_roots)[18] = CMPLX(-0.50000000000000000000000000000E0,-0.217944947177033677611849099193E1);
        (*exact_roots)[19] = CMPLX(-0.50000000000000000000000000000E0,0.217944947177033677611849099193E1);
        (*exact_roots)[20] = CMPLX(-0.50000000000000000000000000000E0,0.217944947177033677611849099193E1);
        (*exact_roots)[21] = CMPLX(-0.50000000000000000000000000000E0,0.217944947177033677611849099193E1);
        (*exact_roots)[22] = CMPLX(-0.481753674101715274987191502872E0,-0.876306680043863587308115903922E0);
        (*exact_roots)[23] = CMPLX(-0.481753674101715274987191502872E0,0.876306680043863587308115903922E0);
        (*exact_roots)[24] = CMPLX(-0.368124552684677959156947147493E0,-0.929776485888251403660942556222E0);
        (*exact_roots)[25] = CMPLX(-0.368124552684677959156947147493E0,0.929776485888251403660942556222E0);
        (*exact_roots)[26] = CMPLX(-0.248689887164854788242283746006E0,-0.968583161128631119490168375465E0);
        (*exact_roots)[27] = CMPLX(-0.248689887164854788242283746006E0,0.968583161128631119490168375465E0);
        (*exact_roots)[28] = CMPLX(-0.125333233564304245373118759817E0,-0.992114701314477831049793042786E0);
        (*exact_roots)[29] = CMPLX(-0.125333233564304245373118759817E0,0.992114701314477831049793042786E0);
        (*exact_roots)[30] = CMPLX(0E0,-0.100000000000000000000000000000E1);
        (*exact_roots)[31] = CMPLX(0E0,0.100000000000000000000000000000E1);
        (*exact_roots)[32] = CMPLX(0.125333233564304245373118759817E0,-0.992114701314477831049793042786E0);
        (*exact_roots)[33] = CMPLX(0.125333233564304245373118759817E0,0.992114701314477831049793042786E0);
        (*exact_roots)[34] = CMPLX(0.248689887164854788242283746006E0,-0.968583161128631119490168375465E0);
        (*exact_roots)[35] = CMPLX(0.248689887164854788242283746006E0,0.968583161128631119490168375465E0);
        (*exact_roots)[36] = CMPLX(0.333333333333333333333333333333E0,0E0);
        (*exact_roots)[37] = CMPLX(0.333333333333333333333333333333E0,0E0);
        (*exact_roots)[38] = CMPLX(0.368124552684677959156947147493E0,-0.929776485888251403660942556222E0);
        (*exact_roots)[39] = CMPLX(0.368124552684677959156947147493E0,0.929776485888251403660942556222E0);
        (*exact_roots)[40] = CMPLX(0.481753674101715274987191502872E0,-0.876306680043863587308115903922E0);
        (*exact_roots)[41] = CMPLX(0.481753674101715274987191502872E0,0.876306680043863587308115903922E0);
        (*exact_roots)[42] = CMPLX(0.587785252292473129168705954639E0,-0.809016994374947424102293417183E0);
        (*exact_roots)[43] = CMPLX(0.587785252292473129168705954639E0,0.809016994374947424102293417183E0);
        (*exact_roots)[44] = CMPLX(0.684547105928688673732283357621E0,-0.728968627421411523146730319055E0);
        (*exact_roots)[45] = CMPLX(0.684547105928688673732283357621E0,0.728968627421411523146730319055E0);
        (*exact_roots)[46] = CMPLX(0.770513242775789230803009636396E0,-0.637423989748689710176712811676E0);
        (*exact_roots)[47] = CMPLX(0.770513242775789230803009636396E0,0.637423989748689710176712811676E0);
        (*exact_roots)[48] = CMPLX(0.844327925502015078548558063967E0,-0.535826794978996618271308767868E0);
        (*exact_roots)[49] = CMPLX(0.844327925502015078548558063967E0,0.535826794978996618271308767868E0);
        (*exact_roots)[50] = CMPLX(0.904827052466019527713668647933E0,-0.425779291565072648862502445744E0);
        (*exact_roots)[51] = CMPLX(0.904827052466019527713668647933E0,0.425779291565072648862502445744E0);
        (*exact_roots)[52] = CMPLX(0.951056516295153572116439333379E0,-0.309016994374947424102293417183E0);
        (*exact_roots)[53] = CMPLX(0.951056516295153572116439333379E0,0.309016994374947424102293417183E0);
        (*exact_roots)[54] = CMPLX(0.982287250728688681085641742865E0,-0.187381314585724630542550734447E0);
        (*exact_roots)[55] = CMPLX(0.982287250728688681085641742865E0,0.187381314585724630542550734447E0);
        (*exact_roots)[56] = CMPLX(0.998026728428271561952336806863E0,-0.62790519529313376076178224566E-1);
        (*exact_roots)[57] = CMPLX(0.998026728428271561952336806863E0,0.62790519529313376076178224566E-1);
        (*exact_roots)[58] = CMPLX(0.100000000000000000000000000000E1,0E0);
        (*exact_roots)[59] = CMPLX(0.100000000000000000000000000000E1,0E0);
        (*exact_roots)[60] = CMPLX(0.100000000000000000000000000000E1,0E0);
        (*exact_roots)[61] = CMPLX(0.100000000000000000000000000000E1,0E0);
		// store Multiplicity II poly, rounded to double precision
        (*poly)[0] = CMPLX(125E0,0E0);
        (*poly)[1] = CMPLX(-1175E0,0E0);
        (*poly)[2] = CMPLX(4215E0,0E0);
        (*poly)[3] = CMPLX(-7444E0,0E0);
        (*poly)[4] = CMPLX(7393E0,0E0);
        (*poly)[5] = CMPLX(-5133E0,0E0);
        (*poly)[6] = CMPLX(3402E0,0E0);
        (*poly)[7] = CMPLX(-1917E0,0E0);
        (*poly)[8] = CMPLX(741E0,0E0);
        (*poly)[9] = CMPLX(-316E0,0E0);
        (*poly)[10] = CMPLX(115E0,0E0);
        (*poly)[11] = CMPLX(-15E0,0E0);
        (*poly)[12] = CMPLX(9E0,0E0);
        (*poly)[13] = CMPLX(0E0,0E0);
        (*poly)[14] = CMPLX(0E0,0E0);
        (*poly)[15] = CMPLX(0E0,0E0);
        (*poly)[16] = CMPLX(0E0,0E0);
        (*poly)[17] = CMPLX(0E0,0E0);
        (*poly)[18] = CMPLX(0E0,0E0);
        (*poly)[19] = CMPLX(0E0,0E0);
        (*poly)[20] = CMPLX(0E0,0E0);
        (*poly)[21] = CMPLX(0E0,0E0);
        (*poly)[22] = CMPLX(0E0,0E0);
        (*poly)[23] = CMPLX(0E0,0E0);
        (*poly)[24] = CMPLX(0E0,0E0);
        (*poly)[25] = CMPLX(0E0,0E0);
        (*poly)[26] = CMPLX(0E0,0E0);
        (*poly)[27] = CMPLX(0E0,0E0);
        (*poly)[28] = CMPLX(0E0,0E0);
        (*poly)[29] = CMPLX(0E0,0E0);
        (*poly)[30] = CMPLX(0E0,0E0);
        (*poly)[31] = CMPLX(0E0,0E0);
        (*poly)[32] = CMPLX(0E0,0E0);
        (*poly)[33] = CMPLX(0E0,0E0);
        (*poly)[34] = CMPLX(0E0,0E0);
        (*poly)[35] = CMPLX(0E0,0E0);
        (*poly)[36] = CMPLX(0E0,0E0);
        (*poly)[37] = CMPLX(0E0,0E0);
        (*poly)[38] = CMPLX(0E0,0E0);
        (*poly)[39] = CMPLX(0E0,0E0);
        (*poly)[40] = CMPLX(0E0,0E0);
        (*poly)[41] = CMPLX(0E0,0E0);
        (*poly)[42] = CMPLX(0E0,0E0);
        (*poly)[43] = CMPLX(0E0,0E0);
        (*poly)[44] = CMPLX(0E0,0E0);
        (*poly)[45] = CMPLX(0E0,0E0);
        (*poly)[46] = CMPLX(0E0,0E0);
        (*poly)[47] = CMPLX(0E0,0E0);
        (*poly)[48] = CMPLX(0E0,0E0);
        (*poly)[49] = CMPLX(0E0,0E0);
        (*poly)[50] = CMPLX(125E0,0E0);
        (*poly)[51] = CMPLX(-1175E0,0E0);
        (*poly)[52] = CMPLX(4215E0,0E0);
        (*poly)[53] = CMPLX(-7444E0,0E0);
        (*poly)[54] = CMPLX(7393E0,0E0);
        (*poly)[55] = CMPLX(-5133E0,0E0);
        (*poly)[56] = CMPLX(3402E0,0E0);
        (*poly)[57] = CMPLX(-1917E0,0E0);
        (*poly)[58] = CMPLX(741E0,0E0);
        (*poly)[59] = CMPLX(-316E0,0E0);
        (*poly)[60] = CMPLX(115E0,0E0);
        (*poly)[61] = CMPLX(-15E0,0E0);
        (*poly)[62] = CMPLX(9E0,0E0);
		break;
		// Multiplicity III
		case 7:
		// set degree
		*deg = 17;
		// allocate space for poly and exact_roots
		*poly = (double complex*)malloc((*deg+1)*sizeof(double complex));
		*exact_roots = (double complex*)malloc((*deg)*sizeof(double complex));
		// store Multiplicity III exact_roots, rounded to double precision
        (*exact_roots)[0] = CMPLX(0.100000000000000000000000000000E1,0E0);
        (*exact_roots)[1] = CMPLX(0.200000000000000000000000000000E1,0E0);
        (*exact_roots)[2] = CMPLX(0.300000000000000000000000000000E1,0E0);
        (*exact_roots)[3] = CMPLX(0.400000000000000000000000000000E1,0E0);
        (*exact_roots)[4] = CMPLX(0.500000000000000000000000000000E1,0E0);
        (*exact_roots)[5] = CMPLX(0.600000000000000000000000000000E1,0E0);
        (*exact_roots)[6] = CMPLX(0.700000000000000000000000000000E1,0E0);
        (*exact_roots)[7] = CMPLX(0.800000000000000000000000000000E1,0E0);
        (*exact_roots)[8] = CMPLX(0.900000000000000000000000000000E1,0E0);
        (*exact_roots)[9] = CMPLX(0.100000000000000000000000000000E2,0E0);
        (*exact_roots)[10] = CMPLX(0.110000000000000000000000000000E2,0E0);
        (*exact_roots)[11] = CMPLX(0.120000000000000000000000000000E2,0E0);
        (*exact_roots)[12] = CMPLX(0.130000000000000000000000000000E2,0E0);
        (*exact_roots)[13] = CMPLX(0.140000000000000000000000000000E2,0E0);
        (*exact_roots)[14] = CMPLX(0.150000000000000000000000000000E2,0E0);
        (*exact_roots)[15] = CMPLX(0.150000000000000000000000000000E2,0E0);
        (*exact_roots)[16] = CMPLX(0.150000000000000000000000000000E2,0E0);
		// store Multiplicity III roots, rounded to double precision
        (*poly)[0] = CMPLX(-294226732800000E0,0E0);
        (*poly)[1] = CMPLX(1015541906400000E0,0E0);
        (*poly)[2] = CMPLX(-1518791527728000E0,0E0);
        (*poly)[3] = CMPLX(1327137724803600E0,0E0);
        (*poly)[4] = CMPLX(-766908691489440E0,0E0);
        (*poly)[5] = CMPLX(313437620164824E0,0E0);
        (*poly)[6] = CMPLX(-94377698961000E0,0E0);
        (*poly)[7] = CMPLX(21485772576905E0,0E0);
        (*poly)[8] = CMPLX(-3758453397270E0,0E0);
        (*poly)[9] = CMPLX(509681511053E0,0E0);
        (*poly)[10] = CMPLX(-53726158200E0,0E0);
        (*poly)[11] = CMPLX(4387265090E0,0E0);
        (*poly)[12] = CMPLX(-274687140E0,0E0);
        (*poly)[13] = CMPLX(12932122E0,0E0);
        (*poly)[14] = CMPLX(-442800E0,0E0);
        (*poly)[15] = CMPLX(10405E0,0E0);
        (*poly)[16] = CMPLX(-150E0,0E0);
        (*poly)[17] = CMPLX(1E0,0E0);
		break;
	}
}
#endif