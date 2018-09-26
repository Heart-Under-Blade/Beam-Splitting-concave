#include "CertainAggregate.h"
#include "Point.h"
#include "global.h"

CertainAggregate::CertainAggregate(const complex &refrIndex, double sizeIndex)
	: Particle(64, refrIndex, true)
{
	SetSymmetry(M_PI, 2*M_PI);
	SetFacetParams();
	isAggregated = true;
	{
	elems[0].origin.arr[0] = Point3f(66.11043637052107, -60.41464034745621, 18.34343028407722);
	elems[0].origin.arr[1] = Point3f(27.05498686460732, -84.3359628839182, 22.63786761597252);
	elems[0].origin.arr[2] = Point3f(-3.484149030904852, -72.88101739791193, 55.07465849713191);
	elems[0].origin.arr[3] = Point3f(5.032164579496744, -37.50474937544366, 83.21701204639601);
	elems[0].origin.arr[4] = Point3f(44.08761408541051, -13.58342683898166, 78.9225747145007);
	elems[0].origin.arr[5] = Point3f(74.62674998092268, -25.03837232498794, 46.48578383334131);
	elems[1].origin.arr[0] = Point3f(66.11043637052107, -60.41464034745621, 18.34343028407722);
	elems[1].origin.arr[1] = Point3f(74.62674998092268, -25.03837232498794, 46.48578383334131);
	elems[1].origin.arr[2] = Point3f(3.484149030904852, 72.88101739791193, -55.07465849713191);
	elems[1].origin.arr[3] = Point3f(-5.032164579496744, 37.50474937544366, -83.21701204639601);
	elems[2].origin.arr[0] = Point3f(74.62674998092268, -25.03837232498794, 46.48578383334131);
	elems[2].origin.arr[1] = Point3f(44.08761408541051, -13.58342683898166, 78.9225747145007);
	elems[2].origin.arr[2] = Point3f(-27.05498686460732, 84.3359628839182, -22.63786761597252);
	elems[2].origin.arr[3] = Point3f(3.484149030904852, 72.88101739791193, -55.07465849713191);
	elems[3].origin.arr[0] = Point3f(44.08761408541051, -13.58342683898166, 78.9225747145007);
	elems[3].origin.arr[1] = Point3f(5.032164579496744, -37.50474937544366, 83.21701204639601);
	elems[3].origin.arr[2] = Point3f(-66.11043637052107, 60.41464034745621, -18.34343028407722);
	elems[3].origin.arr[3] = Point3f(-27.05498686460732, 84.3359628839182, -22.63786761597252);
	elems[4].origin.arr[0] = Point3f(5.032164579496744, -37.50474937544366, 83.21701204639601);
	elems[4].origin.arr[1] = Point3f(-3.484149030904852, -72.88101739791193, 55.07465849713191);
	elems[4].origin.arr[2] = Point3f(-74.62674998092268, 25.03837232498794, -46.48578383334131);
	elems[4].origin.arr[3] = Point3f(-66.11043637052107, 60.41464034745621, -18.34343028407722);
	elems[5].origin.arr[0] = Point3f(-3.484149030904852, -72.88101739791193, 55.07465849713191);
	elems[5].origin.arr[1] = Point3f(27.05498686460732, -84.3359628839182, 22.63786761597252);
	elems[5].origin.arr[2] = Point3f(-44.08761408541051, 13.58342683898166, -78.9225747145007);
	elems[5].origin.arr[3] = Point3f(-74.62674998092268, 25.03837232498794, -46.48578383334131);
	elems[6].origin.arr[0] = Point3f(27.05498686460732, -84.3359628839182, 22.63786761597252);
	elems[6].origin.arr[1] = Point3f(66.11043637052107, -60.41464034745621, 18.34343028407722);
	elems[6].origin.arr[2] = Point3f(-5.032164579496744, 37.50474937544366, -83.21701204639601);
	elems[6].origin.arr[3] = Point3f(-44.08761408541051, 13.58342683898166, -78.9225747145007);
	elems[7].origin.arr[0] = Point3f(-5.032164579496744, 37.50474937544366, -83.21701204639601);
	elems[7].origin.arr[1] = Point3f(3.484149030904852, 72.88101739791193, -55.07465849713191);
	elems[7].origin.arr[2] = Point3f(-27.05498686460732, 84.3359628839182, -22.63786761597252);
	elems[7].origin.arr[3] = Point3f(-66.11043637052107, 60.41464034745621, -18.34343028407722);
	elems[7].origin.arr[4] = Point3f(-74.62674998092268, 25.03837232498794, -46.48578383334131);
	elems[7].origin.arr[5] = Point3f(-44.08761408541051, 13.58342683898166, -78.9225747145007);


	elems[8].origin.arr[0] = Point3f(-50.11391689844849, 124.4704224181641, -88.38614194106512);
	elems[8].origin.arr[1] = Point3f(-32.94480156339847, 159.3110758840723, -78.82840451632649);
	elems[8].origin.arr[2] = Point3f(-22.96537288621924, 166.9368539495522, -40.85132574276704);
	elems[8].origin.arr[3] = Point3f(-30.15505954409002, 139.7219785491239, -12.43198439394624);
	elems[8].origin.arr[4] = Point3f(-47.32417487914003, 104.8813250832157, -21.98972181868488);
	elems[8].origin.arr[5] = Point3f(-57.30360355631927, 97.25554701773579, -59.96680059224432);
	elems[9].origin.arr[0] = Point3f(-50.11391689844849, 124.4704224181641, -88.38614194106512);
	elems[9].origin.arr[1] = Point3f(-57.30360355631927, 97.25554701773579, -59.96680059224432);
	elems[9].origin.arr[2] = Point3f(54.58137288621924, 47.4411460504478, -79.36467425723296);
	elems[9].origin.arr[3] = Point3f(61.77105954409002, 74.65602145087612, -107.7840156060538);
	elems[10].origin.arr[0] = Point3f(-57.30360355631927, 97.25554701773579, -59.96680059224432);
	elems[10].origin.arr[1] = Point3f(-47.32417487914003, 104.8813250832157, -21.98972181868488);
	elems[10].origin.arr[2] = Point3f(64.56080156339847, 55.06692411592767, -41.38759548367351);
	elems[10].origin.arr[3] = Point3f(54.58137288621924, 47.4411460504478, -79.36467425723296);
	elems[11].origin.arr[0] = Point3f(-47.32417487914003, 104.8813250832157, -21.98972181868488);
	elems[11].origin.arr[1] = Point3f(-30.15505954409002, 139.7219785491239, -12.43198439394624);
	elems[11].origin.arr[2] = Point3f(81.72991689844849, 89.90757758183587, -31.82985805893487);
	elems[11].origin.arr[3] = Point3f(64.56080156339847, 55.06692411592767, -41.38759548367351);
	elems[12].origin.arr[0] = Point3f(-30.15505954409002, 139.7219785491239, -12.43198439394624);
	elems[12].origin.arr[1] = Point3f(-22.96537288621924, 166.9368539495522, -40.85132574276704);
	elems[12].origin.arr[2] = Point3f(88.91960355631926, 117.1224529822642, -60.24919940775568);
	elems[12].origin.arr[3] = Point3f(81.72991689844849, 89.90757758183587, -31.82985805893487);
	elems[13].origin.arr[0] = Point3f(-22.96537288621924, 166.9368539495522, -40.85132574276704);
	elems[13].origin.arr[1] = Point3f(-32.94480156339847, 159.3110758840723, -78.82840451632649);
	elems[13].origin.arr[2] = Point3f(78.94017487914003, 109.4966749167843, -98.22627818131511);
	elems[13].origin.arr[3] = Point3f(88.91960355631926, 117.1224529822642, -60.24919940775568);
	elems[14].origin.arr[0] = Point3f(-32.94480156339847, 159.3110758840723, -78.82840451632649);
	elems[14].origin.arr[1] = Point3f(-50.11391689844849, 124.4704224181641, -88.38614194106512);
	elems[14].origin.arr[2] = Point3f(61.77105954409002, 74.65602145087612, -107.7840156060538);
	elems[14].origin.arr[3] = Point3f(78.94017487914003, 109.4966749167843, -98.22627818131511);
	elems[15].origin.arr[0] = Point3f(61.77105954409002, 74.65602145087612, -107.7840156060538);
	elems[15].origin.arr[1] = Point3f(54.58137288621924, 47.4411460504478, -79.36467425723296);
	elems[15].origin.arr[2] = Point3f(64.56080156339847, 55.06692411592767, -41.38759548367351);
	elems[15].origin.arr[3] = Point3f(81.72991689844849, 89.90757758183587, -31.82985805893487);
	elems[15].origin.arr[4] = Point3f(88.91960355631926, 117.1224529822642, -60.24919940775568);
	elems[15].origin.arr[5] = Point3f(78.94017487914003, 109.4966749167843, -98.22627818131511);


	elems[16].origin.arr[0] = Point3f(-32.46675231729333, 120.6181016912531, 46.84750551506697);
	elems[16].origin.arr[1] = Point3f(-6.702378116253382, 115.9591012296608, 56.77175544781602);
	elems[16].origin.arr[2] = Point3f(-3.208231155621423, 100.9744762735481, 80.16517229833511);
	elems[16].origin.arr[3] = Point3f(-25.47845839602941, 90.64885177902769, 93.63433921610513);
	elems[16].origin.arr[4] = Point3f(-51.24283259706935, 95.30785224062005, 83.7100892833561);
	elems[16].origin.arr[5] = Point3f(-54.73697955770132, 110.2924771967328, 60.31667243283701);
	elems[17].origin.arr[0] = Point3f(-32.46675231729333, 120.6181016912531, 46.84750551506697);
	elems[17].origin.arr[1] = Point3f(-54.73697955770132, 110.2924771967328, 60.31667243283701);
	elems[17].origin.arr[2] = Point3f(-50.17376884437857, 45.03552372645191, 17.8348277016649);
	elems[17].origin.arr[3] = Point3f(-27.90354160397059, 55.36114822097229, 4.36566078389486);
	elems[18].origin.arr[0] = Point3f(-54.73697955770132, 110.2924771967328, 60.31667243283701);
	elems[18].origin.arr[1] = Point3f(-51.24283259706935, 95.30785224062005, 83.7100892833561);
	elems[18].origin.arr[2] = Point3f(-46.67962188374662, 30.0508987703392, 41.22824455218398);
	elems[18].origin.arr[3] = Point3f(-50.17376884437857, 45.03552372645191, 17.8348277016649);
	elems[19].origin.arr[0] = Point3f(-51.24283259706935, 95.30785224062005, 83.7100892833561);
	elems[19].origin.arr[1] = Point3f(-25.47845839602941, 90.64885177902769, 93.63433921610513);
	elems[19].origin.arr[2] = Point3f(-20.91524768270667, 25.39189830874685, 51.15249448493303);
	elems[19].origin.arr[3] = Point3f(-46.67962188374662, 30.0508987703392, 41.22824455218398);
	elems[20].origin.arr[0] = Point3f(-25.47845839602941, 90.64885177902769, 93.63433921610513);
	elems[20].origin.arr[1] = Point3f(-3.208231155621423, 100.9744762735481, 80.16517229833511);
	elems[20].origin.arr[2] = Point3f(1.354979557701316, 35.71752280326722, 37.68332756716299);
	elems[20].origin.arr[3] = Point3f(-20.91524768270667, 25.39189830874685, 51.15249448493303);
	elems[21].origin.arr[0] = Point3f(-3.208231155621423, 100.9744762735481, 80.16517229833511);
	elems[21].origin.arr[1] = Point3f(-6.702378116253382, 115.9591012296608, 56.77175544781602);
	elems[21].origin.arr[2] = Point3f(-2.139167402930642, 50.70214775937994, 14.28991071664391);
	elems[21].origin.arr[3] = Point3f(1.354979557701316, 35.71752280326722, 37.68332756716299);
	elems[22].origin.arr[0] = Point3f(-6.702378116253382, 115.9591012296608, 56.77175544781602);
	elems[22].origin.arr[1] = Point3f(-32.46675231729333, 120.6181016912531, 46.84750551506697);
	elems[22].origin.arr[2] = Point3f(-27.90354160397059, 55.36114822097229, 4.36566078389486);
	elems[22].origin.arr[3] = Point3f(-2.139167402930642, 50.70214775937994, 14.28991071664391);
	elems[23].origin.arr[0] = Point3f(-27.90354160397059, 55.36114822097229, 4.36566078389486);
	elems[23].origin.arr[1] = Point3f(-50.17376884437857, 45.03552372645191, 17.8348277016649);
	elems[23].origin.arr[2] = Point3f(-46.67962188374662, 30.0508987703392, 41.22824455218398);
	elems[23].origin.arr[3] = Point3f(-20.91524768270667, 25.39189830874685, 51.15249448493303);
	elems[23].origin.arr[4] = Point3f(1.354979557701316, 35.71752280326722, 37.68332756716299);
	elems[23].origin.arr[5] = Point3f(-2.139167402930642, 50.70214775937994, 14.28991071664391);


	elems[24].origin.arr[0] = Point3f(-142.8371105044146, 9.364161125355935, -41.78242184149515);
	elems[24].origin.arr[1] = Point3f(-105.4913330702954, 36.26900924041865, -28.1654450938232);
	elems[24].origin.arr[2] = Point3f(-89.94694901254758, 34.54209602069939, 17.21505617045101);
	elems[24].origin.arr[3] = Point3f(-111.748342388919, 5.910334685917398, 48.97858068705328);
	elems[24].origin.arr[4] = Point3f(-149.0941198230382, -20.99451342914533, 35.36160393938133);
	elems[24].origin.arr[5] = Point3f(-164.638503880786, -19.26760020942605, -10.01889732489289);
	elems[25].origin.arr[0] = Point3f(-142.8371105044146, 9.364161125355935, -41.78242184149515);
	elems[25].origin.arr[1] = Point3f(-164.638503880786, -19.26760020942605, -10.01889732489289);
	elems[25].origin.arr[2] = Point3f(-86.05305098745242, -112.9220960206994, -40.50105617045102);
	elems[25].origin.arr[3] = Point3f(-64.25165761108104, -84.2903346859174, -72.26458068705328);
	elems[26].origin.arr[0] = Point3f(-164.638503880786, -19.26760020942605, -10.01889732489289);
	elems[26].origin.arr[1] = Point3f(-149.0941198230382, -20.99451342914533, 35.36160393938133);
	elems[26].origin.arr[2] = Point3f(-70.50866692970459, -114.6490092404186, 4.8794450938232);
	elems[26].origin.arr[3] = Point3f(-86.05305098745242, -112.9220960206994, -40.50105617045102);
	elems[27].origin.arr[0] = Point3f(-149.0941198230382, -20.99451342914533, 35.36160393938133);
	elems[27].origin.arr[1] = Point3f(-111.748342388919, 5.910334685917398, 48.97858068705328);
	elems[27].origin.arr[2] = Point3f(-33.16288949558539, -87.74416112535593, 18.49642184149515);
	elems[27].origin.arr[3] = Point3f(-70.50866692970459, -114.6490092404186, 4.8794450938232);
	elems[28].origin.arr[0] = Point3f(-111.748342388919, 5.910334685917398, 48.97858068705328);
	elems[28].origin.arr[1] = Point3f(-89.94694901254758, 34.54209602069939, 17.21505617045101);
	elems[28].origin.arr[2] = Point3f(-11.36149611921401, -59.11239979057395, -13.26710267510711);
	elems[28].origin.arr[3] = Point3f(-33.16288949558539, -87.74416112535593, 18.49642184149515);
	elems[29].origin.arr[0] = Point3f(-89.94694901254758, 34.54209602069939, 17.21505617045101);
	elems[29].origin.arr[1] = Point3f(-105.4913330702954, 36.26900924041865, -28.1654450938232);
	elems[29].origin.arr[2] = Point3f(-26.90588017696184, -57.38548657085467, -58.64760393938133);
	elems[29].origin.arr[3] = Point3f(-11.36149611921401, -59.11239979057395, -13.26710267510711);
	elems[30].origin.arr[0] = Point3f(-105.4913330702954, 36.26900924041865, -28.1654450938232);
	elems[30].origin.arr[1] = Point3f(-142.8371105044146, 9.364161125355935, -41.78242184149515);
	elems[30].origin.arr[2] = Point3f(-64.25165761108104, -84.2903346859174, -72.26458068705328);
	elems[30].origin.arr[3] = Point3f(-26.90588017696184, -57.38548657085467, -58.64760393938133);
	elems[31].origin.arr[0] = Point3f(-64.25165761108104, -84.2903346859174, -72.26458068705328);
	elems[31].origin.arr[1] = Point3f(-86.05305098745242, -112.9220960206994, -40.50105617045102);
	elems[31].origin.arr[2] = Point3f(-70.50866692970459, -114.6490092404186, 4.8794450938232);
	elems[31].origin.arr[3] = Point3f(-33.16288949558539, -87.74416112535593, 18.49642184149515);
	elems[31].origin.arr[4] = Point3f(-11.36149611921401, -59.11239979057395, -13.26710267510711);
	elems[31].origin.arr[5] = Point3f(-26.90588017696184, -57.38548657085467, -58.64760393938133);


	elems[32].origin.arr[0] = Point3f(185.2248364132817, 13.70500759335614, 65.55079690640709);
	elems[32].origin.arr[1] = Point3f(153.1768486251542, -27.67065939766073, 73.91624129272087);
	elems[32].origin.arr[2] = Point3f(107.0718437987181, -20.80496351699799, 99.13906330035027);
	elems[32].origin.arr[3] = Point3f(93.01482676040946, 27.43639935468163, 115.9964409216659);
	elems[32].origin.arr[4] = Point3f(125.0628145485369, 68.8120663456985, 107.6309965353521);
	elems[32].origin.arr[5] = Point3f(171.167819374973, 61.94637046503576, 82.40817452772272);
	elems[33].origin.arr[0] = Point3f(185.2248364132817, 13.70500759335614, 65.55079690640709);
	elems[33].origin.arr[1] = Point3f(171.167819374973, 61.94637046503576, 82.40817452772272);
	elems[33].origin.arr[2] = Point3f(105.9921562012819, 86.96496351699798, -43.53706330035027);
	elems[33].origin.arr[3] = Point3f(120.0491732395905, 38.72360064531837, -60.39444092166589);
	elems[34].origin.arr[0] = Point3f(171.167819374973, 61.94637046503576, 82.40817452772272);
	elems[34].origin.arr[1] = Point3f(125.0628145485369, 68.8120663456985, 107.6309965353521);
	elems[34].origin.arr[2] = Point3f(59.88715137484581, 93.83065939766072, -18.31424129272087);
	elems[34].origin.arr[3] = Point3f(105.9921562012819, 86.96496351699798, -43.53706330035027);
	elems[35].origin.arr[0] = Point3f(125.0628145485369, 68.8120663456985, 107.6309965353521);
	elems[35].origin.arr[1] = Point3f(93.01482676040946, 27.43639935468163, 115.9964409216659);
	elems[35].origin.arr[2] = Point3f(27.83916358671833, 52.45499240664385, -9.948796906407093);
	elems[35].origin.arr[3] = Point3f(59.88715137484581, 93.83065939766072, -18.31424129272087);
	elems[36].origin.arr[0] = Point3f(93.01482676040946, 27.43639935468163, 115.9964409216659);
	elems[36].origin.arr[1] = Point3f(107.0718437987181, -20.80496351699799, 99.13906330035027);
	elems[36].origin.arr[2] = Point3f(41.89618062502694, 4.213629534964245, -26.80617452772272);
	elems[36].origin.arr[3] = Point3f(27.83916358671833, 52.45499240664385, -9.948796906407093);
	elems[37].origin.arr[0] = Point3f(107.0718437987181, -20.80496351699799, 99.13906330035027);
	elems[37].origin.arr[1] = Point3f(153.1768486251542, -27.67065939766073, 73.91624129272087);
	elems[37].origin.arr[2] = Point3f(88.00118545146304, -2.652066345698501, -52.02899653535212);
	elems[37].origin.arr[3] = Point3f(41.89618062502694, 4.213629534964245, -26.80617452772272);
	elems[38].origin.arr[0] = Point3f(153.1768486251542, -27.67065939766073, 73.91624129272087);
	elems[38].origin.arr[1] = Point3f(185.2248364132817, 13.70500759335614, 65.55079690640709);
	elems[38].origin.arr[2] = Point3f(120.0491732395905, 38.72360064531837, -60.39444092166589);
	elems[38].origin.arr[3] = Point3f(88.00118545146304, -2.652066345698501, -52.02899653535212);
	elems[39].origin.arr[0] = Point3f(120.0491732395905, 38.72360064531837, -60.39444092166589);
	elems[39].origin.arr[1] = Point3f(105.9921562012819, 86.96496351699798, -43.53706330035027);
	elems[39].origin.arr[2] = Point3f(59.88715137484581, 93.83065939766072, -18.31424129272087);
	elems[39].origin.arr[3] = Point3f(27.83916358671833, 52.45499240664385, -9.948796906407093);
	elems[39].origin.arr[4] = Point3f(41.89618062502694, 4.213629534964245, -26.80617452772272);
	elems[39].origin.arr[5] = Point3f(88.00118545146304, -2.652066345698501, -52.02899653535212);


	elems[40].origin.arr[0] = Point3f(5.244814931859345, -63.04768028164555, -41.47000914880908);
	elems[40].origin.arr[1] = Point3f(3.601083929062028, -45.19258458605589, -35.18560734022149);
	elems[40].origin.arr[2] = Point3f(11.36318908272452, -40.21598447956706, -18.57286599619336);
	elems[40].origin.arr[3] = Point3f(20.76902523918433, -53.09448006866788, -8.244526460752819);
	elems[40].origin.arr[4] = Point3f(22.41275624198164, -70.94957576425753, -14.5289282693404);
	elems[40].origin.arr[5] = Point3f(14.65065108831915, -75.92617587074636, -31.14166961336853);
	elems[41].origin.arr[0] = Point3f(5.244814931859345, -63.04768028164555, -41.47000914880908);
	elems[41].origin.arr[1] = Point3f(14.65065108831915, -75.92617587074636, -31.14166961336853);
	elems[41].origin.arr[2] = Point3f(60.48281091727549, -62.78401552043294, -56.49313400380664);
	elems[41].origin.arr[3] = Point3f(51.07697476081567, -49.90551993133212, -66.82147353924718);
	elems[42].origin.arr[0] = Point3f(14.65065108831915, -75.92617587074636, -31.14166961336853);
	elems[42].origin.arr[1] = Point3f(22.41275624198164, -70.94957576425753, -14.5289282693404);
	elems[42].origin.arr[2] = Point3f(68.24491607093798, -57.80741541394411, -39.88039265977851);
	elems[42].origin.arr[3] = Point3f(60.48281091727549, -62.78401552043294, -56.49313400380664);
	elems[43].origin.arr[0] = Point3f(22.41275624198164, -70.94957576425753, -14.5289282693404);
	elems[43].origin.arr[1] = Point3f(20.76902523918433, -53.09448006866788, -8.244526460752819);
	elems[43].origin.arr[2] = Point3f(66.60118506814067, -39.95231971835445, -33.59599085119093);
	elems[43].origin.arr[3] = Point3f(68.24491607093798, -57.80741541394411, -39.88039265977851);
	elems[44].origin.arr[0] = Point3f(20.76902523918433, -53.09448006866788, -8.244526460752819);
	elems[44].origin.arr[1] = Point3f(11.36318908272452, -40.21598447956706, -18.57286599619336);
	elems[44].origin.arr[2] = Point3f(57.19534891168085, -27.07382412925364, -43.92433038663147);
	elems[44].origin.arr[3] = Point3f(66.60118506814067, -39.95231971835445, -33.59599085119093);
	elems[45].origin.arr[0] = Point3f(11.36318908272452, -40.21598447956706, -18.57286599619336);
	elems[45].origin.arr[1] = Point3f(3.601083929062028, -45.19258458605589, -35.18560734022149);
	elems[45].origin.arr[2] = Point3f(49.43324375801836, -32.05042423574247, -60.5370717306596);
	elems[45].origin.arr[3] = Point3f(57.19534891168085, -27.07382412925364, -43.92433038663147);
	elems[46].origin.arr[0] = Point3f(3.601083929062028, -45.19258458605589, -35.18560734022149);
	elems[46].origin.arr[1] = Point3f(5.244814931859345, -63.04768028164555, -41.47000914880908);
	elems[46].origin.arr[2] = Point3f(51.07697476081567, -49.90551993133212, -66.82147353924718);
	elems[46].origin.arr[3] = Point3f(49.43324375801836, -32.05042423574247, -60.5370717306596);
	elems[47].origin.arr[0] = Point3f(51.07697476081567, -49.90551993133212, -66.82147353924718);
	elems[47].origin.arr[1] = Point3f(60.48281091727549, -62.78401552043294, -56.49313400380664);
	elems[47].origin.arr[2] = Point3f(68.24491607093798, -57.80741541394411, -39.88039265977851);
	elems[47].origin.arr[3] = Point3f(66.60118506814067, -39.95231971835445, -33.59599085119093);
	elems[47].origin.arr[4] = Point3f(57.19534891168085, -27.07382412925364, -43.92433038663147);
	elems[47].origin.arr[5] = Point3f(49.43324375801836, -32.05042423574247, -60.5370717306596);


	elems[48].origin.arr[0] = Point3f(53.78577786126138, -0.572803734518736, 131.480915294099);
	elems[48].origin.arr[1] = Point3f(83.0022415228523, -17.95804752578802, 131.8702087940059);
	elems[48].origin.arr[2] = Point3f(86.05596890084885, -45.63589073138481, 151.3794820912683);
	elems[48].origin.arr[3] = Point3f(59.8932326172545, -55.92849014571232, 170.4994618886237);
	elems[48].origin.arr[4] = Point3f(30.67676895566357, -38.54324635444303, 170.1101683887168);
	elems[48].origin.arr[5] = Point3f(27.62304157766701, -10.86540314884625, 150.6008950914544);
	elems[49].origin.arr[0] = Point3f(53.78577786126138, -0.572803734518736, 131.480915294099);
	elems[49].origin.arr[1] = Point3f(27.62304157766701, -10.86540314884625, 150.6008950914544);
	elems[49].origin.arr[2] = Point3f(-5.835968900848862, -68.81810926861519, 73.6205179087317);
	elems[49].origin.arr[3] = Point3f(20.3267673827455, -58.52550985428768, 54.50053811137628);
	elems[50].origin.arr[0] = Point3f(27.62304157766701, -10.86540314884625, 150.6008950914544);
	elems[50].origin.arr[1] = Point3f(30.67676895566357, -38.54324635444303, 170.1101683887168);
	elems[50].origin.arr[2] = Point3f(-2.78224152285231, -96.49595247421198, 93.12979120599405);
	elems[50].origin.arr[3] = Point3f(-5.835968900848862, -68.81810926861519, 73.6205179087317);
	elems[51].origin.arr[0] = Point3f(30.67676895566357, -38.54324635444303, 170.1101683887168);
	elems[51].origin.arr[1] = Point3f(59.8932326172545, -55.92849014571232, 170.4994618886237);
	elems[51].origin.arr[2] = Point3f(26.43422213873862, -113.8811962654813, 93.51908470590098);
	elems[51].origin.arr[3] = Point3f(-2.78224152285231, -96.49595247421198, 93.12979120599405);
	elems[52].origin.arr[0] = Point3f(59.8932326172545, -55.92849014571232, 170.4994618886237);
	elems[52].origin.arr[1] = Point3f(86.05596890084885, -45.63589073138481, 151.3794820912683);
	elems[52].origin.arr[2] = Point3f(52.59695842233299, -103.5885968511537, 74.39910490854555);
	elems[52].origin.arr[3] = Point3f(26.43422213873862, -113.8811962654813, 93.51908470590098);
	elems[53].origin.arr[0] = Point3f(86.05596890084885, -45.63589073138481, 151.3794820912683);
	elems[53].origin.arr[1] = Point3f(83.0022415228523, -17.95804752578802, 131.8702087940059);
	elems[53].origin.arr[2] = Point3f(49.54323104433643, -75.91075364555695, 54.8898316112832);
	elems[53].origin.arr[3] = Point3f(52.59695842233299, -103.5885968511537, 74.39910490854555);
	elems[54].origin.arr[0] = Point3f(83.0022415228523, -17.95804752578802, 131.8702087940059);
	elems[54].origin.arr[1] = Point3f(53.78577786126138, -0.572803734518736, 131.480915294099);
	elems[54].origin.arr[2] = Point3f(20.3267673827455, -58.52550985428768, 54.50053811137628);
	elems[54].origin.arr[3] = Point3f(49.54323104433643, -75.91075364555695, 54.8898316112832);
	elems[55].origin.arr[0] = Point3f(20.3267673827455, -58.52550985428768, 54.50053811137628);
	elems[55].origin.arr[1] = Point3f(-5.835968900848862, -68.81810926861519, 73.6205179087317);
	elems[55].origin.arr[2] = Point3f(-2.78224152285231, -96.49595247421198, 93.12979120599405);
	elems[55].origin.arr[3] = Point3f(26.43422213873862, -113.8811962654813, 93.51908470590098);
	elems[55].origin.arr[4] = Point3f(52.59695842233299, -103.5885968511537, 74.39910490854555);
	elems[55].origin.arr[5] = Point3f(49.54323104433643, -75.91075364555695, 54.8898316112832);


	elems[56].origin.arr[0] = Point3f(-31.99942627699921, -194.0207287020714, 104.7597626660042);
	elems[56].origin.arr[1] = Point3f(-63.79327037391881, -165.2479307279772, 107.9656282718477);
	elems[56].origin.arr[2] = Point3f(-55.83310479334075, -126.090958509988, 123.8517004940618);
	elems[56].origin.arr[3] = Point3f(-16.0790951158431, -115.7067842660931, 136.5319071104325);
	elems[56].origin.arr[4] = Point3f(15.71474898107649, -144.4795822401873, 133.3260415045891);
	elems[56].origin.arr[5] = Point3f(7.754583400498438, -183.6365544581765, 117.439969282375);
	elems[57].origin.arr[0] = Point3f(-31.99942627699921, -194.0207287020714, 104.7597626660042);
	elems[57].origin.arr[1] = Point3f(7.754583400498438, -183.6365544581765, 117.439969282375);
	elems[57].origin.arr[2] = Point3f(36.32830479334075, -137.909041490012, -9.589700494061816);
	elems[57].origin.arr[3] = Point3f(-3.425704884156898, -148.2932157339069, -22.26990711043253);
	elems[58].origin.arr[0] = Point3f(7.754583400498438, -183.6365544581765, 117.439969282375);
	elems[58].origin.arr[1] = Point3f(15.71474898107649, -144.4795822401873, 133.3260415045891);
	elems[58].origin.arr[2] = Point3f(44.28847037391881, -98.75206927202279, 6.296371728152337);
	elems[58].origin.arr[3] = Point3f(36.32830479334075, -137.909041490012, -9.589700494061816);
	elems[59].origin.arr[0] = Point3f(15.71474898107649, -144.4795822401873, 133.3260415045891);
	elems[59].origin.arr[1] = Point3f(-16.0790951158431, -115.7067842660931, 136.5319071104325);
	elems[59].origin.arr[2] = Point3f(12.49462627699921, -69.97927129792856, 9.502237333995758);
	elems[59].origin.arr[3] = Point3f(44.28847037391881, -98.75206927202279, 6.296371728152337);
	elems[60].origin.arr[0] = Point3f(-16.0790951158431, -115.7067842660931, 136.5319071104325);
	elems[60].origin.arr[1] = Point3f(-55.83310479334075, -126.090958509988, 123.8517004940618);
	elems[60].origin.arr[2] = Point3f(-27.25938340049844, -80.36344554182351, -3.17796928237496);
	elems[60].origin.arr[3] = Point3f(12.49462627699921, -69.97927129792856, 9.502237333995758);
	elems[61].origin.arr[0] = Point3f(-55.83310479334075, -126.090958509988, 123.8517004940618);
	elems[61].origin.arr[1] = Point3f(-63.79327037391881, -165.2479307279772, 107.9656282718477);
	elems[61].origin.arr[2] = Point3f(-35.21954898107649, -119.5204177598127, -19.06404150458911);
	elems[61].origin.arr[3] = Point3f(-27.25938340049844, -80.36344554182351, -3.17796928237496);
	elems[62].origin.arr[0] = Point3f(-63.79327037391881, -165.2479307279772, 107.9656282718477);
	elems[62].origin.arr[1] = Point3f(-31.99942627699921, -194.0207287020714, 104.7597626660042);
	elems[62].origin.arr[2] = Point3f(-3.425704884156898, -148.2932157339069, -22.26990711043253);
	elems[62].origin.arr[3] = Point3f(-35.21954898107649, -119.5204177598127, -19.06404150458911);
	elems[63].origin.arr[0] = Point3f(-3.425704884156898, -148.2932157339069, -22.26990711043253);
	elems[63].origin.arr[1] = Point3f(36.32830479334075, -137.909041490012, -9.589700494061816);
	elems[63].origin.arr[2] = Point3f(44.28847037391881, -98.75206927202279, 6.296371728152337);
	elems[63].origin.arr[3] = Point3f(12.49462627699921, -69.97927129792856, 9.502237333995758);
	elems[63].origin.arr[4] = Point3f(-27.25938340049844, -80.36344554182351, -3.17796928237496);
	elems[63].origin.arr[5] = Point3f(-35.21954898107649, -119.5204177598127, -19.06404150458911);
	}

	Resize(sizeIndex);

	for (int i = 0; i < nElems; ++i)
	{
		int last = elems[i].origin.nVertices-1;

		for (int j = 0; j < elems[i].origin.nVertices/2; ++j)
		{
			Point3f buf = elems[i].origin.arr[j];
			elems[i].origin.arr[j] = elems[i].origin.arr[last-j];
			elems[i].origin.arr[last-j] = buf;
		}
	}

	SetDefaultNormals();
	SetDefaultCenters();
	Reset();

	for (int i = 0; i < nElems; ++i)
	{
		elems[i].origin.isOverlayedIn = false;
		elems[i].origin.isOverlayedOut = false;
		elems[i].actual.isOverlayedIn = false;
		elems[i].actual.isOverlayedOut = false;
	}
}

void CertainAggregate::Resize(double sizeIndex)
{
	for (int i = 0; i < nElems; ++i)
	{
		for (int j = 0; j < elems[i].origin.nVertices; ++j)
		{
			elems[i].origin.arr[j].cx *= sizeIndex;
			elems[i].origin.arr[j].cy *= sizeIndex;
			elems[i].origin.arr[j].cz *= sizeIndex;
		}
	}
}

void CertainAggregate::SetFacetParams()
{
	for (int i = 0; i < nElems; ++i)
	{
		if (i%8 == 7 || i%8 == 0 || i == 0)
		{
			elems[i].origin.nVertices = 6;
		}
		else
		{
			elems[i].origin.nVertices = 4;
		}
	}

}

void CertainAggregate::GetParticalFacetIdRange(Facet *facet, int &begin, int &end) const
{
	// REF: упростить
	int &index = facet->index;

	if (index < 8)
	{
		begin = 0;
		end = 8;
	}
	else if (index >= 8 && index < 16)
	{
		begin = 8;
		end = 16;
	}
	else if (index >= 16 && index < 24)
	{
		begin = 16;
		end = 24;
	}
	else if (index >= 24 && index < 32)
	{
		begin = 24;
		end = 32;
	}
	else if (index >= 32 && index < 40)
	{
		begin = 32;
		end = 40;
	}
	else if (index >= 40 && index < 48)
	{
		begin = 40;
		end = 48;
	}
	else if (index >= 48 && index < 56)
	{
		begin = 48;
		end = 56;
	}
	else if (index >= 56 && index < 64)
	{
		begin = 56;
		end = 64;
	}
}
