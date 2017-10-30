
objects= anaMuon.exe NuLatAB-Plot.exe NuLatMuon-Plot.exe NuLat-Primary.exe NuLatAB.exe NuLatAB-PSD.exe anaposAB.exe NuLatMuon.exe NuLatABNoRoot.exe  shell-generation.exe Map.exe PSD-Shell.exe event-list.exe NuLatAB-pulse.exe
analysis: $(objects)

anaMuon.exe: anaMuon.cpp anafunction.h
	g++ anaMuon.cpp -L/home/xinjian/program/root/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic -pthread -std=c++11 -m64 -I/home/xinjian/program/root/include -o anaMuon.exe

NuLatAB-Plot.exe: NuLatAB-Plot.cpp anafunction.h
	g++ NuLatAB-Plot.cpp -L/home/xinjian/program/root/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic -pthread -std=c++11 -m64 -I/home/xinjian/program/root/include -o NuLatAB-Plot.exe

NuLatMuon-Plot.exe: NuLatMuon-Plot.cpp anafunction.h
	g++ NuLatMuon-Plot.cpp -L/home/xinjian/program/root/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic -pthread -std=c++11 -m64 -I/home/xinjian/program/root/include -o NuLatMuon-Plot.exe


NuLat-Primary.exe: NuLat-Primary.cpp anafunction.h
	g++ NuLat-Primary.cpp -L/home/xinjian/program/root/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic -pthread -std=c++11 -m64 -I/home/xinjian/program/root/include -o NuLat-Primary.exe

NuLatAB.exe: NuLatAB.cpp anafunction.h
	g++ NuLatAB.cpp -L/home/xinjian/program/root/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic -pthread -std=c++11 -m64 -I/home/xinjian/program/root/include -o NuLatAB.exe

NuLatAB-PSD.exe: NuLatAB-PSD.cpp anafunction.h
	g++ NuLatAB-PSD.cpp -L/home/xinjian/program/root/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic -pthread -std=c++11 -m64 -I/home/xinjian/program/root/include -o NuLatAB-PSD.exe

anaposAB.exe: anaposAB.cpp anafunction.h
	g++ anaposAB.cpp -L/home/xinjian/program/root/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic -pthread -std=c++11 -m64 -I/home/xinjian/program/root/include -o anaposAB.exe

NuLatMuon.exe: NuLatMuon.cpp anafunction.h
	g++ NuLatMuon.cpp -L/home/xinjian/program/root/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic -pthread -std=c++11 -m64 -I/home/xinjian/program/root/include -o NuLatMuon.exe

NuLatABNoRoot.exe: NuLatAB-NoRoot.cpp anafunction.h
	g++ NuLatAB-NoRoot.cpp -L/home/xinjian/program/root/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic -pthread -std=c++11 -m64 -I/home/xinjian/program/root/include -o NuLatABNoRoot.exe

NuLatAB-pulse.exe: pulse-shape-output.cpp anafunction.h
	g++ pulse-shape-output.cpp -L/home/xinjian/program/root/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic -pthread -std=c++11 -m64 -I/home/xinjian/program/root/include -o NuLatAB-pulse.exe


#NuLat-Window.exe: NuLat-Window.cpp anafunction.h
#l	g++ NuLat-Window.cpp -L/home/xinjian/program/root/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic -pthread -std=c++11 -m64 -I/home/xinjian/program/root/include -o NuLat-Window.exe


shell-generation.exe: shell.cpp
	g++ shell.cpp -o shell-generation.exe

Map.exe: NuLatMap.cpp
	g++ NuLatMap.cpp -o Map.exe

PSD-Shell.exe: Shell-PSD.cpp
	g++ Shell-PSD.cpp -o PSD-Shell.exe

event-list.exe: Event-All-list.cpp
	g++ Event-All-list.cpp -o event-list.exe



