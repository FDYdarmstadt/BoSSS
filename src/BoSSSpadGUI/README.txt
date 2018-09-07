############################################################################################
Run from npm: 
############################################################################################
Build C# Solution public/src/Public.sln in Release configuration. Then
1, Install npm on your machine (https://www.npmjs.com/).  (package manager for Node.js)
2, In this folder, open console and run
	npm install 
3, Then run: 
	npm run buildAndStart
4, If this doesnt work, try: 
	1,	Delete folder node_modules
	2,  run: npm install --scripts-prepend-node-path
	3,	run: npm run buildAndStart --scripts-prepend-node-path

Optional:
1, If you want to start without building, run: 
	npm run start
2, If you want to build without starting, run:
	npm run build

############################################################################################
Build Installer with InnoSetup and install as Standalone, windows, 64bit: 
############################################################################################
1, Install npm on your machine (https://www.npmjs.com/).  (package manager for Node.js)
2, Install InnoSetup on your machine (http://www.jrsoftware.org/isinfo.php)
3, In this folder, open console and run
	npm install
	npm run build
	npm run package
4, Open Subfolder InnoSetup and compile BoSSSpad-setup.iss using Innosetup. 
5, Open Subfolder InnoSetup\installer and run BoSSSpadStandalone-setup.exe .