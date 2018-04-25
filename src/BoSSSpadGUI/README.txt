Install and Run: 

Build C# Solution public/src/Public.sln in Release configuration. Then
1, Install npm on your machine (https://www.npmjs.com/). 
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