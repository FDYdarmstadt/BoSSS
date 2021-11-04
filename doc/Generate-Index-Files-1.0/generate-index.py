import os

indexTextStart = """<!DOCTYPE html>
<html>
<head><title>Index of {folderPath}</title></head>
<body>
    <h2>Index of {folderPath}</h2>
    <hr>
    <ul>
		<li>
			<a href='../'>../</a>
		</li>
"""
indexTextEnd = """
	</ul>
</body>
</html>
"""

def index_folder(folderPath):
	print("Indexing: " + folderPath +'/')
	#Getting the content of the folder
	files = os.listdir(folderPath)
	#If Root folder, correcting folder name
	root = folderPath
	if folderPath == '.':
		root = 'Root'
	indexText = indexTextStart.format(folderPath=root)
	for file in files:
		#Avoiding index.html files
		if file != 'index.html':
			indexText += "\t\t<li>\n\t\t\t<a href='" + file + "'>" + file + "</a>\n\t\t</li>\n"
		#Recursive call to continue indexing
		if os.path.isdir(folderPath+'/'+file):
			index_folder(folderPath + '/' + file)
	indexText += indexTextEnd
	#Create or override previous index.html
	index = open(folderPath+'/index.html', "w")
	#Save indexed content to file
	index.write(indexText)

#Indexing root directory (Script position)
index_folder('.')
