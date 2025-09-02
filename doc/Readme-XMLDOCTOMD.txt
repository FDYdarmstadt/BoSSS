File: XMLDOCTOMD.py

Destination: d:\Users\monzer\Documents\BoSSS-master\public\doc

Functions: 
   Converts .XMLDOC to .md documents having 2 output options:
      1. Class: it Converts XMLDOC into class-based Markdown documentation (1 Markdown file per class) used for the gitlab documentation
      2. Namespace: it Converts XMLDOC  into namespace-based Markdown documentation (1 Markdown file per namespace) used for gpt knowledge
   If duplicate XML files are found, the tool processes only one copy and skips the rest to avoid duplicate Markdown output.
   Cross-references work across all generated .md files.


How it works?
   Once you run the tool, you need to choose a conversion mode: (Just enter a number 1 or 2)
      1. Class mode
      2. Namespace mode
   The tool scans the XML documentation files in:[XML_DIR = r'../src'] the directory where all the XMLDOC are placed
   and convert all this XMLDOC to Markdown according to the selected mode and writes the results to: [OUTPUT_DIR = r'./xmldoc-to-md/Trial6]

   If class mode (1) is chosen, conversion of XMLDOC → Markdown based on classes takes place, which is used for the gitlab documentation.
      Procedure to Publish to GitLab Wiki:
        1.1 Clone the project’s Wiki repository.
        1.2 Copy the generated Markdown pages into the wiki working tree.
        1.3 Commits and pushes the changes.
   
   If Namespace mode (2) is chosen, conversion of XMLDOC → Markdown based on namespace takes place, which is used for gpt knowledge.
      Upload the generated Markdown files to your Custom GPT’s knowledge base.
