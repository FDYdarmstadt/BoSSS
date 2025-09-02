File: FormulaReplacement.py

Destination: d:\Users\monzer\Documents\BoSSS-master\public\doc\latex-eq-check


Functions: It have two options

(A) Formula Verification via XML Documentation: 
	It is needed to ensures all formulas are verified through a two-stage pipeline:
		1. Extraction – Collect equations from the XMLDOC.
		2. Compilation – Validate them by compiling with pdflatex.

	How it works? Once you run the tool, choose option 1 [Checking Formulas (scan .XMLDOC and save .tex report)]
	1. It Scans all the .XMLDOC files from the source directory (../../src).
	2. It extracts all the equations with there delimeters along with their associated members (fields/properties) for easy reference.
	3. It Generate a .tex report listing all extracted equations (T1.tex).
	4. Compile the report using pdflatex.
		4.1 If compilation succeeds → formulas are valid.
		4.2 If compilation fails → the LaTeX source code must be manually corrected in the source-code.


(B) Markdown Equation Delimiter Conversion required for GitLab documentation rendering (Overwrites the original .md files, does not save them into different directory)
	It is needed to convert the equation delimeters of the original .md files (extracted from the converted Markdown based on classes from the XMLDOCTOMD.py:
	[OUTPUT_DIR = r'./xmldoc-to-md/Trial6]), to the laTeX equation delimiters required for GitLab documentation rendering.

	How it works? Once you run the tool, choose option 2 [Markdown Equation Delimiter Conversion required for GitLab documentation rendering]
	1. It Process all Markdown (.md) files to fix equation delimiters.
	2. It replace the existing equation delimiters with the proper delimiters required for GitLab documentation rendering.
	3. It overwrite the original .md files with updated versions that follow the correct GitLab documentation format.
	4. It guarantees consistent and correct rendering of inline and block equations across GitLab documentation.


   

