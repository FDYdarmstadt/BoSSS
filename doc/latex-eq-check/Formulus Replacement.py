import os
import re

# Root directory to scan XML files
root_dir = r'd:\Users\monzer\Documents\BoSSS-master\public\src'

# Output report file path
report_path = r'd:\Users\monzer\Documents\BoSSS Formulas\Formulus Replacement\equation_delimiterssss_report.txt'

# Inline math replacements
inline_patterns = [
    # 1. ($`...`$) → ('$...'$) not solved
    (re.compile(r"\(\$`(.*?)`\$\)"), r"('(\$\1'$\)"),
    # 2. \f$...$\f → '$...'$ solved
    (re.compile(r"\\f\$([^$]+?)\\f\$", re.DOTALL), r"'$\1'$"),
    # 3. $'...'$ → '$...'$ solved
    (re.compile(r"\$`(.*?)`\$"), r"'$\1'$"),
    # 4. $`...`$ → '$...'$
    (re.compile(r"\$'([^']+?)'\$"), r"'$\1'$"),
    # 4. $`...$` → '$...'$ not solved
    (re.compile(r"\$'([^']+?)\$'"), r"'$\1'$"),
    # 5. $...$ → '$...'$ solved
    (re.compile(r"\$([^']+?)\$"), r"'$\1'$"),
    ]

# Display math replacements
display_patterns = [
        # \[...\] → ```math ... ```
        (re.compile(r"\\\[(.*?)\\\]", re.DOTALL), r"```\nmath\n\1\n```"),

        # \f[...\f] → ```math ... ```
        (re.compile(r"\\f\[(.*?)\\f\]", re.DOTALL), r"```\nmath\n\1\n```"),

        # $$...$$ → ```math ... ```
        (re.compile(r"\$\$\s*(.*?)\s*\$\$", re.DOTALL), r"```\nmath\n\1\n```"),
    ]

# Collect report entries
report_entries = []

# Process all XML files
for dirpath, _, filenames in os.walk(root_dir):
    for file in filenames:
        if file.endswith(".xml"):
            full_path = os.path.join(dirpath, file)
            with open(full_path, 'r', encoding='utf-8') as f:
                lines = f.readlines()

            for idx, line in enumerate(lines):
                original = line
                modified = line

                # Apply inline replacements
                for pattern, repl in inline_patterns:
                    modified = re.sub(pattern, repl, modified)

                # Apply display replacements
                for pattern, repl in display_patterns:
                    modified = re.sub(pattern, repl, modified)

                # Report any changes
                if original != modified:
                    report_entries.append(
                        f"[{full_path}:{idx+1}]\nORIGINAL:\n{original.rstrip()}\nMODIFIED:\n{modified.rstrip()}\n"
                    )

# Ensure output directory exists
os.makedirs(os.path.dirname(report_path), exist_ok=True)

# Write the report
with open(report_path, 'w', encoding='utf-8') as f:
    f.write('\n'.join(report_entries))

print(f"Scan complete. {len(report_entries)} modified entries written to:\n{report_path}")
