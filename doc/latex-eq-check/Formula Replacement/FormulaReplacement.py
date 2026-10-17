import os
import re

# ------------------ OPTION 1: XML CHECK ------------------
def escape_latex(s):
    """
    Escapes LaTeX special characters in a string, such as underscores or percent signs.
    Used only for escaping member names.
    """
    return s.replace('_', r'\_').replace('%', r'\%').replace('&', r'\&').replace('#', r'\#') \
            .replace('{', r'\{').replace('}', r'\}').replace('~', r'\textasciitilde{}') \
            .replace('^', r'\textasciicircum{}').replace('\\', r'\textbackslash{}')

def check_formulas_xml():
    # === Paths ===
    root_dir = r'd:\Users\monzer\Documents\BoSSS-master\public\src'
    report_path = r'd:\Users\monzer\Documents\BoSSS-master\public\doc\latex-eq-check\Formula Replacement\T1\Trial1.tex'
    #root_dir = r'../../src'
    #report_path = r'./srccodeeqreplacement-2.tex'

    # === Inline formula patterns ===
    inline_patterns = [
        (re.compile(r"\(\$`(.*?)`\$\)"), r"($\1$"),                     
        (re.compile(r"\\f\$([^']+?)\\f\$", re.DOTALL), r"$\1$"),       
        (re.compile(r"\$`([^']+?)`\$"), r"$\1$"),                     
        (re.compile(r"\$`([^']+?)\$`"), r"$\1$"),                     
        (re.compile(r"`\$([^']+?)`\$"), r"$\1$"),                      
        (re.compile(r"`\$([^']+?)\$`"), r"$\1$"),                    
        (re.compile(r"\$([^']+?)\$"), r"$\1$"),                        
    ]

    # === Display formula patterns ===
    display_block_start = re.compile(r"```math\s*")                  # Start of ```math
    display_block_end = re.compile(r"```")                           # End of ``` block
    display_bracket_pattern = re.compile(r"\\f\[(.*?)\\f\]", re.DOTALL)  # \f[...\f] → \[...\]

    report_entries = []
    seen_files = set()  # To avoid processing duplicate filenames

    # === Traverse XML files ===
    for dirpath, _, filenames in os.walk(root_dir):
        for file in filenames:
            if not file.lower().endswith(".xml"):
                continue

            base_name = file.lower()
            if base_name in seen_files:
                continue  # Skip duplicated files
            seen_files.add(base_name)

            full_path = os.path.join(dirpath, file)
            with open(full_path, 'r', encoding='utf-8') as f:
                content = f.read()

            # === Extract <member name="...">...</member> blocks ===
            members = re.findall(r'<member name="([^"]+)">([\s\S]*?)<\/member>', content)

            for member_name, member_block in members:
                lines = member_block.splitlines()
                modified_equations = []
                in_math_block = False
                block_lines = []

                # === Process each line in member block ===
                for line in lines:
                    stripped = line.strip()

                    if not in_math_block:
                        # Check for ```math block start
                        if display_block_start.match(stripped):
                            in_math_block = True
                            block_lines = [line]
                        else:
                            # Apply inline transformations
                            modified = line
                            for pattern, repl in inline_patterns:
                                modified = re.sub(pattern, repl, modified)

                            # Replace \f[...\f] blocks with LaTeX display math
                            if display_bracket_pattern.search(modified):
                                modified = display_bracket_pattern.sub(r"\\[\1\\]", modified)

                            # Extract modified math parts
                            if modified != line:
                                equations = re.findall(r"(\$\s*.*?\s*\$|\\\[.*?\\\])", modified, re.DOTALL)
                                modified_equations.extend(eq.strip() for eq in equations)
                    else:
                        # Inside a ```math block
                        block_lines.append(line)
                        if display_block_end.match(stripped):
                            math_content = '\n'.join(block_lines[1:-1])
                            display_eq = f"\\[\n{math_content.strip()}\n\\]"
                            modified_equations.append(display_eq)
                            in_math_block = False
                            block_lines = []

                if modified_equations:
                    # Escape member name for LaTeX header
                    safe_name = escape_latex(member_name)
                    report_entries.append(f"\\texttt{{{safe_name}}}")
                    report_entries.append("Modified:")
                    for eq in modified_equations:
                        if eq.strip():  # Avoid empty math expressions
                            report_entries.append(eq.strip())
                    report_entries.append("")  # Blank line between members

    # === Write output LaTeX file ===
    os.makedirs(os.path.dirname(report_path), exist_ok=True)
    with open(report_path, 'w', encoding='utf-8') as f:
        # Start LaTeX document
        f.write("\\documentclass{article}\n"
                "\\usepackage{amsmath, amssymb}\n"
                "\\usepackage[T1]{fontenc}\n"
                "\\usepackage[utf8]{inputenc}\n"
                "\\begin{document}\n\n")

        f.write('\n'.join(report_entries))

        # End LaTeX document
        f.write("\n\n\\end{document}")

    print(f"Scan complete. {len(report_entries)} modified entries written to:\n{report_path}")

# ------------------ OPTION 2: MODIFY MARKDOWN FILES ------------------

def fix_md_equations():
    md_root = r'd:\\Users\\monzer\\Documents\\BoSSS-master\\public\\doc\\xmldoc-to-md\\Trial2'

    inline_patterns = [
        (re.compile(r"\(\$`(.*?)`\$\)"), r"($\1$)"),
        (re.compile(r"\\f\$([^']+?)\\f\$", re.DOTALL), r"$\1$"), 
        (re.compile(r"\$`([^']+?)`\$"), r"$\1$"),
        (re.compile(r"\$`([^']+?)\$`"), r"$\1$"), 
        (re.compile(r"`\$([^']+?)`\$"), r"$\1$"), 
        (re.compile(r"`\$([^']+?)\$`"), r"$\1$"),
        (re.compile(r"\$([^']+?)\$"), r"$\1$"),
    ]

    display_patterns = [
        (re.compile(r"\\\[(.*?)\\\]", re.DOTALL), r"\n```math\n\1\n```"),
        (re.compile(r"\\f\[(.*?)\\f\]", re.DOTALL), r"\n```math\n\1\n```"),
        (re.compile(r"\$\$\s*(.*?)\s*\$\$", re.DOTALL), r"\n```math\n\1\n```"),
    ]

    def strip_spaces_in_inline_math(text):
        """
        Remove leading/trailing spaces inside $...$ math expressions.
        """
        return re.sub(r"\$(\s*[^$]*?\s*)\$", lambda m: f"${m.group(1).strip()}$", text)

    modified_count = 0

    for dirpath, _, filenames in os.walk(md_root):
        for file in filenames:
            if file.endswith(".md"):
                full_path = os.path.join(dirpath, file)
                with open(full_path, 'r', encoding='utf-8') as f:
                    content = f.read()

                original = content
                for pattern, repl in inline_patterns:
                    content = re.sub(pattern, repl, content)
                for pattern, repl in display_patterns:
                    content = re.sub(pattern, repl, content)
                content = strip_spaces_in_inline_math(content) # Clean spaces inside math expressions

                if content != original:
                    with open(full_path, 'w', encoding='utf-8') as f:
                        f.write(content)
                    modified_count += 1

    print(f"\n[✔] Markdown cleanup complete. {modified_count} files updated.")

# ------------------ MAIN MENU ------------------

def main():
    print("\nChoose a mode:\n")
    print("1 - Checking Formulas (scan .XMLDOC and save .tex report)")
    print("2 - Markdown Equation Delimiter Conversion required for GitLab documentation rendering)")

    choice = input("\nEnter 1 or 2: ").strip()

    if choice == '1':
        check_formulas_xml()
    elif choice == '2':
        fix_md_equations()
    else:
        print("\n[ERROR] Invalid input. Please enter 1 or 2.")

if __name__ == "__main__":
    main()
