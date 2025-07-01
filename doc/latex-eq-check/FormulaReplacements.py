import os
import re

# ------------------ OPTION 1: XML CHECK ------------------
def check_formulas_xml():
    root_dir = r'..\..\src'
    report_path = r'.\Eq-replacement.tex'

    # === INLINE FORMULA PATTERNS ===
    inline_patterns = [
        (re.compile(r"\(\$`(.*?)`\$\)"), r"($\1$)"),        
        (re.compile(r"\\f\$([^']+?)\\f\$", re.DOTALL), r"$\1$"), 
        (re.compile(r"\$`([^']+?)`\$"), r"$\1$"),             
        (re.compile(r"\$`([^']+?)\$`"), r"$\1$"),             
        (re.compile(r"`\$([^']+?)`\$"), r"$\1$"),           
        (re.compile(r"`\$([^']+?)\$`"), r"$\1$"),              
        (re.compile(r"\$([^']+?)\$"), r"$\1$"),              
    ]

    # === DISPLAY FORMULA PATTERNS ===
    display_block_start = re.compile(r"```math\s*")
    display_block_end = re.compile(r"```")

    display_bracket_pattern = re.compile(r"\\f\[(.*?)\\f\]", re.DOTALL)

    report_entries = []      # To store transformation logs
    seen_files = set()       # To avoid reprocessing duplicated XML files

    for dirpath, _, filenames in os.walk(root_dir):
        for file in filenames:
            if not file.lower().endswith(".xml"):
                continue

            base_name = file.lower()
            if base_name in seen_files:
                continue  # To skip duplicate files
            seen_files.add(base_name)

            full_path = os.path.join(dirpath, file)
            with open(full_path, 'r', encoding='utf-8') as f:
                lines = f.readlines()

            in_math_block = False         
            block_start_idx = None      
            block_lines = []            

            # === PROCESS EACH LINE ===
            for idx, line in enumerate(lines):
                if not in_math_block:
                    if display_block_start.match(line.strip()):
                        # Start of ```math block
                        in_math_block = True
                        block_start_idx = idx
                        block_lines = [line]
                    else:
                        # Process inline patterns and \f[ ... \f]
                        original = line
                        modified = line

                        # Apply inline transformations
                        for pattern, repl in inline_patterns:
                            modified = re.sub(pattern, repl, modified)

                        # Apply display-style \f[ ... \f] transformation
                        if display_bracket_pattern.search(modified):
                            modified = display_bracket_pattern.sub(r"\\[\1\\]", modified)

                        # Report modified equations
                        if original != modified:
                            report_entries.append(
                                f"[{full_path}:{idx + 1}]\nORIGINAL:\n{original.rstrip()}\nMODIFIED:\n{modified.rstrip()}\n"
                            )
                else:
                    # Accumulate lines in a ```math block
                    block_lines.append(line)
                    if display_block_end.match(line.strip()):
                        # End of ```math block reached
                        original_block = ''.join(block_lines)
                        math_content = ''.join(block_lines[1:-1])  # Skip ```math and closing ```
                        modified_block = f"\\[{math_content.strip()}\\]\n"

                        report_entries.append(
                            f"[{full_path}:{block_start_idx + 1}-{idx + 1}]\nORIGINAL:\n{original_block.rstrip()}\nMODIFIED:\n{modified_block.rstrip()}\n"
                        )

                        # Reset block state
                        in_math_block = False
                        block_lines = []
                        
    os.makedirs(os.path.dirname(report_path), exist_ok=True)
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(report_entries))

    print(f"Scan complete. {len(report_entries)} modified entries written to:\n{report_path}")


# ------------------ OPTION 2: MODIFY MARKDOWN FILES ------------------

def fix_md_equations():
    md_root = r'd:\Users\monzer\Documents\.md\REPLACING EQ TRIALS\2'

    inline_patterns = [
        (re.compile(r"\(\$`(.*?)`\$\)"), r"('$\1'$)"),
        (re.compile(r"\\f\$([^']+?)\\f\$", re.DOTALL), r"'$\1'$"), 
        (re.compile(r"\$`([^']+?)`\$"), r"'$\1'$"),
        (re.compile(r"\$`([^']+?)\$`"), r"'$\1'$"), 
        (re.compile(r"`\$([^']+?)`\$"), r"'$\1'$"), 
        (re.compile(r"`\$([^']+?)\$`"), r"'$\1'$"),
        (re.compile(r"\$([^']+?)\$"), r"'$\1'$"),
    ]

    display_patterns = [
        (re.compile(r"\\\[(.*?)\\\]", re.DOTALL), r"\n```math\n\1\n```"),
        (re.compile(r"\\f\[(.*?)\\f\]", re.DOTALL), r"\n```math\n\1\n```"),
        (re.compile(r"\$\$\s*(.*?)\s*\$\$", re.DOTALL), r"\n```math\n\1\n```"),
    ]

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

                if content != original:
                    with open(full_path, 'w', encoding='utf-8') as f:
                        f.write(content)
                    modified_count += 1

    print(f"\n[✔] Markdown cleanup complete. {modified_count} files updated.")

# ------------------ MAIN MENU ------------------

def main():
    print("\nChoose a mode:\n")
    print("1 - Checking Formulas (scan XML and save .tex report)")
    print("2 - Convert XML to Markdown (fix delimiters in .md files)")

    choice = input("\nEnter 1 or 2: ").strip()

    if choice == '1':
        check_formulas_xml()
    elif choice == '2':
        fix_md_equations()
    else:
        print("\n[ERROR] Invalid input. Please enter 1 or 2.")

if __name__ == "__main__":
    main()
