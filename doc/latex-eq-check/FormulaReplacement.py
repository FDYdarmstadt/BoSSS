import os
import re
import argparse
import hashlib
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent

#--------------- Mode 1: Checking Formulas (scan .XMLDOC and save .tex report)------------------

def fmt_member_header(name: str) -> str:
    return r"\texttt{\detokenize{" + name + "}}"

def needs_verbatim_kind(kind: str) -> bool:
    return kind in ("FENCE", "FBRACKET", "FDOLLAR", "DOXYENV")

def extract_equations_all(text: str):
    """
    Patterns:
      FENCE, DOLLAR2, BRACKET, PAREN, FBRACKET, FDOLLAR,
      PARENBT, BT_DLR, DLR_BT, BT_DLR_MIX, DLR_BT_SAME,
      DOLLAR, LATEXENV, DOXYENV
    """
    envs = r"(equation\*?|align\*?|alignat\*?|aligned\*?|gather\*?|multline\*?|flalign\*?|eqnarray\*?)"

    # $...$ matcher to prevent pairing a stray '$' with a far-away '$'
    DOLLAR_MAX = 600
    dollar_body = rf"(?:(?![<>])[\s\S]){{0,{DOLLAR_MAX}}}?"  # anything except '<' or '>', up to limit
    dollar_rx = re.compile(
        rf"(?<!\\)\$(?![\s$])"
        rf"{dollar_body}"   
        rf"(?<!\s)(?<!\\)\$(?!\$)",
        re.DOTALL,
    )

    patterns = [
        ("FENCE", re.compile(r"```[ \t]*math[^\n]*\n[\s\S]*?\n```", re.IGNORECASE)),
        ("DOLLAR2", re.compile(r"(?<!\\)\$\$(?!\$)[\s\S]*?(?<!\\)\$\$(?!\$)")),
        ("BRACKET", re.compile(r"\\\[[\s\S]*?\\\]")),
        ("FBRACKET", re.compile(r"\\f\[[\s\S]*?\\f\]")),
        ("PAREN", re.compile(r"\\\([\s\S]*?\\\)")),
        ("FDOLLAR", re.compile(r"\\f\$[\s\S]*?\\f\$")),
        ("PARENBT", re.compile(r"\(\$`[\s\S]*?`\$\)")),
        ("DLR_BT", re.compile(r"\$`[\s\S]*?`\$")),
        ("DLR_BT_SAME", re.compile(r"\$`[\s\S]*?\$`")),
        ("BT_DLR", re.compile(r"`\$[\s\S]*?`\$")),
        ("BT_DLR_MIX", re.compile(r"`\$[\s\S]*?\$`")),
        ("DOLLAR", dollar_rx),
        ("LATEXENV", re.compile(rf"\\begin{{{envs}}}[\s\S]*?\\end{{\1}}")),
        ("DOXYENV", re.compile(rf"\\f\{{{envs}\}}\{{[\s\S]*?\}}")),
    ]

    # To ensure different dollar delimeters wins over $...$
    PRIORITY = {
        "DOXYENV": 100,
        "FDOLLAR": 100,
        "FBRACKET": 100,
        "FENCE": 95,
        "LATEXENV": 90,
        "DOLLAR2": 80,
        "BRACKET": 80,
        "PAREN": 70,
        "PARENBT": 70,
        "DLR_BT": 70,
        "DLR_BT_SAME": 70,
        "BT_DLR": 70,
        "BT_DLR_MIX": 70,
        "DOLLAR": 10,
    }

    spans = []
    for kind, rx in patterns:
        for m in rx.finditer(text):
            spans.append((m.start(), m.end(), text[m.start():m.end()], kind))

    if not spans:
        return []

    # Sort: higher priority, then longer, then earlier
    spans.sort(key=lambda t: (-PRIORITY.get(t[3], 0), -(t[1] - t[0]), t[0]))

    kept, occupied = [], []
    def overlaps(a, b):
        return not (a[1] <= b[0] or a[0] >= b[1])

    for s, e, frag, kind in spans:
        if any(overlaps((s, e), occ) for occ in occupied):
            continue
        kept.append((s, e, frag, kind))
        occupied.append((s, e))

    kept.sort(key=lambda t: t[0])  # back to source order
    return [(frag, kind) for _, __, frag, kind in kept]


def check_formulas_xml():
    # Relative Path
    root_dir = (SCRIPT_DIR / "../../src").resolve()
    report_path = SCRIPT_DIR / "T1.tex"

    report_entries = []
    seen_hashes = set()  # skip duplicates by file content

    for dirpath, _, filenames in os.walk(root_dir):
        for file in filenames:
            if not file.lower().endswith(".xml"):
                continue

            full_path = os.path.join(dirpath, file)

            # Dedup by exact bytes
            try:
                with open(full_path, "rb") as fb:
                    raw = fb.read()
            except Exception as ex:
                report_entries += [
                    r"\paragraph{Skipped file:} " + fmt_member_header(full_path),
                    r"\begin{verbatim}", f"{ex}", r"\end{verbatim}", ""
                ]
                continue

            file_hash = hashlib.sha1(raw).hexdigest()
            if file_hash in seen_hashes:
                continue
            seen_hashes.add(file_hash)

            content = raw.decode("utf-8", errors="replace")

            # <member ...>...</member>
            members = re.findall(r'<member name="([^"]+)">([\s\S]*?)<\/member>', content)

            for member_name, member_block in members:
                eqs = extract_equations_all(member_block)
                if not eqs:
                    continue

                report_entries.append(fmt_member_header(member_name))
                report_entries.append("Equations:")

                for eq_text, kind in eqs:
                    report_entries.append(eq_text)
                report_entries.append("")



    report_path.parent.mkdir(parents=True, exist_ok=True)
    with open(report_path, 'w', encoding='utf-8', newline='') as f:
        f.write("\\documentclass{article}\n"
                "\\usepackage{amsmath, amssymb}\n"
                "\\usepackage[T1]{fontenc}\n"
                "\\usepackage[utf8]{inputenc}\n"
                "\\begin{document}\n\n")
        f.write('\n'.join(report_entries))
        f.write("\n\\end{document}\n")

    print(f"Wrote {len(report_entries)} lines to: {report_path}")

#-------------- Mode 2: Markdown Equation Delimiter Conversion (UNCHANGED) ------------------

def fix_md_equations():
    md_root = (SCRIPT_DIR / "../xmldoc-to-md/Trial5").resolve()

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
                content = strip_spaces_in_inline_math(content)

                if content != original:
                    with open(full_path, 'w', encoding='utf-8') as f:
                        f.write(content)
                    modified_count += 1

    print(f"\n[✔] Markdown cleanup complete. {modified_count} files updated.")

# ------------------ MAIN MENU ------------------

def main():
    parser = argparse.ArgumentParser(description="Choose a mode of operation")
    parser.add_argument("choice", nargs="?", choices=["1", "2"], help="Select 1 or 2")
    args = parser.parse_args()

    if args.choice:
        choice = args.choice
    else:
        print("\nChoose a mode:\n")
        print("1 - Checking Formulas (scan .XMLDOC and save .tex report)")
        print("2 - Markdown Equation Delimiter Conversion required for GitLab documentation rendering (fixes delimiters in .md files)")
        choice = input("\nEnter 1 or 2: ").strip()

    if choice == '1':
        check_formulas_xml()
    elif choice == '2':
        fix_md_equations()
    else:
        print("\n[ERROR] Invalid input. Please enter 1 or 2.")

if __name__ == "__main__":
    main()
