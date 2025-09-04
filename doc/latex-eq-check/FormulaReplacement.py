import os
import re
from pathlib import Path
import argparse


SCRIPT_DIR = Path(__file__).resolve().parent
#-------------- Mode 1: Checking Formulas (scan .XMLDOC and save .tex report)------------------

# ----- helpers -----
def fmt_member_header(name: str) -> str:
    # print literally inside \texttt (underscores/backslashes safe)
    return r"\texttt{\detokenize{" + name + "}}"

def needs_verbatim_kind(kind: str) -> bool:
    # non-LaTeX delimiters that LaTeX can't typeset
    return kind in ("FENCE", "FBRACKET", "FDOLLAR", "DOXYENV")

# ----- robust extractor -----
def extract_equations_all(text: str):
    """
    Pairs: FENCE, DOLLAR2, BRACKET, PAREN, FBRACKET, FDOLLAR,
           PARENBT, BT_DLR, DLR_BT, BTBT, BT_DLR_MIX, DLR_BT_SAME,
           DOLLAR, LATEXENV, DOXYENV
    """
    eqs = []
    n = len(text)
    i = 0

    fence_re = re.compile(r'```[ \t]*math\b', re.IGNORECASE)

    # token pairs (start_token, end_token)
    PAIRS = [
        ("FBRACKET", "\\f[", "\\f]"),
        ("FDOLLAR",  "\\f$", "\\f$"),
        ("BRACKET",  "\\[",  "\\]"),
        ("PAREN",    "\\(",  "\\)"),
        ("PARENBT",  "(\\$`", "`\\$)"),   # ($`...`$)
        ("DLR_BT",   "$`",   "`$"),       # $`...`$
        ("DLR_BT_SAME", "$`", "$`"),      # $`...$`
        ("BT_DLR",   "`$",   "`$"),       # `$...`$
        ("BT_DLR_MIX", "`$", "$`"),       # `$...$`
    ]

    # LaTeX environments
    envs = r"(equation\*?|align\*?|alignat\*?|aligned\*?|gather\*?|multline\*?|flalign\*?|eqnarray\*?)"
    latex_env_re = re.compile(rf"\\begin{{{envs}}}[\s\S]*?\\end{{\1}}")
    doxy_env_re  = re.compile(rf"\\f\{{{envs}\}}\{{[\s\S]*?\}}")  # \f{env}{ ... }

    def find_single_dollar(s, start):
        idx = s.find('$', start)
        while idx != -1:
            # skip $$ and escaped \$
            if (idx + 1 < n and s[idx + 1] == '$') or (idx > 0 and s[idx - 1] == '\\'):
                idx = s.find('$', idx + 1)
                continue
            return idx
        return -1

    while i < n:
        candidates = []

        # ``` math
        m = fence_re.search(text, i)
        if m:
            candidates.append(("FENCE", m.start(), m.end()))

        # fixed tokens
        for kind, tok, _ in PAIRS:
            j = text.find(tok, i)
            if j != -1:
                candidates.append((kind, j, j + len(tok)))

        # $$ (treat distinctly so we don't mix with single $)
        j = text.find("$$", i)
        if j != -1:
            candidates.append(("DOLLAR2", j, j + 2))

        # single $
        j = find_single_dollar(text, i)
        if j != -1:
            candidates.append(("DOLLAR", j, j + 1))

        # LaTeX environments and Doxygen envs
        m_env  = latex_env_re.search(text, i)
        if m_env:
            candidates.append(("LATEXENV", m_env.start(), m_env.end()))
        m_doxy = doxy_env_re.search(text, i)
        if m_doxy:
            candidates.append(("DOXYENV", m_doxy.start(), m_doxy.end()))

        if not candidates:
            break

        # earliest start wins
        kind, s, start_end = min(candidates, key=lambda t: t[1])

        # find matching close for each kind
        if kind == "FENCE":
            close = text.find("```", start_end)
            e = n if close == -1 else close + 3

        elif kind == "DOLLAR2":
            close = text.find("$$", start_end)
            e = n if close == -1 else close + 2

        elif kind == "FBRACKET":
            close = text.find("\\f]", start_end)
            e = n if close == -1 else close + 3

        elif kind == "FDOLLAR":
            close = text.find("\\f$", start_end)
            e = n if close == -1 else close + 3

        elif kind == "BRACKET":
            close = text.find("\\]", start_end)
            e = n if close == -1 else close + 2

        elif kind == "PAREN":
            close = text.find("\\)", start_end)
            e = n if close == -1 else close + 2

        elif kind == "PARENBT":
            close = text.find("`\\$)", start_end)
            e = n if close == -1 else close + 4

        elif kind == "DLR_BT":
            close = text.find("`$", start_end)
            e = n if close == -1 else close + 2

        elif kind == "DLR_BT_SAME":
            close = text.find("$`", start_end)
            e = n if close == -1 else close + 2

        elif kind == "BT_DLR":
            close = text.find("`$", start_end)
            e = n if close == -1 else close + 2

        elif kind == "BT_DLR_MIX":
            close = text.find("$`", start_end)
            e = n if close == -1 else close + 2

        elif kind == "DOLLAR":
            k = start_end
            while True:
                k = text.find('$', k)
                if k == -1:
                    e = n
                    break
                if (k + 1 < n and text[k + 1] == '$') or (k > 0 and text[k - 1] == '\\'):
                    k += 1
                    continue
                e = k + 1
                break

        elif kind in ("LATEXENV", "DOXYENV"):
            e = start_end  # already points to end via regex

        frag = text[s:e]
        eqs.append((s, e, frag, kind))
        i = e

    # Resolve overlaps: keep longest first, then earliest start
    if not eqs:
        return []
    eqs.sort(key=lambda t: (-(t[1]-t[0]), t[0]))
    kept, occupied = [], []
    def overlaps(a, b):
        return not (a[1] <= b[0] or a[0] >= b[1])
    for s, e, txt, kind in eqs:
        if any(overlaps((s, e), (ks, ke)) for ks, ke in occupied):
            continue
        occupied.append((s, e))
        kept.append((s, e, txt, kind))
    kept.sort(key=lambda t: t[0])  # back to source order
    # De-dup by exact text+kind
    out, seen = [], set()
    for _, __, txt, kind in kept:
        key = (txt.strip(), kind)
        if txt.strip() and key not in seen:
            out.append((txt.strip(), kind))
            seen.add(key)
    return out

# ----- main -----
def check_formulas_xml():
    
    # Absolute Path
    #root_dir = r'd:\Users\monzer\Documents\BoSSS-master\public\src'
    #report_path = r'd:\Users\monzer\Documents\BoSSS-master\public\doc\latex-eq-check\Listingequationstrials\T2\nt.tex'

    # Relative Path
    SCRIPT_DIR = Path(__file__).resolve().parent
    root_dir = (SCRIPT_DIR / "../../src").resolve()
    report_path = SCRIPT_DIR / "T1.tex"

    report_entries = []
    seen_files = set()

    for dirpath, _, filenames in os.walk(root_dir):
        for file in filenames:
            if not file.lower().endswith(".xml"):
                continue
            if file.lower() in seen_files:
                continue
            seen_files.add(file.lower())

            full_path = os.path.join(dirpath, file)
            with open(full_path, 'r', encoding='utf-8') as f:
                content = f.read()

            # <member ...>
            members = re.findall(r'<member name="([^"]+)">([\s\S]*?)<\/member>', content)

            for member_name, member_block in members:
                eqs = extract_equations_all(member_block)
                if not eqs:
                    continue  # only list members that have equations

                report_entries.append(fmt_member_header(member_name))
                report_entries.append("Equations:")
                report_entries.append("")

                for eq_text, kind in eqs:
                    if needs_verbatim_kind(kind):
                        report_entries.append("\\begin{verbatim}")
                        report_entries.append(eq_text)
                        report_entries.append("\\end{verbatim}")
                    else:
                        report_entries.append(eq_text)
                    report_entries.append("")

                report_entries.append("")

    os.makedirs(os.path.dirname(report_path), exist_ok=True)
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write("\\documentclass{article}\n"
                "\\usepackage{amsmath, amssymb}\n"
                "\\usepackage[T1]{fontenc}\n"
                "\\usepackage[utf8]{inputenc}\n"
                "\\begin{document}\n\n")
        f.write('\n'.join(report_entries))
        f.write("\n\\end{document}\n")

    print(f"Wrote {len(report_entries)} lines to: {report_path}")

#if __name__ == '__main__':
   # check_formulas_xml()



#-------------- Mode 2: Markdown Equation Delimiter Conversion required for GitLab documentation rendering (fix delimiters in .md files)------------------

def fix_md_equations():
    # Absolute Path
    # md_root = r'd:\\Users\\monzer\\Documents\\BoSSS-master\\public\\doc\\xmldoc-to-md\\Trial2'
    # Relative Path
    # md_root = r'../../src'
    # SCRIPT_DIR = Path(__file__).resolve().parent
    md_root = (SCRIPT_DIR / "../xmldoc-to-md/Trial6").resolve()

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
    parser = argparse.ArgumentParser(description="Choose a mode of operation")
    parser.add_argument("choice", nargs="?", choices=["1", "2"], help="Select 1 or 2")
    args = parser.parse_args()

    if args.choice:
        choice = args.choice
    else:
        print("\nChoose a mode:\n")
        print("1 - Checking Formulas (scan .XMLDOC and save .tex report)")
        print("2 - Markdown Equation Delimiter Conversion required for GitLab documentation rendering (fix delimiters in .md files)")
        choice = input("\nEnter 1 or 2: ").strip()

    if choice == '1':
        check_formulas_xml()
    elif choice == '2':
        fix_md_equations()
    else:
        print("\n[ERROR] Invalid input. Please enter 1 or 2.")

if __name__ == "__main__":
    main()