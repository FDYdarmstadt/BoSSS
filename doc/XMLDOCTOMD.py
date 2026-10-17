# XMLDOC to Markdown

# Select mode based on number:
#   Class: 1 (Converts XMLDOC into Class-based Markdown documentation [1 Markdown file per class] used for the gitlab documentation)
#   Namespace: 2 (Converts XMLDOC into Namespace-based Markdown documentation (1 Markdown file per namespace) used for gpt knowledge)

import os
import xml.etree.ElementTree as ET
import re
from collections import defaultdict
from tqdm import tqdm
from pathlib import Path  # For path handling
import argparse

# ---------------- Paths ----------------
# Relative Path:
XML_DIR = r'../src'
OUTPUT_DIR = r'./xmldoc-to-md/Trial5'

# Absolute Path:
# XML_DIR = r"d:\\Users\\monzer\\Documents\\BoSSS-master\\public\\src"
# OUTPUT_DIR = r"d:\\Users\\monzer\\Documents\\BoSSS-master\\public\\doc\\xmldoc-to-md\\Trial5"

# Resolve relative paths against this script
BASE_DIR = Path(__file__).resolve().parent
XML_DIR = str((BASE_DIR / XML_DIR).resolve())
OUTPUT_DIR = str((BASE_DIR / OUTPUT_DIR).resolve())

# Ensure output directory exists, if not, it will create it
os.makedirs(OUTPUT_DIR, exist_ok=True)
print(f"Using XML_DIR   = {XML_DIR}")
print(f"Using OUTPUT_DIR= {OUTPUT_DIR}")
if not os.path.exists(XML_DIR):
    raise FileNotFoundError(f"XML_DIR does not exist: {XML_DIR}")

# cross-page links
NS_LINKS_WITH_MD = True

def sanitize_filename(name):  # Cleans file names for saving to disk
    name = re.sub(r'[<>:"/\\|?*]', '_', name)
    name = re.sub(r'[\(\),@]', '_', name)
    return name[:150]

def _ns_page_href(ns_page):
    fname = sanitize_filename(ns_page) + ".md"
    return fname if NS_LINKS_WITH_MD else fname[:-3]

# --------------- CLASS MODE ---------------

class_data = defaultdict(list)
crossref_map = {}
seen_members = set()
known_types = set()  # T: Strict 1-file-per-class)

# ---- Make class-mode anchors match Markdown/GitLab auto heading ids ----
def slugify_heading(text: str) -> str:
    s = text.lower().strip()
    s = re.sub(r'[^\w\s-]', '', s)   # drop punctuation (keep letters/digits/_/space/-)
    s = s.replace('_', '-')          # underscores -> hyphens
    s = re.sub(r'\s+', '-', s)       # spaces -> hyphen
    s = re.sub(r'-{2,}', '-', s)     # collapse multiple hyphens
    return s.strip('-')

def make_anchor(prefix, member_name):  # Creates an anchor that matches the VISIBLE heading text
    visible = f"{prefix.capitalize()}: {member_name}"
    return slugify_heading(visible)

def cref_to_link(cref_name):
    if cref_name in crossref_map:
        filename, anchor = crossref_map[cref_name]
        page = filename.replace(".md", "")
        return f"[{cref_name}]({page}#{anchor})"
    return f"**{cref_name}**"

def extract_text_crossref(element):
    if element is None or not hasattr(element, 'iter'):
        return ""
    result = ""

    def recursive_extract(node):
        nonlocal result
        if node.text:
            result += node.text
        for child in node:
            if child.tag == "see" and "cref" in child.attrib:
                cref = child.attrib["cref"]
                cref_name = cref.split(":")[-1]
                result += cref_to_link(cref_name)
            elif child.tag == "paramref" and "name" in child.attrib:
                param_name = child.attrib["name"]
                result += f"'{param_name}'"
            recursive_extract(child)
        if node.tail:
            result += node.tail

    recursive_extract(element)
    return "\n".join(line.lstrip() for line in result.splitlines())

def extract_all_summaries(member):
    summaries = member.findall("summary")
    all_summaries = []
    for summary in summaries:
        text = extract_text_crossref(summary).strip()
        if text:
            all_summaries.append(text)
    return "\n\n".join(all_summaries)

# ---- helper used ONLY for class-mode crossref building ----
def _split_for_class_crossref(full_name: str):
    
    # 1) original
    m = re.match(r"([\w\.]+)\.(.+)", full_name)
    if m:
        return m.group(1), m.group(2)

    # 2) extended type pattern
    m = re.match(r"^([\w\.\+`]+)\.(.+)$", full_name)
    if m:
        return m.group(1), m.group(2)

    # 3) fallback: last dot before '('
    paren = full_name.find('(')
    cutoff = paren if paren != -1 else len(full_name)
    cls, sep, mem_no_params = full_name[:cutoff].rpartition('.')
    if sep:
        return cls, mem_no_params + full_name[cutoff:]
    return None, None

def collect_crossrefs(xml_files):  # class-mode crossref map (anchors now match heading text)
    for xml_path in xml_files:
        try:
            tree = ET.parse(xml_path)
            root = tree.getroot()
            for member in root.findall("members/member"):
                name = member.attrib.get("name", "")
                if not name:
                    continue
                kind = name[0]
                full_name = name[2:]

                if kind == "T":
                    filename = sanitize_filename(full_name) + ".md"
                    # Anchor must match the auto id for the class heading: "# {full_name}"
                    anchor = slugify_heading(full_name)
                    crossref_map[full_name] = (filename, anchor)

                elif kind in "FPM":
                    class_name, member_signature = _split_for_class_crossref(full_name)
                    if not class_name:
                        continue

                    if class_name not in crossref_map:
                        class_file = sanitize_filename(class_name) + ".md"
                        class_anchor = slugify_heading(class_name)  # matches "# {class_name}"
                        crossref_map[class_name] = (class_file, class_anchor)

                    file_name, _ = crossref_map[class_name]
                    prefix = {"F": "field", "P": "property", "M": "method"}.get(kind, "")
                    anchor = make_anchor(prefix, member_signature)  # matches visible member heading text
                    crossref_map[full_name] = (file_name, anchor)

        except Exception as e:
            print(f"❌ Error reading {xml_path}: {e}")

# ---- helper for grouping ONLY ----
def split_type_and_member(full_name: str):

    paren = full_name.find('(')
    cutoff = paren if paren != -1 else len(full_name)
    class_name, sep, member_no_params = full_name[:cutoff].rpartition('.')
    if not sep:
        return None
    member_signature = member_no_params + full_name[cutoff:]
    return class_name, member_signature

def index_class_types(xml_files):  # Pass 1: one file per T:
    global known_types
    for xml_path in xml_files:
        try:
            tree = ET.parse(xml_path)
            root = tree.getroot()
            for member in root.findall("members/member"):
                name = member.attrib.get("name", "")
                if not name or name[0] != "T":
                    continue
                type_full = name[2:]
                known_types.add(type_full)
                if type_full in class_data:
                    continue
                entry = [f"# {type_full}"]
                summaries_text = extract_all_summaries(member)
                if summaries_text:
                    entry.append(f"**Summary:** {summaries_text}")

                t_remarks = member.find("remarks")
                if t_remarks is not None:
                    remarks_text = extract_text_crossref(t_remarks).lstrip().rstrip()
                    if remarks_text:
                        entry.append("")  # new paragraph before each remark
                        entry.append(f"**Remark:**\n{remarks_text}")

                for param in member.findall("param"):
                    pname = param.attrib.get("name", "")
                    ptext = extract_text_crossref(param).strip()
                    first_line = ptext.lstrip().splitlines()[0] if ptext.strip() else ""
                    entry.append("")  # new paragraph before each parameter
                    if first_line.startswith("-") or first_line.startswith("1st index") or "\n- " in ptext:
                        entry.append(f"**Parameter:** `{pname}` -\n{ptext}")
                    else:
                        entry.append(f"**Parameter:** `{pname}` - {ptext}")

                t_returns = member.find("returns")
                if t_returns is not None:
                    returns_text = extract_text_crossref(t_returns).lstrip().rstrip()
                    entry.append("")  # new paragraph before each return
                    entry.append(f"**Returns:**\n{returns_text}")

                class_data[type_full] = entry
        except Exception as e:
            print(f"❌ Error indexing types in {xml_path}: {e}")

def attach_members_to_types(xml_files):  # Pass 2: attach F/P/M to existing T: pages
    for xml_path in xml_files:
        try:
            tree = ET.parse(xml_path)
            root = tree.getroot()
            for member in root.findall("members/member"):
                name = member.attrib.get("name", "")
                if not name:
                    continue
                kind = name[0]
                if kind not in "FPM":
                    continue
                full_name = name[2:]

                split = split_type_and_member(full_name)
                if not split:
                    continue
                class_name, member_signature = split

                # strictly attach only if that class was created from T:
                if class_name not in class_data:
                    continue

                if kind == "F":
                    heading = f"### Field: {member_signature}"
                elif kind == "P":
                    heading = f"### Property: {member_signature}"
                elif kind == "M":
                    heading = f"## Method: {member_signature}"
                else:
                    heading = f"### {member_signature}"

                section = ["", heading]
                summaries_text = extract_all_summaries(member)
                if summaries_text:
                    section.append(f"**Summary:** {summaries_text}")

                if kind == "M":
                    for param in member.findall("param"):
                        pname = param.attrib.get("name", "")
                        ptext = extract_text_crossref(param).strip()
                        first_line = ptext.lstrip().splitlines()[0] if ptext.strip() else ""
                        section.append("")  # new paragraph before each parameter
                        if first_line.startswith("-") or first_line.startswith("1st index") or "\n- " in ptext:
                            section.append(f"**Parameter:** `{pname}` -\n{ptext}")
                        else:
                            section.append(f"**Parameter:** `{pname}` - {ptext}")
                    returns = member.find("returns")
                    if returns is not None:
                        returns_text = extract_text_crossref(returns).lstrip().rstrip()
                        section.append("")  # new paragraph before each returns
                        section.append(f"**Returns:**\n{returns_text}")
                    remarks = member.find("remarks")
                    if remarks is not None:
                        remarks_text = extract_text_crossref(remarks).lstrip().rstrip()
                        if remarks_text:
                            section.append("")  # new paragraph before each remark
                            section.append(f"**Remark:**\n{remarks_text}")
                elif kind in ("F", "P"):
                    remarks = member.find("remarks")
                    if remarks is not None:
                        remarks_text = extract_text_crossref(remarks).lstrip().rstrip()
                        if remarks_text:
                            section.append("")  # new paragraph before each remark
                            section.append(f"**Remark:**\n{remarks_text}")

                class_data[class_name].extend(section)

        except Exception as e:
            print(f"❌ Error attaching members in {xml_path}: {e}")

def save_class_files():  # Saves class-based Markdown documentation to disk
    for class_name, entries in class_data.items():
        safe_filename = sanitize_filename(class_name) + ".md"
        class_filepath = os.path.join(OUTPUT_DIR, safe_filename)
        with open(class_filepath, "w", encoding="utf-8") as md_file:
            md_file.write("\n".join(entries) + "\n")
        print(f"✅ Saved: {class_filepath}")

# --------------- NAMESPACE MODE ---------------

ns_link_map = {}

def build_namespace_link_map(xml_files):
    """Build a map so namespace pages can cross-link to other namespace pages."""
    tmp = {}
    for xml_path in xml_files:
        try:
            tree = ET.parse(xml_path)
            root = tree.getroot()
            assembly = root.find("assembly/name")
            if assembly is None or not assembly.text:
                continue
            ns_page = assembly.text.strip()
            for member in root.findall("members/member"):
                name = member.attrib.get("name", "")
                if not name:
                    continue
                kind = name[0]
                full_name = name[2:]
                if kind == "T":
                    tmp[full_name] = (ns_page, full_name.lower())
                elif kind in "FPM":
                    tmp[full_name] = (ns_page, full_name.lower())
                    m = re.match(r"^([\w\.\+`]+)\.(.+)$", full_name)
                    if m:
                        cls = m.group(1)
                        tmp.setdefault(cls, (ns_page, cls.lower()))
        except Exception as e:
            print(f"❌ Error building namespace link map from {xml_path}: {e}")
    return tmp

def extract_text_namespace(element, current_namespace):  # Namespace extractor
    if element is None or not hasattr(element, 'iter'):
        return ""
    result = ""

    def recursive_extract(node):
        nonlocal result
        if node.text:
            result += node.text
        for child in node:
            if child.tag == "see" and "cref" in child.attrib:
                cref = child.attrib["cref"]
                cref_name = cref.split(":")[-1]
                target = ns_link_map.get(cref_name)
                if target:
                    page_name, anchor = target
                    if page_name == current_namespace:
                        result += f"[{cref_name}](#{anchor})"
                    else:
                        page_href = _ns_page_href(page_name)
                        result += f"[{cref_name}]({page_href}#{anchor})"
                else:
                    result += f"**{cref_name}**"
            elif child.tag == "paramref" and "name" in child.attrib:
                param_name = child.attrib["name"]
                result += f"'{param_name}'"
            recursive_extract(child)
        if node.tail:
            result += node.tail

    recursive_extract(element)
    return "\n".join(line.lstrip() for line in result.splitlines())

namespace_data = defaultdict(list)
ns_type_started = defaultdict(set)

def process_by_namespace(xml_path):
    try:
        tree = ET.parse(xml_path)
        root = tree.getroot()
        assembly = root.find("assembly/name")
        if assembly is None or not assembly.text:
            return
        namespace = assembly.text.strip()

        members = root.findall("members/member")
        if not members:
            return

        def ensure_type_section(type_full):
            if type_full not in ns_type_started[namespace]:
                namespace_data[namespace].append(
                    f"## Class: {type_full} <a id=\"{type_full.lower()}\"></a>\n"
                )
                ns_type_started[namespace].add(type_full)

        for member in members:
            name = member.attrib.get("name", "")
            if not name or name in seen_members:
                continue
            seen_members.add(name)

            kind = name[0]
            full_name = name[2:]

            if kind == "T":
                ensure_type_section(full_name)

                summaries = member.findall("summary")
                if summaries:
                    all_summaries = []
                    for s in summaries:
                        txt = extract_text_namespace(s, namespace).strip()
                        if txt:
                            all_summaries.append(txt)
                    if all_summaries:
                        namespace_data[namespace].append(f"**Summary:** " + "\n\n".join(all_summaries) + "\n")

                t_remarks = member.find("remarks")
                if t_remarks is not None:
                    rtxt = extract_text_namespace(t_remarks, namespace).strip()
                    if rtxt:
                        namespace_data[namespace].append(f"**Remark:**\n{rtxt}\n")

            elif kind in ("F", "P", "M"):
                match = re.match(r"^([\w\.\+`]+)\.(.+)$", full_name)
                if not match:
                    cls_name, member_sig = full_name, full_name
                else:
                    cls_name, member_sig = match.groups()

                ensure_type_section(cls_name)

                if kind == "F":
                    heading = f"### Field: {full_name} <a id=\"{full_name.lower()}\"></a>"
                elif kind == "P":
                    heading = f"### Property: {full_name} <a id=\"{full_name.lower()}\"></a>"
                elif kind == "M":
                    heading = f"## Method: {full_name} <a id=\"{full_name.lower()}\"></a>"
                else:
                    heading = f"### {full_name} <a id=\"{full_name.lower()}\"></a>"

                section = ["", heading]

                summaries = member.findall("summary")
                if summaries:
                    all_summaries = []
                    for s in summaries:
                        txt = extract_text_namespace(s, namespace).strip()
                        if txt:
                            all_summaries.append(txt)
                    if all_summaries:
                        section.append(f"**Summary:** " + "\n\n".join(all_summaries))

                if kind == "M":
                    for p in member.findall("param"):
                        pname = p.attrib.get("name", "")
                        ptext = extract_text_namespace(p, namespace).strip()
                        first_line = ptext.lstrip().splitlines()[0] if ptext else ""
                        if first_line.startswith("-") or first_line.startswith("1st index") or "\n- " in ptext:
                            section.append(f"**Parameter:** `{pname}` -\n{ptext}")
                        else:
                            section.append(f"**Parameter:** `{pname}` - {ptext}")
                    returns = member.find("returns")
                    if returns is not None:
                        rtxt = extract_text_namespace(returns, namespace).lstrip().rstrip()
                        section.append(f"**Returns:**\n{rtxt}")
                    remarks = member.find("remarks")
                    if remarks is not None:
                        rtxt = extract_text_namespace(remarks, namespace).lstrip().rstrip()
                        if rtxt:
                            section.append(f"**Remark:**\n{rtxt}")
                else:
                    remarks = member.find("remarks")
                    if remarks is not None:
                        rtxt = extract_text_namespace(remarks, namespace).lstrip().rstrip()
                        if rtxt:
                            section.append(f"**Remark:**\n{rtxt}")

                namespace_data[namespace].append("\n".join(section) + "\n")

    except Exception as e:
        print(f"❌ Error processing {xml_path}: {e}")

def save_namespace_files():  # Saves namespace-based Markdown documentation to disk
    for namespace, entries in namespace_data.items():
        namespace_file = os.path.join(OUTPUT_DIR, f"{sanitize_filename(namespace)}.md")
        with open(namespace_file, "w", encoding="utf-8") as md_file:
            md_file.write(f"# Namespace: {namespace}\n\n")
            md_file.write("\n".join(entries))
        print(f"✅ Saved: {namespace_file}")


# ---------------- Main ----------------
def main():  # Entry point: loads XMLs, selects mode, processes accordingly
    parser = argparse.ArgumentParser(description="XMLDOC to Markdown")
    parser.add_argument("choice", nargs="?", choices=["1", "2"], help="Select 1 (Class) or 2 (Namespace)")
    args = parser.parse_args()

    if args.choice:
        mode_input = args.choice
    else:
        print("Select mode:")
        print("1 - Class")
        print("2 - Namespace")
        mode_input = input("Enter 1 or 2: ").strip()

    mode = "class" if mode_input == "1" else "namespace"

    all_files = [
        os.path.join(root, file)
        for root, _, files in os.walk(XML_DIR)
        for file in files
        if file.lower().endswith(".xml")
    ]

    seen_basenames = set()
    xml_files = []

    for path in all_files:
        base_name = os.path.splitext(os.path.basename(path))[0].lower()
        if base_name not in seen_basenames:
            seen_basenames.add(base_name)
            xml_files.append(path)

    print(f"🔍 Found {len(xml_files)} XML files")

    if mode == "class":
        # Strict 1-file-per-class:
        class_data.clear()
        crossref_map.clear()
        seen_members.clear()
        known_types.clear()

        # 1) Build crossrefs
        collect_crossrefs(xml_files)

        # 2) Create pages only from T:
        index_class_types(xml_files)

        # 3) Attach F/P/M to those pages (progress per file)
        for xml_file in tqdm(xml_files, desc="Processing by class (attach members)"):
            attach_members_to_types([xml_file])

        save_class_files()

    else:
        # Namespace: build link map once and process
        namespace_data.clear()
        ns_type_started.clear()
        seen_members.clear()

        global ns_link_map
        ns_link_map = build_namespace_link_map(xml_files)

        for xml_file in tqdm(xml_files, desc="Processing by namespace"):
            process_by_namespace(xml_file)
        save_namespace_files()
if __name__ == "__main__":
    main()