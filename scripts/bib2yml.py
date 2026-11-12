#!/usr/bin/env python3
"""Build _data/publications.yml from _data/publications.bib.

Usage (from the repo root):   python3 scripts/bib2yml.py
Needs Python 3 and PyYAML (pip install pyyaml). Commit both files afterwards.

Fields understood per entry:
  author, title, year          required
  venue                        full venue text as shown on the page (HTML like <i>..</i> allowed); wins over the next row
  journal|booktitle|howpublished, volume, number, pages
                               used to build the venue when `venue` is absent
  doi, eprint (arXiv id), pdf, code     become DOI / arXiv / PDF / Code buttons
  links                        extra buttons:  Label|url ; Label|url
"""
import re, sys, pathlib, yaml

root = pathlib.Path(__file__).resolve().parent.parent
SRC, OUT = root / '_data' / 'publications.bib', root / '_data' / 'publications.yml'

def parse_bib(text):
    entries, i, n = [], 0, len(text)
    while True:
        i = text.find('@', i)
        if i < 0: return entries
        m = re.compile(r'@(\w+)\s*[{(]\s*([^,\s]+)\s*,').match(text, i)
        if not m: i += 1; continue
        kind, key = m.group(1).lower(), m.group(2)
        i, fields = m.end(), {}
        if kind in ('comment', 'string', 'preamble'): continue
        while i < n:
            while i < n and text[i] in ' \t\r\n,': i += 1
            if i >= n or text[i] in '})': i += 1; break
            m = re.compile(r'([\w-]+)\s*=\s*').match(text, i)
            if not m: i += 1; continue
            name, i = m.group(1).lower(), m.end()
            if text[i] == '{':
                depth, j = 1, i + 1
                while depth:
                    depth += {'{': 1, '}': -1}.get(text[j], 0); j += 1
                val, i = text[i + 1:j - 1], j
            elif text[i] == '"':
                j = i + 1
                while text[j] != '"' or text[j - 1] == '\\': j += 1
                val, i = text[i + 1:j], j + 1
            else:
                m = re.compile(r'[^,}\s]+').match(text, i); val, i = m.group(0), m.end()
            fields[name] = val
        entries.append((kind, key, fields))

def tidy(s):
    s = re.sub(r'\s+', ' ', s).strip()
    s = s.replace('\\&', '&amp;').replace('---', '—').replace('--', '–')
    return re.sub(r'(?<!\\)[{}]', '', s)

def name(a):                       # "Last, F. M." or "F. M. Last"  ->  "F. M. Last"
    a = tidy(a)
    if ',' in a:
        last, first = [x.strip() for x in a.split(',', 1)]
        return f'{first} {last}'.strip()
    return a

def authors(s):
    names = [name(a) for a in re.split(r'\s+and\s+', s.strip())]
    return names[0] if len(names) == 1 else ' and '.join(names) if len(names) == 2 else ', '.join(names[:-1]) + ', and ' + names[-1]

def venue(f):
    if 'venue' in f: return re.sub(r'\s+', ' ', f['venue']).strip()
    j = f.get('journal') or f.get('booktitle') or f.get('howpublished') or f.get('note', '')
    s = f'<i>{tidy(j)}</i>' if j else ''
    if f.get('volume'): s += f' {f["volume"]}' + (f'({f["number"]})' if f.get('number') else '')
    if f.get('pages'): s += f', {tidy(f["pages"])}'
    return (s + f' ({f["year"]})').strip()

def links(f):
    out = []
    for lab, key, fmt in [('Journal', 'url', '{}'), ('DOI', 'doi', 'https://doi.org/{}'), ('arXiv', 'eprint', 'https://arxiv.org/abs/{}'),
                          ('PDF', 'pdf', '{}'), ('Code', 'code', '{}')]:
        if f.get(key): out.append({'label': lab, 'url': fmt.format(f[key].strip())})
    for part in f.get('links', '').split(';'):
        if '|' in part:
            lab, url = part.split('|', 1); out.append({'label': lab.strip(), 'url': url.strip()})
    return out

def main():
    entries = parse_bib(SRC.read_text(encoding='utf-8'))
    pubs = []
    for kind, key, f in entries:
        for req in ('author', 'title', 'year'):
            if req not in f: sys.exit(f'Entry {key}: missing "{req}"')
        p = {'year': int(f['year']), 'authors': authors(f['author']), 'title': tidy(f['title']), 'venue': venue(f)}
        if links(f): p['links'] = links(f)
        pubs.append(p)
    pubs.sort(key=lambda p: -p['year'])            # stable: file order kept within a year
    head = '# GENERATED from _data/publications.bib by scripts/bib2yml.py. Edit the .bib file, not this one.\n'
    OUT.write_text(head + yaml.safe_dump(pubs, allow_unicode=True, sort_keys=False, width=200), encoding='utf-8')
    print(f'Wrote {len(pubs)} publications to {OUT.relative_to(root)}')

if __name__ == '__main__': main()
