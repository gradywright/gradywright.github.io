# Grady Wright's website

Source for https://gradywright.github.io. It is a [Jekyll](https://jekyllrb.com) site hosted on GitHub Pages: commit a change to the main branch and GitHub rebuilds the site in a minute or two. This file is not published.

## Where things live

| To change...                  | Edit this                                   |
| ----------------------------- | ------------------------------------------- |
| Home page text, hero, About   | `index.md`                                  |
| Publications                  | `_data/publications.bib`, then run a script |
| Talks, software on Research   | `research.html` (bottom of the file)        |
| Teaching history              | `_data/teaching.yml`                        |
| Students                      | `_data/students.yml` + photos               |
| Top menu links                | `_data/nav.yml`                             |
| Footer, site title, SEO text  | `_layouts/default.html`, `_config.yml`      |
| Colors, fonts, spacing        | `assets/css/style.css` (top of file: tokens) |
| Animations, menu behavior     | `assets/js/site.js`                         |
| CV, PDFs, papers              | `documents/`, `papers/`, `bibtex/`          |

Pages never need HTML for lists of papers, students or courses. Those are filled in from the data files above.

## Preview locally

```
gem install bundler jekyll webrick
jekyll serve
```

Open http://localhost:4000. It rebuilds when you save. (On a Mac you may need a newer Ruby than the system one: `brew install ruby` and add it to your PATH.)

## Adding a publication

Publications are written in BibTeX in `_data/publications.bib`.

1. Paste a new entry into the file. Order doesn't matter, the page sorts by year, newest first.
2. Run this from the repo folder (needs Python 3 and `pip install pyyaml`):
   ```
   python3 scripts/bib2yml.py
   ```
3. Commit both `_data/publications.bib` and `_data/publications.yml`. The second is generated, so don't edit it by hand.

```bibtex
@article{doe2027,
  author  = {Doe, Jane and Wright, Grady B.},
  title   = {Title of the paper},
  journal = {Journal Name},
  volume  = {12}, number = {3}, pages = {100--120}, year = {2027},
  doi     = {10.1000/xyz},
  eprint  = {2701.00001},
  links   = {PDF|/papers/doe2027.pdf ; BibTeX|/bibtex/doe2027.bib},
}
```

- `doi` becomes a DOI button, `eprint` (the arXiv id) becomes an arXiv button, and `url` becomes a Journal button.
- `links` adds any other button, written `Label|url` and separated by `;` (PDF, BibTeX, Code, Slides...).
- `venue = {<i>Journal</i> 12(3), 100-120 (2027)}` overrides the auto-built venue text. Use it for "Submitted (2027)", conference names and so on.
- Put PDFs in `papers/` and `.bib` files in `bibtex/`.
- "G. B. Wright" is bolded automatically in author lists.
- Escape ampersands as `\&`. A new year gets its own heading and a link in the year bar automatically.

## Talks and software (Research page)

Open `research.html` and scroll to the bottom. Talks and software are short lists. Copy an existing `<li>` and edit the text and link.

## Teaching

`_data/teaching.yml`. Newest term first. Copy a block and edit it:

```yaml
- term: Fall 2027
  courses:
  - Numerical Analysis I (<a href="https://math.boisestate.edu/~wright/courses/m565/">MATH 565</a>)
```

Use `note:` instead of `courses:` for things like a sabbatical (see Fall 2013 and Spring 2014). The sections are "Previous courses", Utah, and Colorado. The "I am not teaching at the moment" sentence at the top is `intro:` in the same file. Change or delete it when teaching resumes.

## Students

`_data/students.yml` has two top-level sections, Current students and Former students, each split into Ph.D., M.S. and B.S. groups. To add a student, copy an entry into the right group:

```yaml
- name: Jane Doe
  url: "https://example.com/jane"          # optional, links the name
  photo: jane-doe                          # optional, see below
  details:
    - 'Ph.D. in <a href="https://www.boisestate.edu/computing/emphasis/cmse/">Computational Math, Science &amp; Engineering</a>, 2028'
    - 'Thesis: <a href="https://scholarworks.boisestate.edu/">Title</a>'
    - 'Currently a Postdoc at Somewhere'
```

Each line under `details` is one line on the page. HTML links are fine, and `&` should be written `&amp;`. When a current student graduates, move their entry from Current to Former and update the details.

**Photos:** save a square image (about 256 x 256 pixels) as `images/students/jane-doe.jpg` and also as `images/students/jane-doe.webp`, then set `photo: jane-doe`. The webp file is smaller, and the jpg is the fallback. Only the current Ph.D. group and former Ph.D. entries use photos right now, but any entry can have one.

## Menu links

`_data/nav.yml`. Each item is a `title` and a `url`. Add `external: true` to make a link open in a new tab (used for Software and Google Scholar). A new page, for example `service.html`, needs front matter at the top:

```
---
layout: page
title: Service
---
```

Then add it to `nav.yml` with `url: /service.html`.

## Home page

`index.md` is plain HTML with small blocks: hero (name, role line, buttons), About, Research, Teaching. The role line under the name is where his current title lives. Edit text directly and keep the surrounding tags. The animated sphere is a `<canvas data-sphere>` and needs no changes.

## Files

- **CV:** replace `documents/GradyWrightCV.pdf` with the new file, keeping the same name. The CV link in the menu then needs no change.
- **Papers, BibTeX, other downloads:** `papers/`, `bibtex/`, `documents/`. Link to them as `/papers/filename.pdf`.
- **Old course pages** (`courses/`, `montestigliano/`, `talks/`) are untouched. They are plain static files.

## Design notes

- Colors and fonts are variables at the top of `assets/css/style.css` (`--navy`, `--accent`, `--paper` and so on). Change them once and the whole site follows.
- Fonts (Source Serif 4, Source Sans 3) are included in `assets/fonts/`, so there are no external requests.
- Animations stop automatically for visitors who set "reduce motion" in their system settings.
- The site works without JavaScript. The menu and the sphere are the only things that need it, and the pages stay readable without them.

## If something breaks

- **A page shows raw `{{ }}` text or fails to build:** open the repository's *Actions* tab on GitHub and read the failed "pages build" for the line number. The usual causes are a missing space in a data file (YAML indents with spaces, not tabs) or an unquoted colon in a title (wrap the value in single quotes).
- **Publications didn't update:** you edited `publications.bib` but didn't run `python3 scripts/bib2yml.py` and commit the regenerated `publications.yml`.
- **Changes not appearing:** wait two minutes, then hard-refresh the browser (Cmd+Shift+R or Ctrl+Shift+R).
