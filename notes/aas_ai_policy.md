# AAS journals: AI policy and RNAAS rules (checked 2026-10-12)

## AI guidelines for authors

Source: AAS Publications Committee, "Guidelines on the use of AI / LLMs in preparing manuscripts" (authors), first posted **9 Sep 2026**; announcement https://journals.aas.org/?p=4629 ; full text https://journals.aas.org/author-llm-guidelines ; linked from https://journals.aas.org/policies/ .
- They apply to new submissions "later this fall"; no exact date is given.
- They are a living document. **Re-read them before submitting.**

| Topic | What the guidelines say |
|---|---|
| **Authorship** | AI cannot be an author (it cannot take responsibility for the content). |
| **Disclosure** | Required on every submission: state whether AI/LLMs were used, naming the tool and its purpose. It goes in an "AI disclosure statement included at the end of the article". Optionally, tell the editor which sections used it. |
| **Allowed, with disclosure and supervision** | Language editing, translation (checked by an author), accessibility, data/code formatting. Drafting text is allowed only if the authors read, check and edit all of it ("not generated unsupervised"). Explanatory illustrations: disclose the tool and check for copying. Data visualisations: check faithfulness to the data; the caption names the tool. |
| **Analysis code / data analysis** | The guidelines make no recommendation on AI use in data analysis or in writing analysis code. |
| **Never allowed** | AI content presented as data or as results from non-AI sources. Undisclosed AI use on figures, data, text or references. Plagiarism via AI. |
| **References** | AI must not be the sole source of the reference list. Check that every reference exists and is relevant (download it; use ADS). Fictitious references can lead to retraction. |
| **Responsibility** | The authors are fully accountable for all figures, data, text and references. |

**Consequences for this project**
1. Claude is not an author. Jon McBee is the sole author, so the Note uses "I" (RNAAS style for single-author Notes).
2. The Note ends with an AI disclosure statement. A draft is in `paper/dr4_companions.tex`. It must stay true:
   - Claude wrote most of the analysis code and drafts of text.
   - Jon reviews and edits every sentence and checks every reference.
   - Validation (null and injection calibration, gaiasupdate reproduction, Gaia-4b recovery) is documented in `notes/phase1_validation.md`.
3. **References:** export every BibTeX entry from ADS. Write none by hand or from memory.
4. **Figure:** if the Note uses a figure, its caption names the tool if Claude produced the plotting code. The one-table-or-one-figure limit applies.
5. **The duotransits note (submitted 2026-10-12):** check whether it carries an AI disclosure statement. If the editor asks, add one in the same form.

## RNAAS preparation rules

Source: https://journals.aas.org/research-note-preparation-guidelines

**Length**
- At most 1,500 words, counting:
  - title, headers, captions and references;
  - acknowledgments of people, facilities and software;
  - about 150 words reserved for the required abstract.
- Table contents do not count; captions do.
- The editorial office's word count is final.

**Content**
- One figure **or** one table, not both. A figure plus its caption must fit one page.
- Long tables show 5–10 lines in print; the full table is machine-readable (MRT) or on Zenodo (AAS community).
- The data behind a figure are encouraged.

**Format**
- No line numbers. AASTeX v7+ preferred; Overleaf accepted.
- Single author writes in the first person singular.

**Process**
- Moderated, not peer reviewed; typically published within 72 h.
- Published as submitted, with no copy-editing.
- DOI and ADS indexing; CC-BY 4.0; no charges.
- At most 20 Notes per first author per 12 months.
