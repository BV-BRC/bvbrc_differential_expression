# Plan: Modernize & fix `bvbrc_differential_expression/expression_transform.py`

## Context

The Differential Expression import tool (`expression_transform.py`, called by
`App-DifferentialExpression.pl`) has three problems we want to fix together:

1. **A real mapping failure.** For the `bug.json` run the gene-ID→`feature_id`
   mapping returns zero hits and the job dies with a misleading "PATRIC API is
   down" error. Verified against the live API, two independent causes: (a) the
   `genome_id` in the params (`83333.1297`, substrain BOP27) does not carry the
   `b####` locus tags in the file — they belong to `511145.12`/`83333.111`; and
   (b) the query hard-codes `AND annotation:PATRIC`, but `b####` refseq_locus_tags
   exist **only on `annotation:RefSeq` features**, so even the correct genome
   returns nothing under the PATRIC filter.
2. **A clumsy interface.** The Perl driver writes two throwaway JSON files
   (`ustr.$$`, `sstr.$$`) that the Python re-parses. `sstr` is just a `data_api`
   URL; `ustr` is a flat bag of scalars. Both should be plain CLI flags.
3. **Legacy code.** Python-2 vestiges, deprecated pandas calls (pandas 2.0.3 /
   Python 3.10.13 in the runtime), dead functions, and debug prints.

Outcome: a clean flag-based interface (no JSON temp files), a mapping step that
actually resolves RefSeq-style locus tags and reports failures accurately, and a
fully Python-3 / modern-pandas-clean program.

Decision (confirmed): **clean break** — remove `--ufile/--ustring/--sfile/--sstring`
entirely; the sole caller (`App-DifferentialExpression.pl`) is updated in lockstep.

## Files

- `bvbrc_differential_expression/expression_transform.py` — main rewrite.
- `bvbrc_differential_expression/service-scripts/App-DifferentialExpression.pl` — driver.
- `bvbrc_differential_expression/requirements.txt` — confirm only (already lists
  pandas/numpy/scipy/openpyxl/requests); no change expected.
- `bvbrc_differential_expression/test/expression_transform.py` — check whether it
  is a stale duplicate; update or delete to match, or leave if unrelated.
- App spec (`app_specs/DifferentialExpression.json`) — **unchanged**; it still
  sends `ustring`, which the Perl now decodes into flags. (Replacing the `ustring`
  blob param with typed app-spec params is a separate, out-of-scope service change.)

## Part A — Flag-based interface (no JSON files)

**New argparse interface in `main()`** (replaces the JSON args):

```
--xfile FILE            (required)   comparisons file
--mfile FILE            (optional)   metadata template
--output-path DIR       (required)
--data-api URL          (required)   API base, e.g. https://www.bv-brc.org/api
--source-id-type TYPE   (required)   refseq_locus_tag | alt_locus_tag | feature_id | ...
--data-type TYPE        (required)   Transcriptomics | Proteomics | Phenomics
--title STR             (required)   was experiment_title
--description STR       (required)   was experiment_description
--organism STR          (optional)
--genome-id STR         (optional)
--pmid STR              (optional)
--host                  (flag)       host genome: identity mapping, no API calls
```

- Delete `--ufile/--ustring/--sfile/--sstring`, both JSON `try/except` parse
  blocks, and the hand-rolled `req_info`/`missing` presence check (argparse
  `required=True` now enforces the four required scalars).
- Keep the downstream dict-based signatures; `main()` builds them from args:
  `form_data` = the scalar fields + the injected constants `source_types`
  (`refseq_locus_tag, alt_locus_tag, feature_id, protein_id, patric_id, gene`) and
  `int_types` (`gi, gene_id`); `server_setup = {"data_api": args.data_api}`.
- **Drop `xsetup`** — it is never read (layout is auto-detected in `fix_headers`).
- **Single `--host`** drives both the mapping-skip and the RefSeq switch (today
  `host` is duplicated in the flag and in `ustring.host`).
- `--data-api` is the **base** URL; the code appends `genome_feature/` where it
  POSTs (today the Perl pre-appends it). This also lets
  `diffexp_api.getGenomeIdsNamesByName` reuse the base.

**`App-DifferentialExpression.pl`:** remove the `sstr.$$`/`ustr.$$` file handles and
their `open/print/close`; `decode_json($ustring)` once and pass flags through
`IPC::Run` (array form — already used, no shell quoting):

```perl
my $u = decode_json($ustring);
my @cmd = ("expression_transform",
    "--xfile", $xfile_tmp, @mfile_arg,
    "--output-path", $out,
    "--data-api", Bio::KBase::AppService::AppConfig->data_api_url,  # base, no /genome_feature/
    "--source-id-type", $u->{source_id_type},
    "--data-type", $u->{data_type},
    "--title", $u->{experiment_title},
    "--description", $u->{experiment_description},
    (length($u->{organism}  // "") ? ("--organism",  $u->{organism})  : ()),
    (length($u->{genome_id} // "") ? ("--genome-id", $u->{genome_id}) : ()),
    (length($u->{pmid}      // "") ? ("--pmid",      $u->{pmid})      : ()),
    ($u->{host} ? ("--host") : ()));
```

## Part B — Mapping correctness (`make_map_query` / `place_ids` / `map_gene_ids`)

1. **Don't hard-require `annotation:PATRIC`.** Drop the annotation clause from the
   query (keep the optional `genome_id` filter). Retrieve all matching docs; when a
   source id maps to multiple `feature_id`s, **prefer `annotation:PATRIC`, else
   `RefSeq`**. This is what makes the `b####` RefSeq-only tags resolve. Preserve the
   `--host` path (identity mapping, no API).
2. **Vectorize the mapping.** Replace the per-row loops + chained assignment +
   Series-`in` bug with: build `id_map = {source_id: feature_id}` from the response
   docs (iterate the `source_types`/`int_types` fields present per doc, applying the
   PATRIC-preference), then
   `mapping_table["feature_id"] = mapping_table["exp_locus_tag"].astype(str).map(id_map)`.
   Removes the `SettingWithCopy` chained writes and the index-vs-value confusion.
   `--host` branch becomes `mapping_table["feature_id"] = mapping_table["exp_locus_tag"]`.
3. **Unify `source_types`** between `make_map_query` and `place_ids` (today they
   differ). Query by `--source-id-type` when given; always request `feature_id` +
   annotation in `fl`; de-dup the `fl` list.
4. **Accurate diagnostics.** Only `exit(2)` on a genuine API error (HTTP not ok).
   When 0 of N map, emit e.g.
   `"Mapped 0 of N gene IDs (genome_id=<id>, source_id_type=<type>) — check that the
   genome and ID type match your data"` instead of the "PATRIC API is down" text.
   Report mapped/total counts on success too.
5. Fix `sys.stderr.write(e)` → `sys.stderr.write(str(e))` in the error path.

## Part C — Python 3 + modern pandas (pandas 2.0.3)

- Remove the `sys.version_info < (2,7): raise "must use…"` block (invalid `raise
  "<string>"` under py3), the `#requires 2.7.9` note, the `FutureWarning`
  suppression, and `pd.options.mode.chained_assignment = None` (no longer needed
  once chained assignment is gone).
- `pd.read_table(...)` → `pd.read_csv(..., sep='\t')` (read_table is deprecated).
- `pd.read_excel(target_file, 0, index_col=None)` → `sheet_name=0` (keyword).
- `grouped['log_ratio'].agg([np.mean, np.std])` → `.agg(['mean','std'])` (numpy
  funcs in agg are deprecated); the existing rename→`expmean`/`expstddev` still applies.
- `groupby(["sampleUserGivenId"], ...)` single-key list → `groupby("sampleUserGivenId")`
  (avoids the tuple-key FutureWarning).
- Replace the always-true `'mfile' in map_args` / `'sstring' in map_args` membership
  tests with truthiness (`map_args.mfile`); the `sstring` one disappears with the clean break.
- Keep `numpy` import (`np.nan`, `np.issubdtype`) and `scipy.stats.zscore`; verify
  no new warnings from `transform(stats.zscore)`.

## Part D — Dead code / noise removal

- Delete `gene_list_to_matrix` (never called) and `pretty_print_POST` (only in a
  commented line).
- Simplify `process_table`: drop the `tries`/`next_up` retry machinery (its
  recursive fallback is commented out at lines ~176–178) → detect format (extension
  or `csv.Sniffer`) → read once → clear error on failure.
- Remove unused `user_parse` / `server_parse` / `parse_server` in `main`.
- Remove debug output: `print(query_results.json())`, `print("query:", …)`, the
  "Grouped DataFrame structure" dumps, and the `try/except` around `.agg` that only
  re-raises.

## Verification

Runtime: `source /home/olson/P3/dev-ubuntu/user-env.sh` (Python 3.10.13, pandas
2.0.3). Auth via `KB_AUTH_TOKEN` env or `~/.patric_token`.

1. **Mapping fix + new interface, happy path** (corrected genome):
   ```
   cd modules/bvbrc_differential_expression
   python3 expression_transform.py \
     --xfile Sangurdekar_Microarray_data_all_conditions.txt \
     --output-path out.dir --data-api https://www.bv-brc.org/api \
     --source-id-type refseq_locus_tag --data-type Transcriptomics \
     --title t --description d --genome-id 511145.12
   python3 -c "import json;print(json.load(open('out.dir/mapping.json'))['mapping']['mapped_ids'])"
   ```
   Expect `mapped_ids > 0` (RefSeq fallback resolves `b####`), and all four JSON
   outputs (`experiment/expression/mapping/sample.json`) present and valid.
2. **Accurate error path:** rerun with `--genome-id 83333.1297` (the bug value) →
   expect the new "Mapped 0 of N … check genome/ID type" message and exit 2, not
   the "API is down" text.
3. **No deprecation warnings:** run step 1 with `PYTHONWARNINGS=error::FutureWarning`
   (and `::DeprecationWarning`) and confirm a clean run.
4. **Perl driver:** `perl -c service-scripts/App-DifferentialExpression.pl` (with
   runtime `PERL5LIB`) syntax-OK; confirm no `sstr.*`/`ustr.*` files are created by
   a dry component run (grep the working dir).
5. **Wrapper build:** `make` in the module builds the `expression_transform`
   wrapper via `wrap_python3`; run the wrapper once to confirm arg passing.

## Commit

Single repo `bvbrc_differential_expression` (branch `master`). On confirmation,
commit `expression_transform.py` + `App-DifferentialExpression.pl` (+ requirements/
test file only if changed), directly to `master`, left unpushed. Leave the untracked
scratch files (`bug.json`, `ustr.*`, `sstr.*`, `diffexp-mapping-analysis.slack.md`)
out of the commit.

## Out of scope
- Upstream organism→genome_id selection (web UI/service) that picked the wrong
  `genome_id` — the script-side robustness above mitigates it, but the real fix is
  in the submission layer.
- Replacing the app spec's opaque `ustring` param with typed parameters.
