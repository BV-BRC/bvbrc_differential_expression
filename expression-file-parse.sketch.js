/*
 * expression-file-parse.sketch.js
 *
 * Browser-side sketch of the file reading + layout detection performed by
 * expression_transform.py (process_table + fix_headers + gene_matrix_to_list).
 * Turns a user-uploaded comparisons file (csv/tsv/xls/xlsx) into the canonical
 * long form used by the mapping preview:
 *
 *     [{ exp_locus_tag, sampleUserGivenId, log_ratio }, ...]
 *
 * Pairs with expression-mapping-preview.sketch.js: feed `result.geneIds` into
 * mapGeneIds() to show mapped/unmapped counts before submit.
 *
 * Deps (both standard, browser-ready):
 *   - PapaParse  -> csv/tsv                (https://www.papaparse.com/)
 *   - SheetJS    -> xls/xlsx (import XLSX) (https://sheetjs.com/)
 * Import/inject however the GUI bundles third-party libs.
 */

const CLAMP = 1_000_000;

// Layout detection: a "gene_list" file has exactly these (normalized) columns;
// anything else with a leading gene column is treated as a "gene_matrix".
const LIST_COLUMNS = ["gene_id", "comparison_id", "log_ratio"];

/** Normalize a header the way fix_headers does. */
function fixName(x, allColumns) {
  let n = String(x).trim().replace(/\s+/g, " ").toLowerCase().replace(/ /g, "_");
  // PATRIC's downloadable template isn't consistent about plurals.
  if (n.endsWith("s") && allColumns.has(n.slice(0, -1))) n = n.slice(0, -1);
  return n;
}

/** Read the raw table as {headers:[...], rows:[{header:value}]}, order preserved. */
async function readRaw(file) {
  const ext = (file.name.split(".").pop() || "").toLowerCase();

  if (ext === "xls" || ext === "xlsx") {
    const XLSX = await import("xlsx"); // or use the globally-injected XLSX
    const wb = XLSX.read(await file.arrayBuffer(), { type: "array" });
    const ws = wb.Sheets[wb.SheetNames[0]];
    const aoa = XLSX.utils.sheet_to_json(ws, { header: 1, blankrows: false, defval: "" });
    const headers = (aoa.shift() || []).map(String);
    const rows = aoa.map((r) => Object.fromEntries(headers.map((h, i) => [h, r[i] ?? ""])));
    return { headers, rows };
  }

  // csv/tsv (or unknown -> let PapaParse sniff the delimiter)
  const Papa = await import("papaparse"); // or the globally-injected Papa
  return new Promise((resolve, reject) => {
    Papa.parse(file, {
      header: true,
      skipEmptyLines: true,
      dynamicTyping: false,           // keep strings; we coerce log_ratio ourselves
      delimiter: ext === "tsv" ? "\t" : ext === "csv" ? "," : "", // "" = auto-sniff
      complete: (res) => resolve({ headers: res.meta.fields || [], rows: res.data }),
      error: reject,
    });
  });
}

function toNumberOrNaN(v) {
  if (v === "" || v == null) return NaN;
  const n = Number(v);
  return Number.isFinite(n) ? n : NaN;
}

function clampDrop(records) {
  const out = [];
  for (const r of records) {
    if (!r.exp_locus_tag || r.exp_locus_tag === "-") continue;
    let lr = toNumberOrNaN(r.log_ratio);
    if (Number.isNaN(lr)) continue;                 // drop non-numeric / missing
    if (lr > CLAMP) lr = CLAMP;
    if (lr < -CLAMP) lr = -CLAMP;
    out.push({ ...r, log_ratio: lr });
  }
  return out;
}

/**
 * Parse an uploaded expression comparisons file into canonical long form.
 *
 * @param {File} file  the user's uploaded csv/tsv/xls/xlsx
 * @returns {Promise<{ setup:"gene_list"|"gene_matrix",
 *                      rows: Array<{exp_locus_tag,sampleUserGivenId,log_ratio:number}>,
 *                      geneIds: string[] }>}   geneIds = unique exp_locus_tag
 */
export async function parseExpressionFile(file) {
  const { headers, rows } = await readRaw(file);
  if (!headers.length) throw new Error("no columns found in expression file");

  // Python assumes the first column is the gene id regardless of its label.
  const geneCol = headers[0];
  const allColumns = new Set([...LIST_COLUMNS]);
  const normalized = headers.map((h) => fixName(h, allColumns));

  // gene_list iff every column normalizes into the list column set.
  const isList = normalized.every((n) => LIST_COLUMNS.includes(n));

  let records;
  if (isList) {
    const idx = Object.fromEntries(normalized.map((n, i) => [n, headers[i]]));
    records = rows.map((r) => ({
      exp_locus_tag: String(r[idx.gene_id] ?? "").trim(),
      sampleUserGivenId: String(r[idx.comparison_id] ?? ""),
      log_ratio: r[idx.log_ratio],
    }));
    return finalize("gene_list", records);
  }

  // gene_matrix: first column is the gene, every other column is a sample.
  // Melt wide -> long (equivalent of pd.melt / gene_matrix_to_list).
  const sampleCols = headers.slice(1);
  records = [];
  for (const r of rows) {
    const gene = String(r[geneCol] ?? "").trim();
    if (!gene) continue;
    for (const s of sampleCols) {
      records.push({ exp_locus_tag: gene, sampleUserGivenId: s, log_ratio: r[s] });
    }
  }
  return finalize("gene_matrix", records);
}

function finalize(setup, records) {
  const rows = clampDrop(records);
  if (!rows.length) throw new Error("no usable rows after parsing expression file");
  const geneIds = [...new Set(rows.map((r) => r.exp_locus_tag))];
  return { setup, rows, geneIds };
}

/*
 * End-to-end GUI preview (parse -> map), no submit:
 *
 *   const { setup, rows, geneIds } = await parseExpressionFile(uploadedFile);
 *   const { total, mapped, unmapped } = await mapGeneIds(geneIds, {
 *     dataApi: "https://www.bv-brc.org/api",
 *     sourceIdType: selectedSourceIdType,   // or omit to auto-detect across types
 *     genomeId: selectedGenomeId,
 *     token: window.App.authorizationToken,
 *   });
 *   showInfo(`${setup}: ${rows.length} data points, ` +
 *            `${mapped}/${total} genes mapped (${unmapped} unmapped).`);
 *
 * Not covered here (still needed for a full client-side transform, if that path
 * is chosen): per-sample stats (mean/std/z-score/significance counts), the four
 * output JSONs, and triggering indexing. See PLAN-expression-transform-modernization.md.
 */
