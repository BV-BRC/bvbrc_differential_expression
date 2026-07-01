/*
 * expression-mapping-preview.sketch.js
 *
 * Browser-side sketch of the gene-ID -> feature_id mapping performed by
 * expression_transform.py (make_map_query + place_ids), intended as a live
 * "preview" in the Differential Expression GUI so the user sees how many of
 * their gene IDs resolve *before* submitting a job.
 *
 * Key differences from the current Python (these are the fixes):
 *   - No hard `AND annotation:PATRIC` filter. We request `annotation` in the
 *     field list and, when one source ID resolves to multiple features, PREFER
 *     PATRIC, else fall back to RefSeq. This is what lets RefSeq-only locus tags
 *     (e.g. E. coli b#### numbers) resolve.
 *   - Accurate result: returns mapped/unmapped counts + lists, instead of a
 *     hard exit on zero matches. Caller decides how to surface it.
 *
 * Same-origin note: the GUI is served from www.bv-brc.org and the data API is
 * at www.bv-brc.org/api, so this is a same-origin POST (no CORS). Pass the
 * user's existing auth token string.
 *
 * No build deps for the mapping itself. File parsing (csv/tsv/xlsx) is separate
 * and would use PapaParse / SheetJS.
 */

const DEFAULT_SOURCE_TYPES = [
  "refseq_locus_tag", "alt_locus_tag", "feature_id",
  "protein_id", "patric_id", "gene",
];
const INT_TYPES = ["gi", "gene_id"];
const CHUNK_SIZE = 1000;

/** Escape a value for use inside a Solr term and quote it. */
function solrQuote(v) {
  return '"' + String(v).replace(/(["\\])/g, "\\$1") + '"';
}

function* chunk(arr, size) {
  for (let i = 0; i < arr.length; i += size) yield arr.slice(i, i + size);
}

/**
 * Build the Solr query body for one chunk of IDs.
 * If sourceIdType is given we query only that field; otherwise we OR across all
 * known source types (and numeric-only IDs across the int types).
 */
function buildQueryBody(ids, { sourceIdType, genomeId, host }) {
  const sourceTypes = sourceIdType ? [sourceIdType] : DEFAULT_SOURCE_TYPES;
  const clauses = [];

  if (!sourceIdType) {
    const intIds = ids.filter((id) => /^\d+$/.test(String(id)));
    if (intIds.length) {
      const ored = intIds.map(solrQuote).join(" OR ");
      for (const t of INT_TYPES) clauses.push(`(${t}:(${ored}))`);
    }
  }
  const ored = ids.map(solrQuote).join(" OR ");
  for (const t of sourceTypes) clauses.push(`(${t}:(${ored}))`);

  // NOTE: no `annotation:PATRIC` filter here (the fix). For host genomes the
  // Python restricts to RefSeq; keep that option available.
  let q = `(${clauses.join(" OR ")})`;
  if (host) q += " AND annotation:RefSeq";
  if (genomeId) q += ` AND genome_id:${solrQuote(genomeId)}`;

  // Ask for annotation so we can prefer PATRIC over RefSeq on collisions.
  const fl = ["feature_id", "annotation", ...new Set([...sourceTypes, ...INT_TYPES])].join(",");
  return new URLSearchParams({ q, fl, rows: "25000", wt: "json" });
}

/** Prefer a PATRIC feature over a RefSeq one for the same source ID. */
function preferAnnotation(existing, incoming) {
  if (!existing) return true;
  return existing.annotation !== "PATRIC" && incoming.annotation === "PATRIC";
}

/**
 * Map a list of gene IDs to BV-BRC feature_ids.
 *
 * @param {string[]} geneIds   unique source IDs from the expression file
 * @param {object}   opts
 * @param {string}   opts.dataApi       API base, e.g. "https://www.bv-brc.org/api"
 * @param {string}   [opts.sourceIdType] refseq_locus_tag | alt_locus_tag | ...
 * @param {string}   [opts.genomeId]     restrict to one genome (optional)
 * @param {boolean}  [opts.host]         host genome: identity mapping, no API
 * @param {string}   [opts.token]        BV-BRC auth token string
 * @returns {Promise<{idMap:Map, total:number, mapped:number, unmapped:number,
 *                     unmappedList:string[]}>}
 */
export async function mapGeneIds(geneIds, opts) {
  const ids = [...new Set(geneIds.map(String))].filter((s) => s && s !== "-");
  const total = ids.length;

  // Host genome: identity mapping, no network (mirrors --host).
  if (opts.host) {
    const idMap = new Map(ids.map((id) => [id, id]));
    return { idMap, total, mapped: total, unmapped: 0, unmappedList: [] };
  }

  const base = opts.dataApi.replace(/\/+$/, "");
  const url = `${base}/genome_feature/`;
  const headers = {
    "Content-Type": "application/solrquery+x-www-form-urlencoded",
    Accept: "application/solr+json",
  };
  if (opts.token) headers.Authorization = opts.token;

  // sourceId -> {feature_id, annotation}, applying PATRIC-preference.
  const best = new Map();

  for (const ch of chunk(ids, CHUNK_SIZE)) {
    const body = buildQueryBody(ch, opts);
    let resp;
    try {
      resp = await fetch(url, { method: "POST", headers, body });
    } catch (e) {
      throw new Error(`data API request failed: ${e.message}`);
    }
    if (!resp.ok) {
      throw new Error(`data API error ${resp.status}: ${await resp.text()}`);
    }
    const docs = (await resp.json())?.response?.docs ?? [];
    for (const d of docs) {
      if (!d.feature_id) continue;
      const cand = { feature_id: d.feature_id, annotation: d.annotation };
      for (const t of [...DEFAULT_SOURCE_TYPES, ...INT_TYPES]) {
        const val = d[t];
        if (val == null) continue;
        for (const sv of Array.isArray(val) ? val : [val]) {
          const key = String(sv);
          if (preferAnnotation(best.get(key), cand)) best.set(key, cand);
        }
      }
    }
  }

  const idMap = new Map();
  const unmappedList = [];
  for (const id of ids) {
    const hit = best.get(id);
    if (hit) idMap.set(id, hit.feature_id);
    else unmappedList.push(id);
  }
  const mapped = idMap.size;
  return { idMap, total, mapped, unmapped: total - mapped, unmappedList };
}

/*
 * Example GUI usage (preview before submit):
 *
 *   const { total, mapped, unmapped, unmappedList } =
 *     await mapGeneIds(geneIdsFromFile, {
 *       dataApi: "https://www.bv-brc.org/api",
 *       sourceIdType: "refseq_locus_tag",
 *       genomeId: selectedGenomeId,        // from the genome picker
 *       token: window.App.authorizationToken,
 *     });
 *
 *   if (mapped === 0) {
 *     showError(`0 of ${total} gene IDs matched in genome ${selectedGenomeId}. ` +
 *               `Check that the genome and ID type match your data.`);
 *   } else {
 *     showInfo(`${mapped} of ${total} gene IDs mapped ` +
 *              `(${unmapped} unmapped).`);
 *   }
 *
 * This preview would have caught the bug.json case immediately: genome
 * 83333.1297 -> 0 mapped, surfaced at genome-selection time instead of as a
 * failed job.
 */
