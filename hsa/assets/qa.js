/* Ask the atlas: plain-English questions answered from the sample table in the page.

   There is no model call here and no API key. The page is a static, encrypted
   artefact on GitHub Pages, so anything that reached an API would have to ship a
   key to every reader. What this does instead is resolve the words in a question
   against the catalogue's own vocabulary -- every tissue, assay, disease and
   donor attribute it actually holds -- and then count.

   The rule the answers are built around: a blank field is not a negative.
   Coverage in this catalogue runs from 100% (assay) to 21% (age group), so a
   count over age looks exactly like a count over assay unless the answer says
   how many samples could have answered at all. Every answer states that, and a
   question that leans on a sparse field says so first.

   A term that resolves to nothing is reported as unresolved, never silently
   dropped, because "0 samples" from a typo and "0 samples" from the data are
   the same number on screen and opposite in meaning. */
(function (global) {
  'use strict';

  var D = null, INDEX = null;

  /* Words people use that the catalogue spells differently. Kept small on
     purpose: everything else resolves against the catalogue's own values, so
     this list only has to cover the cases where usage and vocabulary genuinely
     part company. */
  var SYN = {
    'tumour': 'cancer', 'tumor': 'cancer', 'tumours': 'cancer', 'tumors': 'cancer',
    'carcinoma': 'cancer', 'malignant': 'cancer', 'oncology': 'cancer',
    'healthy': 'normal', 'control': 'normal', 'controls': 'normal', 'non diseased': 'normal',
    'woman': 'female', 'women': 'female', 'man': 'male', 'men': 'male',
    'bowel': 'intestine', 'gut': 'intestine', 'colon': 'colon', 'renal': 'kidney',
    'hepatic': 'liver', 'cardiac': 'heart', 'cerebral': 'brain', 'pulmonary': 'lung',
    'gastric': 'stomach', 'cutaneous': 'skin', 'dermal': 'skin',
    'visiumhd': 'visium hd', 'visium-hd': 'visium hd', 'merscope': 'merfish',
    'stereoseq': 'stereo-seq', 'slideseq': 'slide-seq', 'seqfish': 'seqfish',
    'paediatric': 'child', 'pediatric': 'child', 'elderly': 'aged',
    'ffpe': 'ffpe', 'fresh frozen': 'fresh_frozen', 'frozen': 'fresh_frozen'
  };

  /* Words that carry intent, not content. Stripped before term matching so
     "how many lung samples" does not try to match "many" against a tissue. */
  var STOP = ('how many much number count of the a an is are there do we have '
    + 'show me list which what where who all and or in for with from to by per '
    + 'samples sample dataset datasets study studies data atlas hsa please tell '
    + 'about got holds has does any give').split(' ');

  /* Words that name a FIELD. They must never be matched as a value: "which
     tissues do we have for Visium HD" was resolving "tissues" against the
     tissue value "multiple_tissues_19_cancer_types" and returning 0, because
     the question's own grammar had become a filter. */
  var FIELD_WORDS = {
    tissue: 'tissue', tissues: 'tissue', organ: 'tissue', organs: 'tissue',
    assay: 'assay', assays: 'assay', technology: 'assay', technologies: 'assay',
    platform: 'assay', platforms: 'assay',
    disease: 'disease_level1', diseases: 'disease_level1', condition: 'disease_level1',
    conditions: 'disease_level1', sex: 'sex', gender: 'sex', age: 'age_group',
    study: 'accession', studies: 'accession', source: 'source', sources: 'source',
    donor: 'donor_id', donors: 'donor_id', preservation: 'preservation'
  };

  function norm(s) {
    return String(s || '').toLowerCase()
      .replace(/[^a-z0-9+\- ]+/g, ' ').replace(/\s+/g, ' ').trim();
  }

  /* One lookup from a normalised string to every (field, value) that string can
     mean. Built once. A value is indexed under its whole normalised form and
     under each of its words, so "breast" finds "breast" and also
     "breast cancer tissue". */
  function buildIndex() {
    INDEX = {};
    var add = function (key, field, vi, exact) {
      if (!key || key.length < 2) return;
      (INDEX[key] = INDEX[key] || []).push({ f: field, i: vi, exact: !!exact });
    };
    Object.keys(D.cols).forEach(function (f) {
      var vals = D.cols[f].v;
      for (var i = 1; i < vals.length; i++) {         // 0 is the blank
        var n = norm(vals[i]);
        add(n, f, i, true);
        n.split(' ').forEach(function (w) { if (STOP.indexOf(w) < 0) add(w, f, i, false); });
      }
    });
  }

  /* Longest-match term resolution over the question. Returns the matched terms
     and, separately, the words that matched nothing, which the answer reports. */
  function resolve(q) {
    var words = norm(q).split(' ').map(function (w) { return SYN[w] || w; }).filter(Boolean);
    var hits = [], unresolved = [], notes = [], used = new Array(words.length).fill(false);
    for (var len = Math.min(4, words.length); len >= 1; len--) {
      for (var i = 0; i + len <= words.length; i++) {
        if (used.slice(i, i + len).some(Boolean)) continue;
        var phrase = words.slice(i, i + len).join(' ');
        phrase = SYN[phrase] || phrase;
        if (len === 1 && FIELD_WORDS[phrase]) continue;   // a field name, not a value
        var m = INDEX[phrase];
        if (!m || !m.length) continue;
        /* A term must not match a value that excludes it. The catalogue has
           "astrocytoma_excluding_glioblastoma", and matching "glioblastoma"
           into it returns the one category that is defined by not being it. */
        m = m.filter(function (x) {
          if (x.exact) return true;
          var v = norm(D.cols[x.f].v[x.i]);
          return !(new RegExp('(excluding|except|other than|non|not)\\b[a-z ]{0,20}\\b'
                   + phrase.replace(/[.*+?^${}()|[\]\\]/g, '\\$&')).test(v));
        });
        if (!m.length) {
          /* The only values this term matched are ones defined by NOT being it.
             That is a real statement about the catalogue's vocabulary and the
             reader needs it, so it is reported rather than dropped. */
          var excl = INDEX[phrase].map(function (x) { return D.cols[x.f].v[x.i]; })
            .filter(function (x, k, a) { return a.indexOf(x) === k; }).slice(0, 3);
          notes.push('"' + phrase + '" is not a category here. The only values containing '
            + 'that word are defined by excluding it: ' + excl.join(', ') + '.');
          for (var z = i; z < i + len; z++) used[z] = true;
          continue;
        }
        // prefer an exact whole-value match, and the field with fewest values
        // (an assay name is a stronger signal than one word of a tissue label)
        /* An accession is a URL or a GEO id, so any word can appear inside one.
           A match there is the weakest kind and must never outrank a real
           content field: "glioblastoma" meant the disease, not a study whose
           URL happens to spell it. */
        m = m.slice().sort(function (a, b) {
          var aw = a.f === 'accession' ? 1 : 0, bw = b.f === 'accession' ? 1 : 0;
          if (aw !== bw) return aw - bw;
          if (a.exact !== b.exact) return a.exact ? -1 : 1;
          return D.cols[a.f].v.length - D.cols[b.f].v.length;
        });
        if (m[0].f === 'accession' && !m[0].exact) {
          notes.push('"' + phrase + '" matched only inside a study identifier, not any '
            + 'tissue, assay or disease category.');
          for (var z2 = i; z2 < i + len; z2++) used[z2] = true;
          continue;
        }
        var best = m[0];
        // keep every value in the same field that this phrase matched, so
        // "cancer" selects all cancer categories rather than an arbitrary one
        var sameField = m.filter(function (x) { return x.f === best.f; });
        hits.push({ term: phrase, field: best.f, values: sameField.map(function (x) { return x.i; }) });
        for (var k = i; k < i + len; k++) used[k] = true;
      }
    }
    words.forEach(function (w, i) {
      if (!used[i] && STOP.indexOf(w) < 0 && !FIELD_WORDS[w] && w.length > 2) unresolved.push(w);
    });
    return { hits: hits, unresolved: unresolved, notes: notes, words: words };
  }

  /* Which field, if any, the question wants broken down. */
  function groupField(words) {
    var byIdx = words.indexOf('by');
    var want = null;
    if (byIdx >= 0 && byIdx + 1 < words.length) want = FIELD_WORDS[words[byIdx + 1]] || null;
    if (!want) {
      var q = words.join(' ');
      if (/^(which|what)\b/.test(q)) {
        for (var w in FIELD_WORDS) if (words.indexOf(w) >= 0) { want = FIELD_WORDS[w]; break; }
      }
    }
    return want;
  }

  function applyFilters(hits) {
    var n = D.n, keep = new Uint8Array(n).fill(1);
    hits.forEach(function (h) {
      var codes = D.cols[h.field].c, set = {};
      h.values.forEach(function (v) { set[v] = 1; });
      for (var i = 0; i < n; i++) if (keep[i] && !set[codes[i]]) keep[i] = 0;
    });
    return keep;
  }

  function tally(keep, field, limit) {
    var codes = D.cols[field].c, vals = D.cols[field].v, counts = {};
    var known = 0, blank = 0;
    for (var i = 0; i < D.n; i++) {
      if (!keep[i]) continue;
      var c = codes[i];
      if (!c) { blank++; continue; }
      known++; counts[c] = (counts[c] || 0) + 1;
    }
    var rows = Object.keys(counts).map(function (c) { return [vals[c], counts[c]]; })
      .sort(function (a, b) { return b[1] - a[1]; });
    return { rows: limit ? rows.slice(0, limit) : rows, known: known, blank: blank,
             distinct: rows.length };
  }

  function count(keep) { var t = 0; for (var i = 0; i < D.n; i++) t += keep[i]; return t; }

  function pct(a, b) { return b ? (100 * a / b).toFixed(a / b > 0.995 ? 0 : 1) + '%' : '0%'; }

  function answer(q) {
    if (!D) return { error: 'The sample table is not loaded on this page.' };
    var r = resolve(q);
    var group = groupField(r.words);
    var keep = applyFilters(r.hits);
    var n = count(keep);

    var filters = r.hits.map(function (h) {
      var col = D.cols[h.field];
      return {
        term: h.term, field: h.field, label: col.label,
        values: h.values.map(function (v) { return col.v[v]; })
          .filter(function (x, i, a) { return a.indexOf(x) === i; }),
        coverage: D.coverage[h.field], coveragePct: pct(D.coverage[h.field], D.n)
      };
    });

    /* The denominator that matters is not 13,652 but how many samples could
       have answered this question at all: a sample with no sex recorded is not
       a male sample, and it is not a female one either. */
    var answerable = D.n;
    filters.forEach(function (f) { answerable = Math.min(answerable, f.coverage); });
    var sparse = filters.filter(function (f) { return f.coverage / D.n < 0.6; });

    /* "how many female DONORS" is not "how many female samples". When the
       question names a unit the catalogue can count distinctly, count that. */
    var unit = 'samples', unitN = n;
    ['donor_id', 'accession'].forEach(function (f) {
      var w = Object.keys(FIELD_WORDS).filter(function (k) { return FIELD_WORDS[k] === f; });
      if (!D.cols[f] || !w.some(function (k) { return r.words.indexOf(k) >= 0; })) return;
      var seen = {}, codes = D.cols[f].c, blank = 0;
      for (var i = 0; i < D.n; i++) {
        if (!keep[i]) continue;
        if (!codes[i]) { blank++; continue; }
        seen[codes[i]] = 1;
      }
      unit = (f === 'accession' ? 'studies' : D.cols[f].label + 's');
      unitN = Object.keys(seen).length;
      out_unitBlank = blank;
    });
    var out_unitBlank = typeof out_unitBlank === 'undefined' ? 0 : out_unitBlank;

    var out = {
      question: q, n: n, total: D.n, answerable: answerable,
      unit: unit, unitN: unitN, unitBlank: out_unitBlank,
      filters: filters, unresolved: r.unresolved, notes: r.notes, sparse: sparse, group: null,
      studies: null
    };
    if (group && D.cols[group]) out.group = Object.assign(
      { field: group, label: D.cols[group].label }, tally(keep, group, 25));
    if (n && n <= D.n) {
      var st = tally(keep, 'accession', 8);
      out.studies = st.rows;
      out.nStudies = st.distinct;
    }
    return out;
  }

  global.HSAQuery = {
    load: function (payload) { D = payload; buildIndex(); return !!D; },
    ask: answer,
    fields: function () { return D ? Object.keys(D.cols) : []; },
    coverage: function () { return D ? D.coverage : {}; },
    values: function (f) { return D && D.cols[f] ? D.cols[f].v : []; }
  };
})(typeof window !== 'undefined' ? window : globalThis);
