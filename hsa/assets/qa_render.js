/* Rendering for "Ask the atlas". Kept apart from the engine so the engine can be
   tested without a DOM, which is how the three semantic bugs in it were found. */
(function () {
  'use strict';
  var EXAMPLES = ['how many Xenium samples', 'breast cancer by assay',
    'which tissues for Visium HD', 'how many FFPE samples', 'lung samples by disease',
    'how many female donors'];

  function esc(s) {
    return String(s).replace(/[&<>"]/g, function (c) {
      return { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c];
    });
  }
  var fmt = function (n) { return n.toLocaleString('en-GB'); };
  var pretty = function (s) { return String(s).replace(/_/g, ' '); };

  function render(a) {
    var out = document.getElementById('askOut');
    var h = [];

    /* No filter resolved means no answer. Showing the catalogue total here would
       read as "13,652 match", which is the opposite of what happened. */
    if (!a.filters.length && !a.group) {
      h.push('<p class="askbig">No filter<small>nothing in that question matched the '
        + 'catalogue&rsquo;s vocabulary</small></p>');
    } else {
      h.push('<p class="askbig">' + fmt(a.unitN) + ' <small>' + esc(a.unit)
        + (a.unit !== 'samples' ? ', from ' + fmt(a.n) + ' samples' : '') + '</small></p>');
      /* The denominator that matters: how many samples could have answered at
         all, which is the sparsest field the question leaned on. */
      if (a.filters.length) {
        h.push('<p class="asksub">Of ' + fmt(a.total) + ' human samples, <b>'
          + fmt(a.answerable) + '</b> have every field this question used filled in. The other '
          + fmt(a.total - a.answerable) + ' are silent on it, which is not the same as not '
          + 'matching.</p>');
      }
    }

    if (a.filters.length) {
      h.push('<div class="askfilters">' + a.filters.map(function (f) {
        return '<span>' + esc(f.term) + ' <i>&rarr;</i> ' + esc(f.label) + ': '
          + esc(f.values.slice(0, 3).map(pretty).join(', '))
          + (f.values.length > 3 ? ' <i>+' + (f.values.length - 3) + '</i>' : '')
          + ' <i>(' + f.coveragePct + ' filled)</i></span>';
      }).join('') + '</div>');
    }

    (a.sparse || []).forEach(function (f) {
      h.push('<p class="askwarn"><b>' + esc(f.label) + ' is recorded for only '
        + f.coveragePct + ' of samples</b> (' + fmt(f.coverage) + ' of ' + fmt(a.total)
        + '). This answer counts only those. The rest have no value for it, so they are '
        + 'neither counted for nor against.</p>');
    });
    (a.notes || []).forEach(function (n) { h.push('<p class="asknote">' + esc(n) + '</p>'); });
    if (a.unresolved && a.unresolved.length) {
      h.push('<p class="asknote">Not matched to anything in the catalogue: <b>'
        + esc(a.unresolved.join(', ')) + '</b>. Those words were ignored, so the count above '
        + 'does not use them.</p>');
    }

    if (a.group && a.group.rows.length) {
      var max = a.group.rows[0][1];
      h.push('<p class="asksub" style="margin-top:14px"><b>By ' + esc(a.group.label) + '</b> &mdash; '
        + fmt(a.group.distinct) + ' distinct, ' + fmt(a.group.known) + ' samples assigned'
        + (a.group.blank ? ', ' + fmt(a.group.blank) + ' with no value' : '') + '</p>');
      h.push('<div class="askgrid">' + a.group.rows.map(function (r) {
        return '<div class="askrow"><span>' + esc(pretty(r[0]))
          + '<i style="width:' + (100 * r[1] / max).toFixed(1) + '%"></i></span><b>'
          + fmt(r[1]) + '</b></div>';
      }).join('') + '</div>');
    }

    if (a.studies && a.studies.length && a.filters.length) {
      h.push('<p class="asksub" style="margin-top:14px"><b>Studies</b> &mdash; '
        + fmt(a.nStudies) + ' contribute, the largest are</p>');
      h.push('<div class="askgrid">' + a.studies.map(function (r) {
        return '<div class="askrow"><span>' + esc(r[0]) + '</span><b>' + fmt(r[1]) + '</b></div>';
      }).join('') + '</div>');
    }

    out.innerHTML = h.join('');
    out.hidden = false;
  }

  function wire() {
    var form = document.getElementById('askForm');
    if (!form || typeof HSAQuery === 'undefined') return;
    var input = document.getElementById('askInput');
    var chips = document.getElementById('askChips');
    var tag = document.getElementById('hsa-data');
    var P = tag ? JSON.parse(tag.textContent) : {};
    var ok = P.qa && HSAQuery.load(P.qa);
    var countEl = document.getElementById('askCount');
    if (!ok) {
      document.getElementById('askPanel').hidden = true;
      return;
    }
    if (countEl) countEl.textContent = fmt(P.qa.n) + ' human samples, '
      + HSAQuery.fields().length + ' fields, answered in your browser';
    chips.innerHTML = EXAMPLES.map(function (e) {
      return '<button type="button">' + esc(e) + '</button>';
    }).join('');
    chips.addEventListener('click', function (e) {
      if (e.target.tagName !== 'BUTTON') return;
      input.value = e.target.textContent;
      render(HSAQuery.ask(input.value));
    });
    form.addEventListener('submit', function (e) {
      e.preventDefault();
      if (input.value.trim()) render(HSAQuery.ask(input.value));
    });
  }
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', wire);
  else wire();
})();
