/* Renders every number, table and chart on the site from window.RESULTS (written by src/evaluate.py),
   and runs the in-browser split sandbox. No libraries. */
(function () {
  'use strict';
  var R = window.RESULTS;
  if (!R) { document.body.insertAdjacentHTML('afterbegin', '<p style="padding:1rem">results.js is missing. Run <code>python src/evaluate.py</code> first.</p>'); return; }
  var S = R.summary, C = R.config, UNIT = 'kcal/mol', REF = R.ref;
  var M = {}; R.models.forEach(function (m) { M[m.name] = m; });
  var SEL = M[R.selected], MEAN = M['Mean'], RM = M[REF], REP = R.repeats, NREP = C.n_repeats;
  var AAS = 'ACDEFGHIKLMNPQRSTVWY', CATS = AAS.split('').concat(['EDGE']);
  var LEARNED = R.models.filter(function (m) { return m.kind === 'model'; });

  // ------------------------------------------------------------------ helpers
  var DUMMY = document.createElement('div');
  function $(id) { return document.getElementById(id) || DUMMY; }
  function has(id) { return !!document.getElementById(id); }
  function esc(s) { return String(s).replace(/[&<>"]/g, function (c) { return { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]; }); }
  function f(x, d) { return (x == null || isNaN(x)) ? '–' : Number(x).toFixed(d == null ? 2 : d).replace('-', '−'); }
  function f3(x) { return f(x, 3); }
  function sg(x, d) { return (x == null || isNaN(x)) ? '–' : (x >= 0 ? '+' : '−') + Math.abs(x).toFixed(d == null ? 3 : d); }
  function pct(x, d) { return (100 * x).toFixed(d || 0) + '%'; }
  function n0(x) { return Number(x).toLocaleString('en-US'); }
  function get(path) { return path.split('.').reduce(function (o, k) { return o == null ? o : o[k]; }, R); }
  function mean(a) { var s = 0; for (var i = 0; i < a.length; i++) s += a[i]; return s / a.length; }
  function maeOf(a, b) { var s = 0; for (var i = 0; i < a.length; i++) s += Math.abs(a[i] - b[i]); return s / a.length; }
  function niceStep(span, n) { var raw = span / n, p = Math.pow(10, Math.floor(Math.log10(raw))), r = raw / p; return (r < 1.5 ? 1 : r < 3.5 ? 2 : r < 7.5 ? 5 : 10) * p; }
  function widthOf(el) { return Math.max(280, Math.floor(el.getBoundingClientRect().width || 640)); }
  function tipAttr(html) { return ' data-tip="' + esc(html) + '"'; }
  function call(ci) { return (ci[0] < 0 && ci[1] > 0) ? 'same' : (ci[1] <= 0 ? 'better' : 'worse'); }
  function civ(ci) { return sg(ci[0]) + ' to ' + sg(ci[1]); }
  function pills(m) { return (m.kind === 'baseline' ? '<span class="pill">baseline</span>' : '') + (m.name === R.selected ? '<span class="pill sel">selected</span>' : ''); }

  document.querySelectorAll('[data-v]').forEach(function (el) {
    var v = get(el.getAttribute('data-v')), fmt = el.getAttribute('data-f');
    if (v == null) { el.textContent = '–'; return; }
    if (fmt === 'pct0') el.textContent = pct(v);
    else if (typeof v === 'number') el.textContent = Number.isInteger(v) ? n0(v) : f(v, fmt ? +fmt : 2);
    else el.textContent = v;
  });

  // ------------------------------------------------------------------ headline
  var VM = call(SEL.diff_vs_mean_ci), VR = call(SEL.diff_vs_ref_ci);
  var isPreview = /PREVIEW/.test(C.dataset_label);
  var lead = '<strong>' + { better: 'Partly.', same: 'Not reliably.', worse: 'No.' }[VM] + '</strong> On ' + n0(S.n_test) + ' test mutations from ' + S.c_test + ' protein clusters the models never saw, the selected model (' + esc(SEL.name) + ') was off by ' + f3(SEL.test_mae) + ' ' + UNIT + ' on average (MAE), against ' + f3(MEAN.test_mae) + ' for always predicting the average. ' +
    'A plain table of the average effect of each amino-acid swap gets ' + f3(RM.test_mae) + ': ' +
    { same: 'the selected model cannot be told apart from it.', better: 'the selected model beats it by a small but measurable margin.', worse: 'the table does better than the selected model.' }[VR];
  $('answer').innerHTML = lead;
  $('facts').innerHTML =
    '<div><dt>Mean baseline, test MAE</dt><dd>' + f3(MEAN.test_mae) + '<small>' + UNIT + '</small></dd></div>' +
    '<div><dt>Substitution table, test MAE</dt><dd>' + f3(RM.test_mae) + '<small>' + UNIT + '</small></dd></div>' +
    '<div><dt>' + esc(SEL.name) + ', test MAE</dt><dd>' + f3(SEL.test_mae) + '<small>vs mean ' + sg(SEL.diff_vs_mean) + ' (95% interval ' + civ(SEL.diff_vs_mean_ci) + ')</small></dd></div>';
  $('byline').textContent = (C.author ? 'By ' + C.author + '. ' : '') + 'Results generated ' + R.generated + (isPreview ? ' from a 10,000-row preview sample, not the full data.' : '.');

  // ------------------------------------------------------------------ data
  var q = [];
  if (S.n_raw) q.push('<strong>Most rows are not single substitutions.</strong> Of ' + n0(S.n_raw) + ' rows in the file, ' + n0(S.n_single) + ' are single substitutions. Wild types, insertions, deletions and multiple substitutions are a different question and are left out.');
  if (S.n_no_target) q.push('<strong>' + pct(S.n_no_target / S.n_single) + ' of single substitutions have no ΔΔG.</strong> ' + n0(S.n_no_target) + ' rows hold <code>-</code> in <code>ddG_ML</code>, mostly because the variant or its wild type fell outside the range the assay can measure. Mutations that unfold the protein completely are the most likely to be missing, so very destabilizing mutations are under-represented.');
  q.push('<strong>The target is noisy.</strong> ΔΔG is the difference of two stabilities, each inferred from sequencing counts through a fitted model, not read off an instrument.');
  q.push('<strong>Most effects are small.</strong> ' + pct(S.share_small) + ' of kept mutations have |ΔΔG| under 0.5 ' + UNIT + ' and ' + pct(S.share_destabilizing) + ' are destabilizing. Mean ' + f(S.target_mean) + ', standard deviation ' + f(S.target_sd) + '.');
  q.push('<strong>Rows are not independent.</strong> The ' + n0(S.n_kept) + ' mutations come from ' + S.n_proteins + ' proteins in ' + S.n_clusters + ' clusters. The largest cluster (' + esc(S.largest_cluster) + ') alone holds ' + n0(S.largest_cluster_n) + ' rows, ' + pct(S.largest_cluster_n / S.n_kept) + ' of the data.');
  $('quality').innerHTML = q.map(function (t) { return '<li>' + t + '</li>'; }).join('');

  // ------------------------------------------------------------------ preprocessing
  var prep = [
    ['Keep rows whose <code>mut_type</code> matches <code>^[A-Z]\\d+[A-Z]$</code>', 'Single substitutions only: the question is about one swap at a time.', 'Fixed rule'],
    ['Convert <code>ddG_ML</code> to a number and drop rows where that fails', 'Those rows have no target to learn from or score against.', 'Fixed rule'],
    ['Split <code>mut_type</code> into original residue, position and new residue; add sequence length', 'These are the F1 features.', 'Nothing learned'],
    ['Check that every sequence uses the 20 amino-acid letters, the position is inside it, and the new residue sits at that position', 'Catches mislabelled rows. Every row checked passed.', 'Nothing learned'],
    ['Add the two residues on each side of the mutation, with <code>EDGE</code> past the ends of the chain', 'These four extra columns make F2. <code>EDGE</code> avoids missing values at the ends.', 'Nothing learned'],
    ['One-hot encode the amino acids and standard-scale position and length', 'Models need numbers. Scaling makes ridge\'s penalty fair across columns.', 'Training rows only: both sit inside the model pipeline'],
    ['No imputation, no clipping', 'No feature is missing, and ΔΔG is used as provided.', '–']
  ];
  $('prepRows').innerHTML = prep.map(function (r) { return '<tr><td>' + r[0] + '</td><td>' + r[1] + '</td><td>' + r[2] + '</td></tr>'; }).join('');

  // ------------------------------------------------------------------ methods
  var INFO = {
    'Mean': ['Predict the training average for every mutation.', 'The floor. Anything useful must beat it.'],
    'Mutant residue mean': ['Predict the training average for the amino acid being introduced.', 'Checks how much is just "proline is bad".'],
    'Substitution table': ['Predict the training average for each original → new pair (380 pairs).', 'A strong, cheap baseline. A model must add something beyond which swap it is.'],
    'Ridge': ['Linear model with a penalty on the size of the weights.', 'The simplest learned model. Additive effects only.'],
    'Decision tree': ['One tree of yes/no questions; predicts the average in each leaf.', 'The simplest non-linear model.'],
    'Random forest': ['Average of 100 trees, each on a bootstrap sample of rows.', 'Captures interactions; averaging lowers variance.'],
    'Extra trees': ['Like a random forest, but split points are chosen at random.', 'More randomness, often less overfitting.'],
    'AdaBoost': ['Shallow trees in sequence, each focusing on the rows the earlier ones got wrong.', 'The classic boosting method.'],
    'HistGradientBoosting': ['Trees added one at a time, each fitted to the remaining error, on binned features.', 'Usually the strongest method on tables like this.'],
    'XGBoost': ['Gradient boosting with regularization and row subsampling.', 'A widely used strong learner, for comparison.']
  };
  var algs = []; R.models.forEach(function (m) { if (algs.indexOf(m.algorithm) < 0) algs.push(m.algorithm); });
  $('methodRows').innerHTML = algs.map(function (a) {
    var i = INFO[a] || ['', ''], base = M[a];
    return '<tr><td><strong>' + esc(a) + '</strong>' + (base && base.kind === 'baseline' ? '<span class="pill">baseline</span>' : '') + '</td><td>' + i[0] + '</td><td>' + i[1] + '</td><td>' + (base ? '–' : 'F1 and F2') + '</td></tr>';
  }).join('');
  $('featRows').innerHTML = '<tr><td>F1</td><td>' + S.features.F1.map(function (x) { return '<code>' + esc(x) + '</code>'; }).join(', ') + '</td><td>Which swap, where, and how long the protein is</td></tr>' +
    '<tr><td>F2</td><td>F1 + ' + S.features.F2.filter(function (x) { return S.features.F1.indexOf(x) < 0; }).map(function (x) { return '<code>' + esc(x) + '</code>'; }).join(', ') + '</td><td>Adds the neighbouring residues</td></tr>';

  // F1 vs F2
  var pairs = algs.filter(function (a) { return M[a + ' (F1)'] && M[a + ' (F2)']; });
  var f2wins = pairs.filter(function (a) { return M[a + ' (F2)'].val_mae < M[a + ' (F1)'].val_mae; });
  var maxGap = Math.max.apply(null, pairs.map(function (a) { return Math.abs(M[a + ' (F2)'].val_mae - M[a + ' (F1)'].val_mae); }));
  $('f1f2Note').innerHTML = 'F2 had the lower validation MAE for ' + f2wins.length + ' of ' + pairs.length + ' algorithms. The largest gap is ' + f3(maxGap) + ' ' + UNIT + '. For scale, the selected model\'s test MAE has a 95% interval ' + f3(SEL.test_mae_ci[1] - SEL.test_mae_ci[0]) + ' wide, so these differences are small next to the uncertainty of a single split.';

  // ------------------------------------------------------------------ leakage
  var lk = R.leak.slice().sort(function (a, b) { return (b.grouped - b.random) / b.grouped - (a.grouped - a.random) / a.grouped; })[0];
  var lkGain = (lk.grouped - lk.random) / lk.grouped;
  $('leakTitle').textContent = 'Random folds make ' + lk.name + ' look ' + pct(lkGain) + ' better than it is';
  $('leakNote').innerHTML = 'Same rows, same models, ' + C.leak_folds + '-fold cross-validation on the training and validation rows. Random folds put mutations of the same protein on both sides, so a flexible model can memorise each protein. The baselines barely move, because a table of averages has little to memorise. The question is about unseen families, so the grouped number is the one that answers it.';

  // ------------------------------------------------------------------ results
  var rs = REP[SEL.name];
  var vtxt = esc(SEL.name) + ' had the lowest validation MAE of the ' + LEARNED.length + ' learned models, so it was selected before the test set was scored. On the test set its MAE is ' + f3(SEL.test_mae) + ' ' + UNIT + ' (95% interval ' + f3(SEL.test_mae_ci[0]) + ' to ' + f3(SEL.test_mae_ci[1]) + '), R² ' + f(SEL.test_r2) + '. ';
  vtxt += 'Against the mean baseline the difference is ' + sg(SEL.diff_vs_mean) + ' (' + civ(SEL.diff_vs_mean_ci) + '), ' + { better: 'a real improvement', same: 'within the uncertainty', worse: 'worse' }[VM] + '. ';
  vtxt += 'Against the substitution table it is ' + sg(SEL.diff_vs_ref) + ' (' + civ(SEL.diff_vs_ref_ci) + '), ' + { better: 'a measurable improvement', same: 'within the uncertainty', worse: 'so the table is better' }[VR] + '.';
  if (rs) vtxt += ' Across ' + NREP + ' repeated grouped splits of the training and validation rows, it beat the mean in ' + rs.n_better_mean + ' and the table in ' + rs.n_better_ref + '.';
  ['verdict', 'briefText', 'tourVerdict'].forEach(function (id) { $(id).innerHTML = vtxt; });
  $('resultRows').innerHTML = R.models.map(function (m) {
    return '<tr' + (m.name === R.selected ? ' class="pick"' : '') + '><td>' + esc(m.name) + pills(m) + '</td><td class="r">' + f3(m.val_mae) + '</td><td class="r">' + f3(m.test_mae) + '</td><td class="r">' + f3(m.test_mae_ci[0]) + ' to ' + f3(m.test_mae_ci[1]) + '</td><td class="r">' + f3(m.test_rmse) + '</td><td class="r">' + f(m.test_r2) + '</td><td class="r">' + (m.test_spearman == null ? 'n/a' : f(m.test_spearman)) + '</td><td class="r">' + (m.test_within_spearman == null ? 'n/a' : f(m.test_within_spearman)) + '</td></tr>';
  }).join('');
  $('resultCaption').textContent = 'Validation: ' + n0(S.n_val) + ' mutations in ' + S.c_val + ' clusters. Test: ' + n0(S.n_test) + ' mutations in ' + S.c_test + ' clusters. MAE and RMSE in ' + UNIT + '. Intervals come from resampling whole test clusters (' + C.n_bootstrap + ' draws). Spearman ranks all test mutations together; within-protein Spearman is the average inside each test protein with at least 20 mutations, which is the ranking a designer choosing among mutants of one protein would use. The mean predicts a constant, so it has no ranking.';
  var anyRef = LEARNED.filter(function (m) { return m.diff_vs_ref_ci[1] < 0; });
  $('forestTitle').textContent = LEARNED.some(function (m) { return m.diff_vs_mean_ci[1] < 0; }) ? 'Most models beat predicting the mean' : 'No model is reliably better than predicting the mean';
  $('forestTitle2').textContent = anyRef.length ? 'Beating the substitution table: ' + anyRef.map(function (m) { return m.name; }).join(', ') : 'No learned model reliably beats the substitution table';
  $('contextNote').innerHTML = 'For scale, ΔΔG has a standard deviation of ' + f(S.target_sd) + ' ' + UNIT + '. The selected model\'s predictions have a standard deviation of ' + f(R.residuals.pred_sd) + ' on the test set, against ' + f(R.residuals.y_sd) + ' for the measured values: it hedges toward the average and rarely predicts large effects.';
  var imp = R.importance;
  $('importanceNote').innerHTML = 'Shuffling <code>' + esc(imp[0].feature) + '</code> raises the selected model\'s validation MAE by ' + f3(imp[0].mae_increase) + ' ' + UNIT + ', more than any other input. ' + (imp.length > 1 ? 'Next is <code>' + esc(imp[1].feature) + '</code> (' + sg(imp[1].mae_increase) + '). ' : '') + 'An input whose shuffling changes nothing is not being used. This shows what the model relies on, not what causes stability.';

  // ------------------------------------------------------------------ limitations and conclusion
  var lim = [
    '<strong>The features are crude.</strong> They know the swap, its position and its neighbours, but not the 3D structure, so they cannot tell whether a residue is buried. The result is about these features, not about what sequence could predict.',
    '<strong>Only measurable mutations are in the data.</strong> ' + (S.n_no_target ? pct(S.n_no_target / S.n_single) + ' of single substitutions have no ΔΔG. ' : '') + 'The models are scored on in-range mutations; their errors on strongly destabilizing ones are unknown.',
    '<strong>The test set is smaller than it looks.</strong> ' + n0(S.n_test) + ' mutations but only ' + S.c_test + ' independent clusters, so the intervals are wide and one large cluster can move the result.',
    '<strong>The split proportions are uneven.</strong> Grouping whole clusters gave ' + pct(S.n_train / S.n_kept) + ' / ' + pct(S.n_val / S.n_kept) + ' / ' + pct(S.n_test / S.n_kept) + ' instead of 70 / 15 / 15, because cluster sizes are very unequal.',
    '<strong>No hyperparameter tuning.</strong> Each algorithm used one fixed setting. Selection was only between algorithms and feature sets, on validation MAE.',
    '<strong>The test set had been scored once before</strong>, for an earlier selected model, before this evaluation script existed. Selection here uses validation only, and the repeated splits do not use the test set.',
    '<strong>Column meanings are inferred.</strong> The file has no data dictionary; the sign of ΔΔG and the designed/natural reading of cluster labels come from the data and the paper.',
    '<strong>Prediction is not causation.</strong> An association between a residue and stability does not show that introducing it causes the change.'
  ];
  if (isPreview) lim.unshift('<strong>These numbers are a preview.</strong> They come from a 10,000-row sample. Run <code>src/evaluate.py</code> on the full splits to replace every number and chart.');
  $('limits').innerHTML = lim.map(function (t) { return '<li>' + t + '</li>'; }).join('');
  $('conclusionText').innerHTML = 'Using only which amino acid is swapped, where, and its neighbours, ' + esc(SEL.name) + ' predicts ΔΔG for unseen protein families with MAE ' + f3(SEL.test_mae) + ' ' + UNIT + ', ' + sg(SEL.diff_vs_mean) + ' relative to the mean (' + civ(SEL.diff_vs_mean_ci) + '). ' +
    (VR === 'same' ? 'It does not measurably beat a table of average effects per swap, so most of what it learned is which swap is being made.' : VR === 'better' ? 'It beats a table of average effects per swap by ' + sg(-SEL.diff_vs_ref) + ', so position or context adds a little.' : 'A table of average effects per swap did better.') +
    ' The neighbouring-residue features (F2) made no difference larger than the noise.';
  $('repoLine').innerHTML = 'All code, the frozen splits\' recipe and the full results are in the ' + (C.repo_url ? '<a href="' + esc(C.repo_url) + '">source repository</a>.' : 'source repository.');
  var ver = R.versions;
  $('envLine').textContent = 'Python ' + ver.python + ', scikit-learn ' + ver['scikit-learn'] + ', pandas ' + ver.pandas + ', numpy ' + ver.numpy + ', scipy ' + ver.scipy + '.';
  $('foot').textContent = 'Built from src/evaluate.py output generated ' + R.generated + '. Data: Tsuboyama et al. (2023), CC BY 4.0.';

  // ------------------------------------------------------------------ charts
  function chartHist(el) {
    var H = S.hist, w = widthOf(el), h = 230, m = { l: 48, r: 8, t: 8, b: 38 }, bw = H.step, nb = H.designed.length, x0 = H.lo, x1 = x0 + nb * bw, a = H.designed, b = H.natural;
    var ymax = Math.max.apply(null, a.map(function (v, i) { return v + b[i]; })), ys = niceStep(ymax, 4); ymax = Math.ceil(ymax / ys) * ys;
    var X = function (v) { return m.l + (v - x0) / (x1 - x0) * (w - m.l - m.r); }, Y = function (v) { return h - m.b - v / ymax * (h - m.t - m.b); };
    var s = '<svg viewBox="0 0 ' + w + ' ' + h + '" role="img" aria-label="Histogram of ΔΔG">';
    for (var t = 0; t <= ymax; t += ys) s += '<line class="grid" x1="' + m.l + '" x2="' + (w - m.r) + '" y1="' + Y(t) + '" y2="' + Y(t) + '"/><text class="muted" x="' + (m.l - 6) + '" y="' + (Y(t) + 4) + '" text-anchor="end">' + n0(t) + '</text>';
    var pw = Math.max(1, X(x0 + bw) - X(x0) - 1);
    for (var i = 0; i < nb; i++) {
      if (!(a[i] + b[i])) continue;
      var lo = x0 + i * bw;
      s += '<g' + tipAttr('<b>' + f(lo) + ' to ' + f(lo + bw) + ' ' + UNIT + '</b><br>' + n0(a[i]) + ' designed, ' + n0(b[i]) + ' natural') + '><rect class="f1" x="' + (X(lo) + .5) + '" y="' + Y(a[i]) + '" width="' + pw + '" height="' + (Y(0) - Y(a[i])) + '"/>' +
        (b[i] ? '<rect class="f2" x="' + (X(lo) + .5) + '" y="' + Y(a[i] + b[i]) + '" width="' + pw + '" height="' + Math.max(0, Y(a[i]) - Y(a[i] + b[i]) - 1) + '"/>' : '') + '</g>';
    }
    var xs = (x1 - x0) > 10 ? 2 : 1;
    for (var v = Math.ceil(x0 / xs) * xs; v <= x1; v += xs) s += '<text class="muted" x="' + X(v) + '" y="' + (h - m.b + 16) + '" text-anchor="middle">' + f(v, 0) + '</text>';
    s += '<line class="axis" x1="' + m.l + '" x2="' + (w - m.r) + '" y1="' + Y(0) + '" y2="' + Y(0) + '"/><text x="' + (m.l + (w - m.l - m.r) / 2) + '" y="' + (h - 4) + '" text-anchor="middle">Measured ΔΔG (' + UNIT + '), positive = more stable</text></svg>';
    el.innerHTML = s;
  }
  function chartBars(el, rows, xmax, axisLabel) {
    var w = widthOf(el), narrow = w < 460, lw = narrow ? 0 : 170, bh = 11, bgap = 3;
    var nb = rows[0].bars.length, rh = (narrow ? 18 : 0) + nb * (bh + bgap) + 12, m = { t: 4, b: axisLabel ? 34 : 6, r: 64 };
    var h = m.t + rows.length * rh + m.b, X = function (v) { return lw + v / xmax * (w - lw - m.r); };
    var s = '<svg viewBox="0 0 ' + w + ' ' + h + '" role="img">', step = niceStep(xmax, narrow ? 3 : 5);
    if (axisLabel) for (var t = 0; t <= xmax + 1e-9; t += step) s += '<line class="grid" x1="' + X(t) + '" x2="' + X(t) + '" y1="' + m.t + '" y2="' + (h - m.b) + '"/><text class="muted" x="' + X(t) + '" y="' + (h - m.b + 15) + '" text-anchor="middle">' + (+t.toFixed(2)) + '</text>';
    rows.forEach(function (r, i) {
      var y = m.t + i * rh;
      s += narrow ? '<text x="0" y="' + (y + 13) + '">' + esc(r.label) + '</text>' : '<text x="' + (lw - 10) + '" y="' + (y + 6 + nb * (bh + bgap) / 2) + '" text-anchor="end">' + esc(r.label) + '</text>';
      r.bars.forEach(function (b, j) {
        var by = y + (narrow ? 18 : 2) + j * (bh + bgap);
        s += '<g' + (b.tip ? tipAttr(b.tip) : '') + '><rect class="' + b.cls + '" x="' + X(0) + '" y="' + by + '" width="' + Math.max(1, X(b.v) - X(0)) + '" height="' + bh + '" rx="2"/><text class="val" x="' + (X(b.v) + 6) + '" y="' + (by + bh - 1.5) + '">' + esc(b.text) + '</text></g>';
      });
    });
    s += '<line class="axis" x1="' + X(0) + '" x2="' + X(0) + '" y1="' + m.t + '" y2="' + (h - m.b) + '"/>';
    if (axisLabel) s += '<text x="' + (lw + (w - lw - m.r) / 2) + '" y="' + (h - 3) + '" text-anchor="middle">' + esc(axisLabel) + '</text>';
    el.innerHTML = s + '</svg>';
  }
  function chartLeak(el) {
    var xmax = Math.max.apply(null, R.leak.map(function (r) { return Math.max(r.grouped, r.random); })) * 1.05;
    chartBars(el, R.leak.map(function (r) { return { label: r.name, bars: [
      { v: r.grouped, cls: 'f1', text: f3(r.grouped), tip: '<b>' + esc(r.name) + '</b><br>Folds keep clusters together: ' + f3(r.grouped) },
      { v: r.random, cls: 'f2', text: f3(r.random), tip: '<b>' + esc(r.name) + '</b><br>Random folds: ' + f3(r.random) }] }; }), xmax, 'Cross-validated MAE (' + UNIT + ')');
  }
  function chartF1F2(el) {
    var xmax = Math.max.apply(null, pairs.map(function (a) { return Math.max(M[a + ' (F1)'].val_mae, M[a + ' (F2)'].val_mae); })) * 1.05;
    chartBars(el, pairs.map(function (a) { var x = M[a + ' (F1)'], y = M[a + ' (F2)']; return { label: a, bars: [
      { v: x.val_mae, cls: 'f1', text: f3(x.val_mae), tip: '<b>' + esc(a) + ', F1</b><br>Validation MAE ' + f3(x.val_mae) },
      { v: y.val_mae, cls: 'f2', text: f3(y.val_mae), tip: '<b>' + esc(a) + ', F2</b><br>Validation MAE ' + f3(y.val_mae) }] }; }), xmax, 'Validation MAE (' + UNIT + ')');
  }
  function chartImportance(el) {
    var xmax = Math.max.apply(null, R.importance.map(function (d) { return d.mae_increase + d.sd; })) * 1.1 || 0.01;
    chartBars(el, R.importance.map(function (d) { return { label: d.feature, bars: [{ v: Math.max(0, d.mae_increase), cls: 'f1', text: sg(d.mae_increase), tip: '<b>' + esc(d.feature) + '</b><br>MAE increase ' + sg(d.mae_increase, 4) + ' ± ' + f(d.sd, 4) }] }; }), xmax, 'Increase in validation MAE when shuffled (' + UNIT + ')');
  }
  function chartForest(el, base) {
    var isMean = base === 'Mean', dk = isMean ? 'diff_vs_mean' : 'diff_vs_ref', rk = isMean ? 'Mean' : REF;
    var w = widthOf(el), narrow = w < 460, lw = narrow ? 0 : 190, rh = narrow ? 52 : 34, m = { t: 22, b: 40, r: 14 };
    var rows = R.models.filter(function (x) { return x.name !== base; }), all = [0];
    rows.forEach(function (r) { all.push(r[dk + '_ci'][0], r[dk + '_ci'][1]); if (REP[r.name]) REP[r.name].mae.forEach(function (v, i) { all.push(v - REP[rk].mae[i]); }); });
    var lo = Math.min.apply(null, all), hi = Math.max.apply(null, all), pad = (hi - lo) * 0.06; lo -= pad; hi += pad;
    var h = m.t + rows.length * rh + m.b, X = function (v) { return lw + 8 + (v - lo) / (hi - lo) * (w - lw - 8 - m.r); };
    var s = '<svg viewBox="0 0 ' + w + ' ' + h + '" role="img" aria-label="Difference in test MAE">', step = niceStep(hi - lo, narrow ? 4 : 6);
    for (var t = Math.ceil(lo / step) * step; t <= hi; t += step) { if (Math.abs(t) < 1e-9) continue; s += '<line class="grid" x1="' + X(t) + '" x2="' + X(t) + '" y1="' + m.t + '" y2="' + (h - m.b) + '"/><text class="muted" x="' + X(t) + '" y="' + (h - m.b + 15) + '" text-anchor="middle">' + sg(t, 2) + '</text>'; }
    s += '<line class="zero" x1="' + X(0) + '" x2="' + X(0) + '" y1="' + (m.t - 4) + '" y2="' + (h - m.b) + '"/><text x="' + X(0) + '" y="' + (h - m.b + 15) + '" text-anchor="middle">0</text>';
    s += '<text class="muted" x="' + (X(0) - 8) + '" y="12" text-anchor="end">← better</text><text class="muted" x="' + (X(0) + 8) + '" y="12">worse →</text>';
    rows.forEach(function (r, i) {
      var y = m.t + i * rh, cy = y + (narrow ? 28 : 12), c = r[dk + '_ci'];
      s += narrow ? '<text x="0" y="' + (y + 12) + '">' + esc(r.name) + '</text>' : '<text x="' + (lw - 10) + '" y="' + (cy + 4) + '" text-anchor="end">' + esc(r.name) + '</text>';
      if (REP[r.name]) REP[r.name].mae.forEach(function (v, k) { var d = v - REP[rk].mae[k]; s += '<circle class="rep" cx="' + X(d) + '" cy="' + (cy + 10) + '" r="2.5"' + tipAttr('<b>' + esc(r.name) + '</b><br>Repeat ' + (k + 1) + ': ' + sg(d) + ' ' + UNIT) + '/>'; });
      s += '<g' + tipAttr('<b>' + esc(r.name) + '</b><br>Test difference ' + sg(r[dk]) + ' ' + UNIT + '<br>95% interval ' + civ(c)) + '><line class="ci" x1="' + X(c[0]) + '" x2="' + X(c[1]) + '" y1="' + cy + '" y2="' + cy + '"/><circle class="est" cx="' + X(r[dk]) + '" cy="' + cy + '" r="4.5"/></g>';
    });
    s += '<text x="' + (lw + (w - lw - m.r) / 2) + '" y="' + (h - 4) + '" text-anchor="middle">Test MAE minus the ' + esc(isMean ? 'mean baseline' : 'substitution table') + '\'s (' + UNIT + ')</text></svg>';
    el.innerHTML = s;
  }

  // ------------------------------------------------------------------ error explorer
  var exModel = $('exModel'), exStratum = $('exStratum');
  var shownModels = Object.keys(R.test_rows[0].pred);
  exModel.innerHTML = shownModels.map(function (n) { return '<option' + (n === R.selected ? ' selected' : '') + '>' + esc(n) + '</option>'; }).join('');
  var strataNames = Array.from(new Set(R.test_rows.map(function (r) { return r.stratum; }))).sort();
  exStratum.innerHTML = '<option value="">All substitution types</option>' + strataNames.map(function (s) { return '<option>' + esc(s) + '</option>'; }).join('');
  function exRows() { var st = exStratum.value; return R.test_rows.filter(function (r) { return !st || r.stratum === st; }); }
  function chartScatter(el) {
    var name = exModel.value, rows = exRows(), w = widthOf(el), m = { l: 44, r: 10, t: 8, b: 40 }, size = Math.min(w - m.l - m.r, 420), h = size + m.t + m.b;
    var vals = []; R.test_rows.forEach(function (r) { vals.push(r.y, r.pred[name]); });
    var lo = Math.floor(Math.min.apply(null, vals)), hi = Math.ceil(Math.max.apply(null, vals));
    var X = function (v) { return m.l + (v - lo) / (hi - lo) * size; }, Y = function (v) { return m.t + size - (v - lo) / (hi - lo) * size; };
    var s = '<svg viewBox="0 0 ' + w + ' ' + h + '" role="img" aria-label="Predicted against measured ΔΔG">', step = (hi - lo) > 8 ? 2 : 1;
    for (var t = lo; t <= hi; t += step) s += '<line class="grid" x1="' + X(t) + '" x2="' + X(t) + '" y1="' + m.t + '" y2="' + (m.t + size) + '"/><line class="grid" x1="' + m.l + '" x2="' + (m.l + size) + '" y1="' + Y(t) + '" y2="' + Y(t) + '"/><text class="muted" x="' + X(t) + '" y="' + (m.t + size + 16) + '" text-anchor="middle">' + f(t, 0) + '</text><text class="muted" x="' + (m.l - 6) + '" y="' + (Y(t) + 4) + '" text-anchor="end">' + f(t, 0) + '</text>';
    s += '<line class="diag" x1="' + X(lo) + '" y1="' + Y(lo) + '" x2="' + X(hi) + '" y2="' + Y(hi) + '"/>';
    rows.forEach(function (r) {
      var p = r.pred[name];
      s += '<circle class="pt" cx="' + X(r.y) + '" cy="' + Y(p) + '" r="3.5" style="fill:var(--' + (r.designed ? 's1' : 's2') + ');fill-opacity:.7"' + tipAttr('<b>' + esc(r.mut) + ' in ' + esc(r.wt) + '</b><br>' + esc(r.stratum) + ', cluster ' + esc(r.group) + '<br>Measured ' + f(r.y) + '<br>Predicted ' + f(p) + '<br>Error ' + sg(p - r.y, 2)) + '/>';
    });
    el.innerHTML = s + '<text x="' + (m.l + size / 2) + '" y="' + (h - 4) + '" text-anchor="middle">Measured ΔΔG (' + UNIT + ')</text><text transform="translate(12 ' + (m.t + size / 2) + ') rotate(-90)" text-anchor="middle">Predicted (' + UNIT + ')</text></svg>';
  }
  function renderExplorer() {
    if (!has('chartScatter')) return;
    var name = exModel.value, rows = exRows();
    chartScatter($('chartScatter'));
    var y = rows.map(function (r) { return r.y; });
    $('exReadout').innerHTML = rows.length ? '<strong>' + n0(rows.length) + ' test mutations shown</strong> (a random ' + n0(R.test_rows.length) + ' of ' + n0(S.n_test) + '). MAE of ' + esc(name) + ' here: ' + f3(maeOf(y, rows.map(function (r) { return r.pred[name]; }))) + '. Mean baseline on the same rows: ' + f3(maeOf(y, rows.map(function (r) { return r.pred['Mean']; }))) + '.' : 'No test mutations in this group.';
    var worst = rows.slice().sort(function (a, b) { return Math.abs(b.pred[name] - b.y) - Math.abs(a.pred[name] - a.y); }).slice(0, 8);
    $('worstRows').innerHTML = worst.map(function (r) { return '<tr><td>' + esc(r.mut) + ' in ' + esc(r.wt) + '</td><td>' + esc(r.stratum) + '</td><td>' + esc(r.group) + '</td><td class="r">' + f(r.y) + '</td><td class="r">' + f(r.pred[name]) + '</td><td class="r">' + sg(r.pred[name] - r.y, 2) + '</td></tr>'; }).join('');
  }
  exModel.addEventListener('change', renderExplorer); exStratum.addEventListener('change', renderExplorer);
  function sliceRows(list) { return list.map(function (s) { return '<tr><td>' + esc(s.slice) + '</td><td class="r">' + n0(s.n) + '</td><td class="r">' + f3(s.mae) + '</td><td class="r">' + sg(s.bias) + '</td><td class="r">' + f3(s.mae_ref) + '</td><td class="r">' + f3(s.mae_mean) + '</td></tr>'; }).join(''); }
  var RS = R.residuals;
  $('sliceRows').innerHTML = '<tr class="sub"><td colspan="6">By the amino acid introduced</td></tr>' + sliceRows(RS.by_type) + '<tr class="sub"><td colspan="6">By origin of the protein</td></tr>' + sliceRows(RS.by_origin) + '<tr class="sub"><td colspan="6">By size of the real effect</td></tr>' + sliceRows(RS.by_size);
  $('clusterRows').innerHTML = RS.by_cluster.map(function (s) { return '<tr><td>' + esc(s.slice) + '</td><td>' + s.origin + '</td><td class="r">' + s.proteins + '</td><td class="r">' + n0(s.n) + '</td><td class="r">' + f3(s.mae) + '</td><td class="r">' + sg(s.bias) + '</td><td class="r">' + f3(s.mae_mean) + '</td></tr>'; }).join('');
  var bigSize = RS.by_size.filter(function (s) { return s.slice === 'over 2'; })[0], smallSize = RS.by_size.filter(function (s) { return /under/.test(s.slice); })[0];
  var worstType = RS.by_type.slice().sort(function (a, b) { return b.mae - a.mae; })[0];
  var bigC = RS.by_cluster[0];
  $('errorNote').innerHTML = (bigSize && smallSize ? 'Errors grow with the size of the real effect: MAE ' + f3(smallSize.mae) + ' for mutations with |ΔΔG| under 0.5, and ' + f3(bigSize.mae) + ' for those over 2, where the model predicts too ' + (bigSize.bias > 0 ? 'high' : 'low') + ' by ' + f(Math.abs(bigSize.bias)) + ' on average. ' : '') +
    'The hardest substitution type is mutations ' + esc(worstType.slice) + ' (MAE ' + f3(worstType.mae) + '). The largest test cluster, ' + esc(bigC.slice) + ', holds ' + pct(bigC.n / S.n_test) + ' of the test set, so it weighs heavily in every test score.';

  // ------------------------------------------------------------------ sandbox (F1 features + ridge, in the browser)
  function rng(seed) { var a = (seed >>> 0) + 0x9e3779b9; return function () { a |= 0; a = a + 0x6D2B79F5 | 0; var t = Math.imul(a ^ a >>> 15, 1 | a); t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t; return ((t ^ t >>> 14) >>> 0) / 4294967296; }; }
  function shuffle(arr, r) { for (var i = arr.length - 1; i > 0; i--) { var j = Math.floor(r() * (i + 1)), t = arr[i]; arr[i] = arr[j]; arr[j] = t; } return arr; }
  function features(seq, mut) {   // F1, same column order as design_f1 in src/evaluate.py
    var mm = /^([A-Z])(\d+)([A-Z])$/.exec(mut), x = [];
    CATS.forEach(function (a) { x.push(a === mm[1] ? 1 : 0); }); CATS.forEach(function (a) { x.push(a === mm[3] ? 1 : 0); });
    x.push(+mm[2], seq.length); return x;
  }
  var SCALE = CATS.map(function () { return false; }).concat(CATS.map(function () { return false; })).concat([true, true]);
  function solve(A, b) {
    var n = b.length, i, j, k;
    for (i = 0; i < n; i++) {
      var p = i; for (k = i + 1; k < n; k++) if (Math.abs(A[k][i]) > Math.abs(A[p][i])) p = k;
      var t = A[i]; A[i] = A[p]; A[p] = t; var tb = b[i]; b[i] = b[p]; b[p] = tb;
      for (k = i + 1; k < n; k++) { var fc = A[k][i] / A[i][i]; if (!fc) continue; for (j = i; j < n; j++) A[k][j] -= fc * A[i][j]; b[k] -= fc * b[i]; }
    }
    var x = new Array(n); for (i = n - 1; i >= 0; i--) { var s = b[i]; for (j = i + 1; j < n; j++) s -= A[i][j] * x[j]; x[i] = s / A[i][i]; }
    return x;
  }
  function ridge(X, y, tr, alpha) {  // centre every column, scale only position and length (as StandardScaler does in the pipeline)
    var p = X[0].length, n = tr.length, mu = new Array(p).fill(0), sd = new Array(p).fill(1), ym = 0, j, k;
    tr.forEach(function (r) { ym += y[r]; for (j = 0; j < p; j++) mu[j] += X[r][j]; }); ym /= n; for (j = 0; j < p; j++) mu[j] /= n;
    for (j = 0; j < p; j++) if (SCALE[j]) { var v = 0; tr.forEach(function (r) { v += (X[r][j] - mu[j]) * (X[r][j] - mu[j]); }); sd[j] = Math.sqrt(v / n) || 1; }
    var A = [], b = new Array(p).fill(0), z = new Array(p); for (j = 0; j < p; j++) { A.push(new Array(p).fill(0)); A[j][j] = alpha; }
    tr.forEach(function (r) { for (j = 0; j < p; j++) z[j] = (X[r][j] - mu[j]) / sd[j]; var yc = y[r] - ym; for (j = 0; j < p; j++) { if (!z[j]) continue; b[j] += z[j] * yc; for (k = j; k < p; k++) A[j][k] += z[j] * z[k]; } });
    for (j = 0; j < p; j++) for (k = 0; k < j; k++) A[j][k] = A[k][j];
    var w = solve(A, b);
    return function (r) { var s = ym; for (var jj = 0; jj < p; jj++) s += w[jj] * (X[r][jj] - mu[jj]) / sd[jj]; return s; };
  }
  function tableMean(D, tr, keyFn) {
    var s = {}, c = {}, ym = 0;
    tr.forEach(function (r) { var k = keyFn(r); ym += D.y[r]; s[k] = (s[k] || 0) + D.y[r]; c[k] = (c[k] || 0) + 1; }); ym /= tr.length;
    return function (r) { var k = keyFn(r); return c[k] ? s[k] / c[k] : ym; };
  }
  function foldsGrouped(groups, k, r) {
    var by = {}; groups.forEach(function (g, i) { (by[g] = by[g] || []).push(i); });
    var keys = shuffle(Object.keys(by), r).sort(function (a, b) { return by[b].length - by[a].length; }), size = new Array(k).fill(0), fold = new Array(groups.length);
    keys.forEach(function (g) { var m = 0; for (var i = 1; i < k; i++) if (size[i] < size[m]) m = i; size[m] += by[g].length; by[g].forEach(function (i) { fold[i] = m; }); });
    return fold;
  }
  function foldsRandom(n, k, r) { var idx = shuffle(Array.from({ length: n }, function (_, i) { return i; }), r), fold = new Array(n); idx.forEach(function (row, i) { fold[row] = i % k; }); return fold; }
  function cv(D, fold, k, alpha) {
    var out = { mean: [], table: [], ridge: [] };
    for (var fo = 0; fo < k; fo++) {
      var tr = [], te = []; fold.forEach(function (v, i) { (v === fo ? te : tr).push(i); });
      if (!tr.length || !te.length) continue;
      var yt = te.map(function (i) { return D.y[i]; }), ym = mean(tr.map(function (i) { return D.y[i]; }));
      out.mean.push(maeOf(yt, te.map(function () { return ym; })));
      out.table.push(maeOf(yt, te.map(tableMean(D, tr, function (r) { return D.pair[r]; }))));
      out.ridge.push(maeOf(yt, te.map(ridge(D.X, D.y, tr, alpha))));
    }
    return { mean: mean(out.mean), table: mean(out.table), ridge: mean(out.ridge) };
  }
  function makeData(rows, label) {
    return { X: rows.map(function (r) { return features(r.seq, r.mut); }), y: rows.map(function (r) { return r.y; }), g: rows.map(function (r) { return r.group; }),
      pair: rows.map(function (r) { return r.mut[0] + r.mut[r.mut.length - 1]; }), label: label, n: rows.length };
  }
  var SB = R.sandbox.rows, DEFAULT = SB.length ? makeData(SB, 'a random ' + n0(SB.length) + ' rows of the project data') : null, DATA = DEFAULT, K = 5;
  function renderSandbox() {
    if (!has('sbAlpha') || !DATA) return;
    var a = Math.pow(10, +$('sbAlpha').value), seed = Math.max(0, parseInt($('sbSeed').value, 10) || 0);
    $('sbAlphaVal').textContent = a >= 10 ? a.toFixed(0) : a.toFixed(a >= 1 ? 1 : 2);
    var g = cv(DATA, foldsGrouped(DATA.g, K, rng(seed)), K, a), r = cv(DATA, foldsRandom(DATA.n, K, rng(seed)), K, a);
    var xmax = Math.max(g.mean, g.table, g.ridge, r.mean, r.table, r.ridge) * 1.05;
    chartBars($('chartSandbox'), [['Mean', 'mean'], ['Substitution table', 'table'], ['Ridge (F1)', 'ridge']].map(function (x) {
      return { label: x[0], bars: [{ v: g[x[1]], cls: 'f1', text: f3(g[x[1]]) }, { v: r[x[1]], cls: 'f2', text: f3(r[x[1]]) }] };
    }), xmax, 'Cross-validated MAE');
    var dg = (g.mean - g.ridge) / g.mean, dr = (r.mean - r.ridge) / r.mean;
    $('sbReadout').innerHTML = '<strong>' + esc(DATA.label) + ':</strong> ' + n0(DATA.n) + ' mutations in ' + new Set(DATA.g).size + ' groups. With folds that keep groups together, ridge is ' + pct(Math.abs(dg), 1) + (dg >= 0 ? ' better' : ' worse') + ' than the mean; with random folds, ' + pct(Math.abs(dr), 1) + (dr >= 0 ? ' better' : ' worse') + '.';
  }
  $('sbAlpha').addEventListener('input', renderSandbox); $('sbSeed').addEventListener('change', renderSandbox);
  $('sbShuffle').addEventListener('click', function () { $('sbSeed').value = (parseInt($('sbSeed').value, 10) || 0) + 1; renderSandbox(); });
  function parseCSV(text) {
    var rows = [], row = [], cur = '', inQ = false, i = text.charCodeAt(0) === 0xFEFF ? 1 : 0, ch;
    for (; i < text.length; i++) {
      ch = text[i];
      if (inQ) { if (ch === '"') { if (text[i + 1] === '"') { cur += '"'; i++; } else inQ = false; } else cur += ch; }
      else if (ch === '"') inQ = true; else if (ch === ',') { row.push(cur); cur = ''; }
      else if (ch === '\n' || ch === '\r') { if (ch === '\r' && text[i + 1] === '\n') i++; row.push(cur); cur = ''; if (row.length > 1 || row[0] !== '') rows.push(row); row = []; }
      else cur += ch;
    }
    if (cur !== '' || row.length) { row.push(cur); rows.push(row); }
    return rows;
  }
  var loaded = null;
  function fail(m) { $('sbError').textContent = m; }
  function guess(h, res) { for (var k = 0; k < res.length; k++) { var i = h.findIndex(function (x) { return res[k].test(x); }); if (i >= 0) return i; } return 0; }
  $('sbFile').addEventListener('change', function (ev) {
    var file = ev.target.files[0]; fail(''); if (!file) return;
    if (file.size > 60 * 1024 * 1024) { fail('This file is ' + Math.round(file.size / 1048576) + ' MB, too large for the browser. Use src/evaluate.py instead.'); return; }
    var fr = new FileReader();
    fr.onload = function () {
      var rows = parseCSV(String(fr.result));
      if (rows.length < 3) { fail('No data rows found.'); return; }
      loaded = { header: rows[0], rows: rows.slice(1), name: file.name };
      var opts = loaded.header.map(function (h, i) { return '<option value="' + i + '">' + esc(h) + '</option>'; }).join('');
      ['sbSeqCol', 'sbMutCol', 'sbTargetCol', 'sbGroupCol'].forEach(function (id) { $(id).innerHTML = opts; });
      $('sbSeqCol').value = guess(loaded.header, [/^aa_seq$/i, /seq/i]); $('sbMutCol').value = guess(loaded.header, [/^mut_type$/i, /mut/i]);
      $('sbTargetCol').value = guess(loaded.header, [/^ddg_ml$/i, /ddg/i, /target/i]); $('sbGroupCol').value = guess(loaded.header, [/cluster/i, /wt_name/i, /group/i]);
      $('sbMap').hidden = false;
    };
    fr.readAsText(file);
  });
  $('sbRun').addEventListener('click', function () {
    if (!loaded) return; fail('');
    var si = +$('sbSeqCol').value, mi = +$('sbMutCol').value, ti = +$('sbTargetCol').value, gi = +$('sbGroupCol').value, rows = [], bad = 0;
    loaded.rows.forEach(function (r) {
      var seq = (r[si] || '').trim().toUpperCase(), mut = (r[mi] || '').trim().toUpperCase(), y = parseFloat(r[ti]), g = (r[gi] || '').trim(), mm = /^([A-Z])(\d+)([A-Z])$/.exec(mut);
      if (!/^[ACDEFGHIKLMNPQRSTVWY]+$/.test(seq) || !mm || !isFinite(y) || !g || AAS.indexOf(mm[1]) < 0 || AAS.indexOf(mm[3]) < 0) { bad++; return; }
      rows.push({ seq: seq, mut: mut, y: y, group: g });
    });
    if (rows.length > 6000) { shuffle(rows, rng(0)); rows = rows.slice(0, 6000); }
    if (rows.length < 50) { fail('Only ' + rows.length + ' usable rows. The mutation column needs single substitutions like W34D and the target column numbers.'); return; }
    if (new Set(rows.map(function (r) { return r.group; })).size < K) { fail('The group column needs at least ' + K + ' distinct values.'); return; }
    DATA = makeData(rows, loaded.name + (bad ? ' (' + bad + ' rows skipped)' : '')); renderSandbox();
  });
  $('sbReset').addEventListener('click', function () { DATA = DEFAULT; fail(''); renderSandbox(); });

  // ------------------------------------------------------------------ claims, repository lines, contribution
  function claimBox(yes, no) { return '<div><h4>This supports</h4><ul>' + yes.map(function (t) { return '<li>' + t + '</li>'; }).join('') + '</ul></div><div><h4>This does not show</h4><ul>' + no.map(function (t) { return '<li>' + t + '</li>'; }).join('') + '</ul></div>'; }
  var CLAIMS = {
    results: claimBox(['On ' + S.c_test + ' unseen clusters, ' + esc(SEL.name) + ' differed from the mean by ' + sg(SEL.diff_vs_mean) + ' ' + UNIT + ' (' + civ(SEL.diff_vs_mean_ci) + ') and from the substitution table by ' + sg(SEL.diff_vs_ref) + ' (' + civ(SEL.diff_vs_ref_ci) + ').', 'How the split is made changes the answer: random folds made ' + esc(lk.name) + ' look ' + pct(lkGain) + ' better.'],
      ['That sequence cannot predict stability better. Only simple features and no 3D structure were tried.', 'That the result holds for mutations outside the assay\'s range, or for other datasets.', 'That any feature causes stability.']),
    split: claimBox(['Random folds overstate accuracy when relatives of the test rows are in training.', 'Folds that keep clusters together give a fairer estimate for unseen families.'],
      ['Anything about the tree models: this panel runs only the baselines and ridge.', 'That grouped CV is precise: with few clusters it is still noisy.']),
    errors: claimBox(['Which kinds of mutation the selected model misses most, and in which direction.', 'Which inputs the model relies on.'],
      ['That small slices have reliable errors.', 'That all error is model error: part is measurement noise.']),
    compare: claimBox(['How the methods rank on validation MAE, which is how the model was chosen.', 'Whether adding neighbouring residues (F2) helps each algorithm.'],
      ['That a gap smaller than the test intervals is real.', 'That the ranking holds with tuned hyperparameters or other features.'])
  };
  document.querySelectorAll('[data-claims]').forEach(function (el) { el.innerHTML = CLAIMS[el.getAttribute('data-claims')] || ''; });
  document.querySelectorAll('.repo').forEach(function (el) {
    var base = (C.repo_url || '').replace(/\/$/, '');
    el.innerHTML = 'In the repository: ' + el.getAttribute('data-items').split(';').map(function (p) { var t = '<code>' + esc(p) + '</code>'; return base ? '<a href="' + esc(base) + '/blob/main/' + esc(p) + '">' + t + '</a>' : t; }).join(', ') + '.';
  });
  if (has('myContribution')) $('myContribution').innerHTML = C.contribution ? '<p>' + esc(C.contribution) + '</p>' : '<p class="note"><strong>Not written yet.</strong> Add <code>contribution</code> in <code>config.json</code> and rerun <code>src/evaluate.py</code>.</p>';

  // ------------------------------------------------------------------ compare page
  function renderCompare() {
    if (!has('cmpA')) return;
    var a = M[$('cmpA').value], b = M[$('cmpB').value];
    $('cmpHA').textContent = a.name; $('cmpHB').textContent = b.name;
    function row(l, x, y, d) { return '<tr><td>' + l + '</td><td class="r">' + x + '</td><td class="r">' + y + '</td><td class="r">' + (d == null ? '' : sg(d)) + '</td></tr>'; }
    $('cmpRows').innerHTML = row('Validation MAE', f3(a.val_mae), f3(b.val_mae), b.val_mae - a.val_mae) + row('Validation R²', f(a.val_r2), f(b.val_r2), b.val_r2 - a.val_r2) +
      row('Test MAE', f3(a.test_mae), f3(b.test_mae), b.test_mae - a.test_mae) + row('95% interval of test MAE', f3(a.test_mae_ci[0]) + ' to ' + f3(a.test_mae_ci[1]), f3(b.test_mae_ci[0]) + ' to ' + f3(b.test_mae_ci[1]), null) +
      row('Test RMSE', f3(a.test_rmse), f3(b.test_rmse), b.test_rmse - a.test_rmse) + row('Test R²', f(a.test_r2), f(b.test_r2), b.test_r2 - a.test_r2) +
      row('Within-protein Spearman', f(a.test_within_spearman), f(b.test_within_spearman), (a.test_within_spearman == null || b.test_within_spearman == null) ? null : b.test_within_spearman - a.test_within_spearman);
    var ra = REP[a.name], rb = REP[b.name], txt = '';
    if (ra && rb && a.name !== b.name) { var w = 0; ra.mae.forEach(function (v, i) { if (v < rb.mae[i]) w++; }); txt = ' Across ' + NREP + ' repeated splits, ' + esc(a.name) + ' had the lower MAE in ' + w + ' and ' + esc(b.name) + ' in ' + (NREP - w) + '.'; }
    else if (a.name !== b.name) txt = ' Repeated splits were run only for the baselines, ridge and the selected model.';
    var ov = !(a.test_mae_ci[1] < b.test_mae_ci[0] || b.test_mae_ci[1] < a.test_mae_ci[0]);
    $('cmpReadout').innerHTML = a.name === b.name ? 'Pick two different methods.' : 'Their test MAE intervals ' + (ov ? 'overlap, so this test cannot separate them.' : 'do not overlap.') + txt;
  }
  if (has('cmpA')) {
    var opts = R.models.map(function (m) { return '<option>' + esc(m.name) + '</option>'; }).join('');
    $('cmpA').innerHTML = opts; $('cmpB').innerHTML = opts; $('cmpA').value = REF; $('cmpB').value = R.selected;
    $('cmpA').addEventListener('change', renderCompare); $('cmpB').addEventListener('change', renderCompare);
    $('rankRows').innerHTML = R.models.slice().sort(function (a, b) { return a.val_mae - b.val_mae; }).map(function (m, i) {
      return '<tr' + (m.name === R.selected ? ' class="pick"' : '') + '><td class="r">' + (i + 1) + '</td><td>' + esc(m.name) + pills(m) + '</td><td class="r">' + f3(m.val_mae) + '</td><td class="r">' + f3(m.test_mae) + '</td><td class="r">' + (REP[m.name] ? f3(REP[m.name].mae_mean) : '–') + '</td></tr>';
    }).join('');
    renderCompare();
  }

  // ------------------------------------------------------------------ tour
  if (has('tourFill')) {
    var sampleRows = function () {
      var idx = [], n = SB.length; while (idx.length < Math.min(5, n)) { var i = Math.floor(Math.random() * n); if (idx.indexOf(i) < 0) idx.push(i); }
      $('tourRows').innerHTML = idx.map(function (i) { var r = SB[i], p = +r.mut.slice(1, -1) - 1, lo = Math.max(0, p - 6), hi = Math.min(r.seq.length, p + 7); return '<tr><td><code>' + esc(r.mut) + '</code></td><td><code>' + (lo ? '…' : '') + esc(r.seq.slice(lo, p)) + '<b>' + esc(r.seq[p]) + '</b>' + esc(r.seq.slice(p + 1, hi)) + (hi < r.seq.length ? '…' : '') + '</code></td><td class="r">' + f(r.y) + '</td><td>' + esc(r.group) + '</td></tr>'; }).join('');
    };
    sampleRows(); $('tourShuffle').addEventListener('click', sampleRows);
    $('tourBase').innerHTML = R.models.filter(function (m) { return m.kind === 'baseline'; }).concat([SEL]).map(function (m) { return '<tr><td>' + esc(m.name) + pills(m) + '</td><td class="r">' + f3(m.test_mae) + '</td><td class="r">' + f3(m.test_mae_ci[0]) + ' to ' + f3(m.test_mae_ci[1]) + '</td></tr>'; }).join('');
    var steps = Array.from(document.querySelectorAll('.tstep'));
    var showStep = function () {
      var mm = /^#step-(\d+)$/.exec(location.hash), n = mm ? Math.min(steps.length, Math.max(1, +mm[1])) - 1 : 0;
      steps.forEach(function (st, i) { st.classList.toggle('active', i === n); });
      $('tourFill').style.width = ((n + 1) / steps.length * 100) + '%'; $('tourCount').textContent = 'Step ' + (n + 1) + ' of ' + steps.length;
      var pv = $('tourPrev'), nx = $('tourNext');
      pv.href = '#step-' + Math.max(1, n); pv.style.visibility = n === 0 ? 'hidden' : 'visible';
      nx.href = n === steps.length - 1 ? 'results.html' : '#step-' + (n + 2); nx.textContent = n === steps.length - 1 ? 'See all results' : 'Next';
      renderAll(); window.scrollTo(0, 0);
    };
    window.addEventListener('hashchange', showStep);
  }

  // ------------------------------------------------------------------ tooltip, nav, render
  var tip = $('tip');
  function showTip(t, x, y) { tip.innerHTML = t.getAttribute('data-tip'); tip.style.display = 'block'; var r = tip.getBoundingClientRect(); tip.style.left = Math.max(8, Math.min(window.innerWidth - r.width - 8, x + 14)) + 'px'; tip.style.top = Math.max(8, (y + 16 + r.height > window.innerHeight ? y - r.height - 12 : y + 16)) + 'px'; }
  document.addEventListener('pointermove', function (e) { var t = e.target.closest && e.target.closest('[data-tip]'); if (t) showTip(t, e.clientX, e.clientY); else tip.style.display = 'none'; });
  document.addEventListener('pointerdown', function (e) { var t = e.target.closest && e.target.closest('[data-tip]'); if (t) showTip(t, e.clientX, e.clientY); else tip.style.display = 'none'; });
  window.addEventListener('scroll', function () { tip.style.display = 'none'; }, { passive: true });
  document.querySelectorAll('.topbar nav a').forEach(function (a) { var here = a.getAttribute('data-page') === document.body.getAttribute('data-page'); a.classList.toggle('here', here); if (here) a.setAttribute('aria-current', 'page'); });
  function renderAll() {
    var map = { chartHist: chartHist, chartLeak: chartLeak, chartF1F2: chartF1F2, chartImportance: chartImportance };
    Object.keys(map).forEach(function (id) { if (has(id)) map[id]($(id)); });
    if (has('chartForest')) chartForest($('chartForest'), 'Mean');
    if (has('chartForestRef')) chartForest($('chartForestRef'), REF);
    renderExplorer(); renderSandbox();
  }
  if (has('tourFill')) showStep(); else renderAll();
  var lastW = window.innerWidth, timer;
  window.addEventListener('resize', function () { if (window.innerWidth === lastW) return; lastW = window.innerWidth; clearTimeout(timer); timer = setTimeout(renderAll, 150); });
  window.__sandbox = { features: features, ridge: ridge, maeOf: maeOf, data: DEFAULT };   // used by tests/check_js.mjs
})();
