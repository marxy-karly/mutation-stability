"""Writes the seven HTML pages in docs/. Every number is filled in by docs/app.js from docs/results.js.
Run once after editing page text:  python tools/build_pages.py"""
from pathlib import Path

DOCS = Path(__file__).resolve().parent.parent / "docs"
NAV = [("index", "Home"), ("tour", "Guided tour"), ("split", "Split"), ("errors", "Errors"), ("compare", "Compare"), ("results", "Results"), ("reproduce", "Reproduce")]


def page(pid, title, body, desc="Predicting the stability change of single amino-acid mutations in unseen protein families."):
    nav = "".join(f'<a href="{p}.html" data-page="{p}">{t}</a>' for p, t in NAV)
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<meta name="description" content="{desc}">
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans:wght@400;500;600&family=Newsreader:opsz,wght@6..72,400;6..72,500;6..72,600&display=swap" rel="stylesheet">
<link rel="stylesheet" href="style.css">
</head>
<body data-page="{pid}">
<header class="topbar"><div class="topbar-in">
  <a class="brand" href="index.html">Mutation Stability</a>
  <nav aria-label="Pages">{nav}</nav>
</div></header>
<div class="shell">
<main>
{body}
</main>
<footer id="foot"></footer>
</div>
<div id="tip" role="status"></div>
<script src="results.js"></script>
<script src="about.js"></script>
<script src="app.js"></script>
</body>
</html>
"""


def head(kicker, h1, answer):
    return f'<header class="pagehead"><p class="kicker">{kicker}</p><h1>{h1}</h1><p class="answer">{answer}</p></header>'


LEAK_FIG = """<div class="figure">
    <h4 id="leakTitle"></h4>
    <p class="sub">Cross-validated MAE on training and validation rows (kcal/mol, lower is better), for the same methods under two ways of making folds.</p>
    <div class="legend"><span><i class="s1"></i>Folds keep clusters together</span><span><i class="s2"></i>Random folds</span></div>
    <div id="chartLeak"></div>
  </div>
  <p id="leakNote"></p>"""

HIST_FIG = """<div class="figure">
    <h4>Most mutations change stability only a little, and destabilizing ones outnumber stabilizing ones</h4>
    <p class="sub">Number of kept single substitutions by measured ΔΔG. Positive means more stable.</p>
    <div class="legend"><span><i class="s1"></i>Designed proteins</span><span><i class="s2"></i>Natural proteins</span></div>
    <div id="chartHist"></div>
  </div>"""

FS_FIG = """<div class="figure">
    <h4 id="fsTitle"></h4>
    <p class="sub">Validation MAE for each algorithm on each feature set it was run with (kcal/mol, lower is better).</p>
    <div class="legend" id="fsLegend"></div>
    <div id="chartFS"></div>
  </div>
  <p id="fsNote"></p>"""

FOREST = """<div class="figure">
    <h4 id="forestTitle"></h4>
    <p class="sub">Test MAE minus the mean baseline's test MAE. Left of zero is better.</p>
    <div class="legend"><span><i class="line"></i>Test set, with 95% interval</span><span><i class="dot"></i>One of <span data-v="config.n_repeats"></span> repeated splits (baselines, ridge, the selected model, and its base algorithm on each feature set)</span></div>
    <div id="chartForest"></div>
  </div>
  <div class="figure">
    <h4 id="forestTitle2"></h4>
    <p class="sub">The same comparison against the substitution table, which already knows which swap is being made.</p>
    <div id="chartForestRef"></div>
  </div>"""

INDEX = """<section class="hero">
  <p class="kicker">ES 335 Machine Learning · course project</p>
  <h1>Can we predict how much one amino-acid change shifts a protein's stability, for protein families the model has never seen?</h1>
  <p class="answer" id="answer"></p>
  <div class="cta"><a class="btn" href="tour.html">Start the guided tour</a><span class="cta-note">9 steps, about five minutes</span><a class="btn ghost" href="results.html">Go to the results</a></div>
  <ol class="flow" aria-label="How the project is built">
    <li><b>One mutation</b><span>A single swap, such as W34D, in a small protein</span></li>
    <li><b>Features</b><span>F1 the swap and position; F2 its neighbours; F3 sequence-derived numbers; F4 burial in an AlphaFold model</span></li>
    <li><b>Split by family</b><span>Whole clusters go to train, validation or test</span></li>
    <li><b>Models</b><span>Three baselines and the learned models, on up to four feature sets</span></li>
    <li><b>Select, then test</b><span>Chosen on validation MAE, scored once on test</span></li>
  </ol>
  <dl class="facts" id="facts"></dl>
  <p class="byline" id="byline"></p>
</section>
<section id="question">
  <h2><span class="n">01</span>Problem formulation</h2>
  <p>Protein engineers often ask whether a single amino-acid change will make a protein more or less stable. Testing that takes one experiment per mutation. A model that ranks candidate mutations first would save most of that work, but only if it works on proteins it was not trained on.</p>
  <dl class="frame">
    <dt>Question</dt><dd>Given a protein's sequence and one amino-acid substitution in it, how well can the change in folding stability (ΔΔG) be predicted when the protein comes from a family the model never saw?</dd>
    <dt>One row</dt><dd>One single-substitution variant of one small protein domain (<span data-v="summary.len_min"></span> to <span data-v="summary.len_max"></span> residues), for example W34D: the tryptophan at position 34 replaced by aspartate.</dd>
    <dt>Target</dt><dd><code>ddG_ML</code>, the variant's stability minus the unmutated protein's, in kcal/mol. Positive means the mutation stabilizes the protein. A regression task.</dd>
    <dt>When the prediction is made</dt><dd>Before any experiment, when only the sequence and the proposed swap exist.</dd>
    <dt>Available then</dt><dd>The swap, its position, the protein's sequence, and a predicted structure of the unmutated protein (AlphaFold), which needs no experiment. Every measurement column in the file is an output of the experiment and is not used.</dd>
    <dt>Limits of the claim</dt><dd>Single substitutions with a measurable ΔΔG, on this assay. Nothing about multiple mutations, insertions or deletions, and nothing about structures of the mutated protein.</dd>
  </dl>
  <p class="repo" data-items="src/preprocess.py"></p>
</section>
<section id="data">
  <h2><span class="n">02</span>Data</h2>
  <dl class="frame">
    <dt>Source</dt><dd><code>K50_dG_Dataset1_Dataset2.csv</code> in <code>Processed_K50_dG_datasets.zip</code>, from the data release of Tsuboyama et al. (2023), <a href="https://zenodo.org/records/7401275">Zenodo record 7401275</a>. The release accompanies the bioRxiv preprint. AlphaFold models of the unmutated proteins come from <code>AlphaFold_model_PDBs.zip</code> in the same record (13.7 MB).</dd>
    <dt>Size</dt><dd><span data-v="summary.n_raw"></span> rows and 38 columns in the file; <span data-v="summary.n_single"></span> single substitutions; <span data-v="summary.n_kept"></span> with a usable target, from <span data-v="summary.n_proteins"></span> proteins in <span data-v="summary.n_clusters"></span> clusters.</dd>
    <dt>Access</dt><dd>Open, CC BY 4.0: reuse needs attribution, and the record asks users to register their use. No personal data. The file is too large for GitHub, so the repository explains how to download it.</dd>
  </dl>
  <h3>Important variables</h3>
  <div class="scroll"><table>
    <thead><tr><th>Column</th><th>Meaning</th><th>Role here</th></tr></thead>
    <tbody>
      <tr><td><code>mut_type</code></td><td>The change, such as <code>W34D</code>; also <code>wt</code>, insertions, deletions and multiple substitutions</td><td>Filter; source of the swap and position</td></tr>
      <tr><td><code>aa_seq</code></td><td>Sequence of the variant</td><td>Source of length and neighbouring residues</td></tr>
      <tr><td><code>ddG_ML</code></td><td>Stability change caused by the mutation, kcal/mol</td><td>Target</td></tr>
      <tr><td><code>WT_name</code></td><td>The unmutated protein the variant comes from</td><td>Error analysis</td></tr>
      <tr><td><code>WT_cluster</code></td><td>Family of related proteins: a number for natural ones, a fold label such as <code>EEHH</code> for designed ones</td><td>Group for splitting</td></tr>
      <tr><td><code>deltaG</code>, <code>dG_ML</code>, <code>log10_K50_*</code> and the confidence columns</td><td>The experiment's own measurements</td><td>Never used: the target is computed from them</td></tr>
    </tbody>
    <caption>The file has no data dictionary. Meanings come from the column names, the paper, and checks on the rows.</caption>
  </table></div>
  <h3>Quality issues relevant to the question</h3>
  <ul id="quality"></ul>
  """ + HIST_FIG + """
</section>
<section id="preprocessing">
  <h2><span class="n">03</span>Preprocessing</h2>
  <div class="scroll"><table><thead><tr><th>Step</th><th>Why</th><th>Learned from</th></tr></thead><tbody id="prepRows"></tbody></table></div>
  <h3>Leakage checks</h3>
  <ul>
    <li><strong>Target in disguise.</strong> <code>deltaG</code>, <code>dG_ML</code> and the protease columns are what ΔΔG is computed from. None is a feature.</li>
    <li><strong>Relatives across the split.</strong> Every protein appears hundreds of times, once per mutation, and related proteins share a cluster. Whole clusters are assigned to one split, and a check confirms no cluster is in two.</li>
    <li><strong>Fitted steps.</strong> One-hot encoding and scaling sit inside each model's pipeline, so they are fitted on training rows only.</li>
    <li><strong>One split, frozen.</strong> The split files were made once and every experiment uses the same ones.</li>
  </ul>
  <p class="repo" data-items="src/preprocess.py;src/add_sequence_features.py;src/split_data.py;src/features.py"></p>
</section>
<section id="methods">
  <h2><span class="n">04</span>Baselines and methods</h2>
  <p>Three baselines set the bar. A model is useful only if it beats the mean. It shows that position and context matter only if it also beats a table of the average effect of each swap, which knows nothing about position.</p>
  <h3>Feature sets</h3>
  <div class="scroll"><table><thead><tr><th>Set</th><th>Columns</th><th>Idea</th></tr></thead><tbody id="featRows"></tbody></table></div>
  <h3>How the features were chosen</h3>
  <p>The feature sets were built in rounds. Each round asked what the previous one was missing. Decisions below rest on validation scores and on the repeated grouped splits, which never touch the test set.</p>
  <ol class="story" id="fsStory"></ol>
  <h3>What each column means</h3>
  <div class="scroll"><table><thead><tr><th>Column</th><th>Set</th><th>Meaning</th></tr></thead><tbody id="glossRows"></tbody>
  <caption>Amino-acid columns are one-hot encoded; every numeric column is standard-scaled inside the model pipeline, using training rows only.</caption></table></div>
  <h3>Methods</h3>
  <div class="scroll"><table><thead><tr><th>Method</th><th>Core idea</th><th>Why include it</th><th>Feature sets</th></tr></thead><tbody id="methodRows"></tbody>
  <caption>Most learned models use one fixed setting (listed in src/evaluate.py). HistGradientBoosting and XGBoost were also tuned over a small grid on the richest feature set, by validation MAE (see Compare). The ensemble's members were fixed before any result was seen. Selection is between all of these, on validation MAE.</caption></table></div>
  """ + FS_FIG + """
  <p class="repo" data-items="src/train_models.py;src/compare_models.py;src/features.py;src/evaluate.py"></p>
</section>
<section id="evaluation">
  <h2><span class="n">05</span>Validation protocol</h2>
  <dl class="frame">
    <dt>Split</dt><dd><span data-v="config.split_note"></span>. Train <span data-v="summary.n_train"></span> mutations (<span data-v="summary.c_train"></span> clusters), validation <span data-v="summary.n_val"></span> (<span data-v="summary.c_val"></span>), test <span data-v="summary.n_test"></span> (<span data-v="summary.c_test"></span>). No cluster is in two splits.</dd>
    <dt>Selection</dt><dd>Every model is trained on the training rows. The learned model with the lowest validation MAE is selected.</dd>
    <dt>Final test</dt><dd>All models are scored on the test set once. The test set plays no part in selection.</dd>
    <dt>Uncertainty</dt><dd>A 95% interval for every test score, from resampling whole test clusters <span data-v="config.n_bootstrap"></span> times. Plus <span data-v="config.n_repeats"></span> repeated grouped splits of the training and validation rows, without touching the test set.</dd>
  </dl>
  <h3>Metrics</h3>
  <ul>
    <li><strong>MAE</strong> (headline and selection): the typical miss in kcal/mol.</li>
    <li><strong>RMSE</strong>: punishes large misses more.</li>
    <li><strong>R²</strong>: share of the variance explained; 0 means no better than the mean.</li>
    <li><strong>Spearman</strong>: does the model rank mutations correctly, overall and within each protein? That is what a screening step needs.</li>
  </ul>
  """ + LEAK_FIG + """
  <h3>Why split by cluster and not by row?</h3>
  <p>Mutations of one protein share its backbone, and related proteins share a cluster. With relatives on both sides, a model can recognise the family and the test stops measuring performance on new proteins. The chart above shows how big that effect is.</p>
  <h3>Why a substitution table as a baseline?</h3>
  <p>Much of what is known about mutations is which swap they are: anything to proline, or a large residue replaced by a small one, tends to hurt. A table of those averages is cheap. A model that cannot beat it has not used position or context.</p>
  <h3>Why intervals?</h3>
  <p>The test set has <span data-v="summary.c_test"></span> independent clusters, far fewer than its rows suggest. A different set of clusters would give a different score; the interval shows how different.</p>
  <p class="repo" data-items="src/split_data.py;src/evaluate.py"></p>
</section>
<section id="brief">
  <h2><span class="n">06</span>Results in brief</h2>
  <p id="briefText"></p>
  <div class="claims" data-claims="results"></div>
  <p><a href="results.html">All results</a> · <a href="compare.html">Compare methods</a> · <a href="errors.html">Where the predictions fail</a></p>
</section>
<section id="ai">
  <h2><span class="n">07</span>AI use and contributions</h2>
  <p>See the <a href="https://github.com/marxy-karly/mutation-stability#ai-use-disclosure">AI-use disclosure in the README</a>.</p>
</section>
<section id="explore-cards"><h2>Explore</h2>
<div class="cards">
  <a class="card" href="tour.html"><b>Guided tour</b><span>Nine steps from the question to what the results can claim.</span></a>
  <a class="card" href="split.html"><b>Split</b><span>See how random folds flatter the models. Rerun it in your browser.</span></a>
  <a class="card" href="errors.html"><b>Errors</b><span>Where the selected model fails, and what it relies on.</span></a>
  <a class="card" href="compare.html"><b>Compare</b><span>Any two methods side by side, and the full ranking.</span></a>
  <a class="card" href="results.html"><b>Results</b><span>Every score, interval and limitation.</span></a>
  <a class="card" href="reproduce.html"><b>Reproduce</b><span>Setup, data access and the order to run the scripts.</span></a>
</div></section>"""

TOUR = head("Guided tour", "Nine steps through the project", "Each step shows one piece of evidence from the same pipeline output as the rest of the site.") + """
<div class="tour-bar"><div class="tour-track"><div class="tour-fill" id="tourFill"></div></div><p class="tour-count" id="tourCount"></p></div>
<section class="tstep" id="step-1"><h2><span class="n">01</span>The question</h2><p>Testing whether one amino-acid change helps or hurts a protein takes an experiment. A model that ranks candidate mutations first would save work, if it works on new proteins.</p>
<dl class="frame"><dt>One row</dt><dd>One single substitution in one protein, such as W34D.</dd><dt>Target</dt><dd>ΔΔG in kcal/mol. Positive means more stable.</dd><dt>Inputs</dt><dd>The swap, its position, the sequence around it, and how buried the residue is in a predicted structure.</dd><dt>Test</dt><dd>Protein families the model never saw.</dd></dl></section>
<section class="tstep" id="step-2"><h2><span class="n">02</span>One row of the data</h2><p>Five random rows. The sequence shows six residues either side of the mutated one, in bold.</p>
<div class="scroll"><table><thead><tr><th>Mutation</th><th>Sequence around it</th><th class="r">ΔΔG</th><th>Cluster</th></tr></thead><tbody id="tourRows"></tbody></table></div><button type="button" class="ghost" id="tourShuffle">Show five other rows</button></section>
<section class="tstep" id="step-3"><h2><span class="n">03</span>What the target looks like</h2>""" + HIST_FIG + """</section>
<section class="tstep" id="step-4"><h2><span class="n">04</span>What is wrong with the data</h2><ul id="quality"></ul></section>
<section class="tstep" id="step-5"><h2><span class="n">05</span>Why the split follows protein families</h2>""" + LEAK_FIG + """</section>
<section class="tstep" id="step-6"><h2><span class="n">06</span>The bar to beat</h2><p>Three baselines that use almost no information, and the selected model, on the test set.</p>
<div class="scroll"><table><thead><tr><th>Method</th><th class="r">Test MAE</th><th class="r">95% interval</th></tr></thead><tbody id="tourBase"></tbody></table></div></section>
<section class="tstep" id="step-7"><h2><span class="n">07</span>The models against the baselines</h2><p id="tourVerdict"></p>""" + FOREST + """</section>
<section class="tstep" id="step-8"><h2><span class="n">08</span>Where the predictions fail</h2><p>Each point is a test mutation. A perfect model would put every point on the dashed line.</p>
<div class="figure"><div class="controls"><label for="exModel">Method<select id="exModel"></select></label></div><div class="legend"><span><i class="s1"></i>Designed</span><span><i class="s2"></i>Natural</span></div><div id="chartScatter"></div><p class="readout" id="exReadout"></p></div></section>
<section class="tstep" id="step-9"><h2><span class="n">09</span>What the results can and cannot claim</h2><p id="conclusionText"></p><div class="claims" data-claims="results"></div><p><a href="results.html">Full results</a> · <a href="reproduce.html">Reproduce it</a></p></section>
<div class="tour-nav"><a class="btn ghost" id="tourPrev" href="#step-1">Back</a><a class="btn" id="tourNext" href="#step-2">Next</a></div>"""

SPLIT = head("Interactive", "Why the split matters", "The same models, scored two ways. Most of the apparent skill under random folds comes from relatives of the test rows sitting in training.") + """
<section id="leak"><h2>Cross-validation on the training data</h2>""" + LEAK_FIG + """<p class="repo" data-items="src/evaluate.py"></p></section>
<section id="sandbox">
  <h2>Run it yourself</h2>
  <p>This panel reruns the comparison in your browser on a random <span data-v="config.sandbox_rows"></span> rows: the mean, the substitution table and ridge on F1, with 5-fold cross-validation. Change the penalty, reshuffle the folds, or load your own CSV.</p>
  <div class="figure">
    <div class="controls">
      <label for="sbAlpha">Ridge penalty α: <span class="num" id="sbAlphaVal"></span><input type="range" id="sbAlpha" min="-1" max="4" step="0.25" value="1"></label>
      <label for="sbSeed">Shuffle seed<input type="number" id="sbSeed" value="0" min="0" step="1" style="width:6rem"></label>
      <button type="button" id="sbShuffle">Reshuffle folds</button>
    </div>
    <div class="legend"><span><i class="s1"></i>Folds keep groups together</span><span><i class="s2"></i>Random folds</span></div>
    <div id="chartSandbox"></div>
    <p class="readout" id="sbReadout"></p>
    <details><summary>Use a different dataset</summary>
      <p>Any CSV with a sequence column, a mutation column (single substitutions such as <code>W34D</code>), a numeric target and a group column. The file stays in your browser. Up to 6,000 rows are used.</p>
      <div class="controls"><label for="sbFile">CSV file<input type="file" id="sbFile" accept=".csv,text/csv"></label></div>
      <div class="controls" id="sbMap" hidden>
        <label for="sbSeqCol">Sequence<select id="sbSeqCol"></select></label><label for="sbMutCol">Mutation<select id="sbMutCol"></select></label>
        <label for="sbTargetCol">Target<select id="sbTargetCol"></select></label><label for="sbGroupCol">Group<select id="sbGroupCol"></select></label>
        <button type="button" id="sbRun">Run on this file</button><button type="button" class="ghost" id="sbReset">Back to the project data</button>
      </div>
      <p class="error" id="sbError" role="alert"></p>
    </details>
  </div>
  <p>Things to try: reshuffle a few times and watch how much the grouped scores move. That spread is why the results come with intervals.</p>
</section>
<section id="claims"><h2>What this supports</h2><div class="claims" data-claims="split"></div></section>"""

ERRORS = head("Interactive", "Where the predictions fail", "Pick a method and a kind of substitution. Below are the error slices and what the selected model relies on.") + """
<section id="explorer">
  <div class="figure">
    <div class="controls"><label for="exModel">Method<select id="exModel"></select></label><label for="exStratum">Substitution type<select id="exStratum"></select></label></div>
    <div class="legend"><span><i class="s1"></i>Designed</span><span><i class="s2"></i>Natural</span></div>
    <div id="chartScatter"></div>
    <p class="readout" id="exReadout"></p>
  </div>
  <h3>Largest misses in this view</h3>
  <div class="scroll"><table><thead><tr><th>Mutation in protein</th><th>Type</th><th>Cluster</th><th class="r">Measured</th><th class="r">Predicted</th><th class="r">Error</th></tr></thead><tbody id="worstRows"></tbody></table></div>
</section>
<section id="slices">
  <h2>Where the selected model's error concentrates</h2>
  <div class="scroll"><table><thead><tr><th>Slice of the test set</th><th class="r">Mutations</th><th class="r">MAE</th><th class="r">Mean signed error</th><th class="r">Table MAE</th><th class="r">Mean MAE</th></tr></thead><tbody id="sliceRows"></tbody>
  <caption>Signed error is predicted minus measured: negative means the model predicts too low. "Table" is the substitution-table baseline.</caption></table></div>
  <p id="errorNote"></p>
  <h3>By test cluster</h3>
  <div class="scroll"><table><thead><tr><th>Cluster</th><th>Origin</th><th class="r">Proteins</th><th class="r">Mutations</th><th class="r">MAE</th><th class="r">Signed error</th><th class="r">Mean MAE</th></tr></thead><tbody id="clusterRows"></tbody></table></div>
</section>
<section id="importance">
  <h2>What the selected model relies on</h2>
  <div class="figure"><h4>Permutation importance by kind of information</h4><p class="sub">How much the validation MAE rises when every column of one kind is shuffled together (3 repeats).</p><div id="chartImportanceGroups"></div></div>
  <p id="importanceGroupNote"></p>
  <div class="figure"><h4>Permutation importance by column</h4><p class="sub">How much the validation MAE rises when one input column is shuffled (3 repeats).</p><div id="chartImportance"></div></div>
  <p id="importanceNote"></p>
  <p class="repo" data-items="src/evaluate.py"></p>
</section>
<section id="claims"><h2>What this supports</h2><div class="claims" data-claims="errors"></div></section>"""

COMPARE = head("Interactive", "Compare methods", "Pick any two methods, see every metric side by side, and see the ranking the selection was based on.") + """
<section id="compare">
  <div class="figure">
    <div class="controls"><label for="cmpA">Method A<select id="cmpA"></select></label><label for="cmpB">Method B<select id="cmpB"></select></label></div>
    <div class="scroll"><table><thead><tr><th>Measure</th><th class="r" id="cmpHA">A</th><th class="r" id="cmpHB">B</th><th class="r">B minus A</th></tr></thead><tbody id="cmpRows"></tbody>
    <caption>MAE and RMSE in kcal/mol. For errors, a positive "B minus A" means B is worse.</caption></table></div>
    <p class="readout" id="cmpReadout"></p>
  </div>
  <h2>Ranking by validation MAE</h2>
  <p>This is the ranking the model was selected from. Test MAE is shown for reference only.</p>
  <div class="scroll"><table><thead><tr><th class="r">Rank</th><th>Method</th><th class="r">Validation MAE</th><th class="r">Test MAE</th><th class="r">Mean MAE over repeats</th></tr></thead><tbody id="rankRows"></tbody></table></div>
</section>
<section id="tuning"><h2>Tuning grid</h2>
  <p>Each setting was trained on the training rows and scored on validation. The one with the lowest validation MAE was kept. The test set was not used.</p>
  <div class="scroll"><table><thead><tr><th>Algorithm</th><th>Setting</th><th>Features</th><th class="r">Validation MAE</th></tr></thead><tbody id="tuneRows"></tbody></table></div>
</section>
<section id="features"><h2>What each feature set adds</h2>""" + FS_FIG + """</section>
<section id="claims"><h2>What this supports</h2><div class="claims" data-claims="compare"></div></section>"""

RESULTS = head("Results", "Every result in one place", "Test scores with intervals, the repeated splits, the limitations and the conclusion.") + """
<section id="results">
  <h2>Results</h2>
  <p id="verdict"></p>
  <div class="scroll"><table>
    <thead><tr><th>Method</th><th class="r">Val MAE</th><th class="r">Test MAE</th><th class="r">95% interval</th><th class="r">Test RMSE</th><th class="r">Test R²</th><th class="r">Spearman</th><th class="r">Within-protein Spearman</th></tr></thead>
    <tbody id="resultRows"></tbody><caption id="resultCaption"></caption></table></div>
  """ + FOREST + """
  <p id="contextNote"></p>
  <p class="repo" data-items="src/evaluate.py;results/evaluation.json;results/test_predictions.csv.gz"></p>
</section>
<section id="limitations"><h2>Limitations and sources of prediction error</h2>
  <p id="limIntro"></p>
  <ol class="story" id="limits"></ol>
  <p id="limNote"></p>
  <h3>Other limitations of the method</h3>
  <ul id="limOther"></ul>
</section>
<section id="conclusion">
  <h2>Conclusion</h2>
  <p id="conclusionText"></p>
  <p>Two findings stand out. First, information about where the residue sits in the folded protein mattered far more than extra sequence detail or tuning. Second, the split decides how good a model looks: the same models look better under random folds than under folds that respect protein families.</p>
  <p id="nextStep"></p>
  <div class="claims" data-claims="results"></div>
</section>"""

REPRODUCE = head("Reproduce", "Check the work", "The repository holds the code and the full results. This page says how to rerun them.") + """
<section id="reproduce">
  <p id="repoLine"></p>
  <pre>pip install -r requirements.txt
python src/preprocess.py              # single substitutions with a numeric ddG_ML
python src/split_data.py              # grouped train / validation / test split
python src/train_models.py            # first baselines (F1)
python src/add_sequence_features.py   # neighbouring residues (F2)
python src/update_splits.py           # same split for the F2 files
python src/compare_models.py          # F1 vs F2 comparison
python src/features.py                # optional: check how many proteins have a structure
python src/evaluate.py                # everything on this site (adds F3 and F4)
python -m http.server -d docs         # open http://localhost:8000</pre>
  <dl class="frame">
    <dt>Data</dt><dd>Download <code>Processed_K50_dG_datasets.zip</code> from <a href="https://zenodo.org/records/7401275">Zenodo record 7401275</a> and unzip <code>K50_dG_Dataset1_Dataset2.csv</code> where <code>src/preprocess.py</code> expects it (see README). It is not committed: GitHub rejects files over 100 MB.</dd>
    <dt>Structures</dt><dd>Download <code>AlphaFold_model_PDBs.zip</code> from the same record and unzip it into <code>data/structures/</code>. <code>src/features.py</code> matches models to proteins by sequence. Without it, F4 is skipped and everything else still runs.</dd>
    <dt>Split</dt><dd><span data-v="config.split_note"></span>. Cluster membership of every row is in the split files.</dd>
    <dt>Environment</dt><dd id="envLine"></dd>
    <dt>Citation</dt><dd>Tsuboyama, K. et al. Mega-scale experimental analysis of protein folding stability in biology and design. <em>Nature</em> 620, 434–444 (2023). Data: doi:10.5281/zenodo.7401275.</dd>
    <dt>AI use</dt><dd>See the README.</dd>
  </dl>
</section>"""

# author, contribution, AI-use text and repository link come from config.json, so editing them needs
# only this script, not a new run of src/evaluate.py
import json
_cfg = json.loads((DOCS.parent / "config.json").read_text(encoding="utf-8"))
_about = {k: _cfg.get(k, "") for k in ("author", "contribution", "ai_use", "repo_url")}
(DOCS / "about.js").write_text("window.ABOUT = " + json.dumps(_about, ensure_ascii=False) + ";\n", encoding="utf-8")

PAGES = {"index": ("Mutation Stability", INDEX), "tour": ("Guided Tour", TOUR), "split": ("The Split", SPLIT), "errors": ("Error Analysis", ERRORS),
         "compare": ("Compare Methods", COMPARE), "results": ("Results", RESULTS), "reproduce": ("Reproduce", REPRODUCE)}
for pid, (title, body) in PAGES.items():
    (DOCS / f"{pid}.html").write_text(page(pid, title, body), encoding="utf-8")
print("wrote", len(PAGES), "pages")
