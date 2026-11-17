// Checks that the website's in-browser F1 features and ridge match src/evaluate.py.  Run: node tests/check_js.mjs
import fs from 'node:fs';
import vm from 'node:vm';
const stub = () => ({ innerHTML: '', textContent: '', value: '1', style: {}, hidden: false, addEventListener() {}, getAttribute() { return null; },
  getBoundingClientRect() { return { width: 640, height: 10 }; }, classList: { toggle() {} }, insertAdjacentHTML() {} });
const els = {};
const ctx = { location: { hash: '' }, window: { addEventListener() {}, innerWidth: 1000, innerHeight: 800, scrollTo() {} }, console, Set, Math, Array, Object, Number, String, parseInt, parseFloat, isFinite, isNaN,
  document: { title: '', createElement: () => stub(), getElementById: id => (els[id] ??= stub()), querySelectorAll: () => [], addEventListener() {}, body: stub() } };
ctx.window.window = ctx.window;
vm.createContext(ctx);
vm.runInContext(fs.readFileSync(new URL('../docs/results.js', import.meta.url), 'utf8'), ctx);
vm.runInContext('var RESULTS = window.RESULTS;' + fs.readFileSync(new URL('../docs/app.js', import.meta.url), 'utf8'), ctx);
const R = ctx.window.RESULTS, sb = ctx.window.__sandbox, rows = R.sandbox.rows;
let bad = 0, maxd = 0;
R.sandbox.x_check.forEach((x, i) => { const f = sb.features(rows[i].seq, rows[i].mut); if (f.length !== x.length) bad++; f.forEach((v, j) => { maxd = Math.max(maxd, Math.abs(v - x[j])); }); });
console.log('features: max |browser - python| over', R.sandbox.x_check.length, 'rows =', maxd);
if (maxd > 1e-9) { console.error('FEATURE MISMATCH'); bad++; }
const tr = [], te = []; rows.forEach((r, i) => (r.test ? te : tr).push(i));
const fit = sb.ridge(sb.data.X, sb.data.y, tr, R.sandbox.check.alpha);
const js = sb.maeOf(te.map(i => sb.data.y[i]), te.map(fit));
console.log('ridge alpha', R.sandbox.check.alpha, '| scikit-learn test MAE', R.sandbox.check.test_mae.toFixed(6), '| browser', js.toFixed(6));
if (Math.abs(js - R.sandbox.check.test_mae) > 1e-5) { console.error('RIDGE MISMATCH'); bad++; }
if (bad) process.exit(1);
console.log('OK: browser features and ridge match scikit-learn');
