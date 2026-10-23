/* Verify the selected published HTML snapshot without rebuilding or changing it. */
'use strict';
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const crypto = require('node:crypto');
const vm = require('node:vm');
const G = require('./graph-model.js');
const root = path.join(__dirname, '..');
const lockPath = path.join(__dirname, 'published-viewer.json');
const hash = bytes => crypto.createHash('sha256').update(bytes).digest('hex');

function inspect() {
  const filename = 'InterpreterPLATO-Visualizer.html';
  const bytes = fs.readFileSync(path.join(root, filename)), html = bytes.toString('utf8');
  assert.equal(/<script[^>]+src=|<link[^>]+stylesheet/i.test(html), false, 'The published HTML must be self-contained');
  const scripts = [...html.matchAll(/<script>\n([\s\S]*?)\n<\/script>/g)];
  assert.equal(scripts.length, 5, 'Expected the selected five-case viewer snapshot');
  scripts.forEach(([, source]) => new vm.Script(source));
  const match = html.match(/^window\.InterpreterRunExamples = (.*);\r?$/m);
  assert.ok(match, 'The selected snapshot must contain its embedded run suite');
  const runs = JSON.parse(match[1]);
  assert.equal(runs.length, 5);
  assert.deepEqual(runs.map(run => run.source.case_id), ['SCORERI-DEV-01', 'SCORERI-DEV-02', 'SCORERI-DEV-03', 'SCORERI-DEV-05', 'SCORERI-DEV-06']);
  const cases = runs.map(run => {
    assert.equal(run.schema, 'interpreterplato-local-run/v1');
    assert.match(run.source.sha256, /^[a-f0-9]{64}$/);
    assert.equal(run.payload.run_mode, 'test');
    assert.equal(run.payload.generated_from_test_fixture, false);
    assert.notEqual(run.payload.formal_result_eligible, true);
    const model = G.build(run.payload);
    for (const id of run.presentation.pathIds) {
      assert.equal(model.paths.find(item => item.id === id)?.valid, true, 'A selected path must exist in its supplied graph');
    }
    return { case_id: run.source.case_id, source_sha256: run.source.sha256,
      payload_sha256: hash(JSON.stringify(run.payload)), nodes: model.nodeCount,
      edges: model.elements.filter(item => item.group === 'edges').length };
  });
  return { schema: 'interpreterplato-published-viewer/v1', filename,
    bytes: bytes.length, sha256: hash(bytes), cases };
}
const snapshot = inspect();
if (process.argv.length === 3 && process.argv[2] === '--update-lock') {
  fs.writeFileSync(lockPath, JSON.stringify(snapshot, null, 2) + '\n');
} else {
  assert.equal(process.argv.length, 2, 'Usage: node visualizer/verify-published.cjs [--update-lock]');
  assert.deepEqual(snapshot, JSON.parse(fs.readFileSync(lockPath, 'utf8')), 'The published viewer differs from its reviewed snapshot');
}
console.log(JSON.stringify({ verified: true, ...snapshot }, null, 2));
