/* Synthetic structural fixture for frontend tests; never distributed as a run example. */
'use strict';
const fs = require('node:fs');
const os = require('node:os');
const path = require('node:path');
const assert = require('node:assert/strict');
function runFixture() {
  return {
    result_schema_version: '5.4', drug1_name: 'Fixture-A', drug2_name: 'Fixture-B', cell_line: 'TEST-CELL',
    run_mode: 'test', generated_from_test_fixture: false, formal_result_eligible: false,
    mechanism_type: 'network_close_distinct_processes', predicted_synergy: 12.34567,
    disease_context_source: 'Test-only gene set', disease_context_labels: ['fixture_context'], disease_context_size: 123,
    timestamp: '2026-10-13T00:00:00Z', audit: { engineering_case: { case_id: 'TEST-CASE' } },
    network_statistics: { valid: true, quality_status: 'adequate', statistically_supported: false, p_A_adjusted: 0.8, p_B_adjusted: 1, n_randomizations: 2000 },
    warning_flags: ['fixture_warning'],
    evidence_graphs: { explanation_evidence_graph: {
      nodes: [
        { id: 'drug:A', type: 'drug', name: 'Fixture-A' }, { id: 'drug:B', type: 'drug', name: 'Fixture-B' },
        ...['A0', 'A1', 'B0', 'B1'].map(name => ({ id: 'gene:' + name, type: 'gene', name })),
        { id: 'disease_module:fixture', type: 'disease_module', name: 'fixture_context' },
        { id: 'pathway:fixture', type: 'pathway', name: 'fixture pathway' },
        { id: 'mechanism:fixture', type: 'mechanism', name: 'fixture report' },
      ], edges: [
        { source: 'drug:A', target: 'gene:A0', type: 'PLATO_predicted_target' },
        { source: 'drug:B', target: 'gene:B0', type: 'PLATO_predicted_target' },
        { source: 'gene:A1', target: 'gene:A0', type: 'PPI', directed: false },
        { source: 'gene:B0', target: 'gene:B1', type: 'PPI', directed: false },
        { source: 'gene:A1', target: 'disease_module:fixture', type: 'belongs_to_disease_module' },
        { source: 'gene:B1', target: 'disease_module:fixture', type: 'belongs_to_disease_module' },
        { source: 'gene:A0', target: 'pathway:fixture', type: 'belongs_to_pathway' },
        { source: 'drug:A', target: 'mechanism:fixture', type: 'supports_prediction' },
      ],
    } },
    local_evidence: { top_paths: [{ drug: 'Fixture-A', path: ['A0', 'A1'], length: 1 }, { drug: 'Fixture-B', path: ['B0', 'B1'], length: 1 }] },
  };
}
function withRunFile(callback) {
  const temporary = fs.mkdtempSync(path.join(os.tmpdir(), 'plato-local-run-test-'));
  const file = path.join(temporary, 'explanation.json'), payload = runFixture();
  fs.writeFileSync(file, JSON.stringify(payload));
  try { return callback({ temporary, file, payload }); }
  finally {
    assert.equal(path.dirname(path.resolve(temporary)), path.resolve(os.tmpdir()));
    assert.ok(path.basename(temporary).startsWith('plato-local-run-test-'));
    fs.rmSync(temporary, { recursive: true, force: true });
  }
}
module.exports = { runFixture, withRunFile };
