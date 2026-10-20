/* Hand-authored synthetic concept fixture. No biological computation or real evidence. */
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory();
  else root.InterpreterExample = factory();
})(typeof globalThis !== 'undefined' ? globalThis : this, function () {
  'use strict';
  const geneNames = ['SYN_A', 'SYN_B', 'SYN_SHARED', 'SYN_BRIDGE_1', 'SYN_BRIDGE_2', 'SYN_DISEASE', 'SYN_UNMAPPED', 'SYN_ISOLATED'];
  const nodes = [
    { id: 'drug:Demo-A', type: 'drug', name: 'Demo-A' }, { id: 'drug:Demo-B', type: 'drug', name: 'Demo-B' },
    ...geneNames.map(name => ({ id: 'gene:' + name, type: 'gene', name, synthetic_fixture: true })),
    { id: 'pathway:SYN_PATHWAY', type: 'pathway', name: '合成通路集合' },
    { id: 'disease_module:SYN_CONTEXT', type: 'disease_module', name: '合成疾病基因集合' },
    { id: 'evidence:fixture', type: 'prediction_evidence', name: '合成证据占位' },
    { id: 'mechanism:baseline', type: 'mechanism', name: '报告摘要结构' },
  ];
  const ppiPairs = [['SYN_A', 'SYN_BRIDGE_1'], ['SYN_BRIDGE_1', 'SYN_SHARED'], ['SYN_SHARED', 'SYN_B'],
    ['SYN_B', 'SYN_BRIDGE_2'], ['SYN_BRIDGE_2', 'SYN_DISEASE'], ['SYN_BRIDGE_1', 'SYN_DISEASE']];
  const targetSets = [['Demo-A', ['SYN_A', 'SYN_SHARED', 'SYN_UNMAPPED']], ['Demo-B', ['SYN_B', 'SYN_SHARED', 'SYN_ISOLATED']]];
  const edges = [
    ...ppiPairs.map(([left, right]) => ({ source: 'gene:' + left, target: 'gene:' + right, type: 'PPI', directed: false,
      source_dataset: 'synthetic_fixture_only_not_observed_ppi', confidence: null, confidence_available: false, synthetic_fixture: true })),
    ...targetSets.flatMap(([drug, genes]) => genes.map(gene => ({ source: 'drug:' + drug, target: 'gene:' + gene,
      type: 'PLATO_predicted_target', provenance: { synthetic_fixture: true, actual_plato_prediction: false } }))),
    ...['SYN_A', 'SYN_B', 'SYN_SHARED'].map(gene => ({ source: 'gene:' + gene, target: 'pathway:SYN_PATHWAY', type: 'belongs_to_pathway' })),
    { source: 'gene:SYN_DISEASE', target: 'disease_module:SYN_CONTEXT', type: 'belongs_to_disease_module' },
    { source: 'evidence:fixture', target: 'gene:SYN_A', type: 'supports_prediction' },
    ...['Demo-A', 'Demo-B'].map(drug => ({ source: 'mechanism:baseline', target: 'drug:' + drug, type: 'summarizes_combination_member', evidence_type: 'report_structure_only' })),
  ];
  return {
    result_schema_version: '5.4', drug1_name: 'Demo-A', drug2_name: 'Demo-B', cell_line: 'SYNTHETIC_CELL',
    run_mode: 'test', generated_from_test_fixture: true, formal_result_eligible: false,
    mechanism_type: 'not_computed', predicted_synergy: null,
    metadata: { generated_from_test_fixture: true, conceptual_example: true, calculation_performed: false,
      notice: '概念示意（非真实计算结果）。全部节点和关系为手工构造的合成测试数据，不能用于科研证据。' },
    drug1_targets: { targets: targetSets[0][1].map(target_id => ({ target_id, evidence_status: 'synthetic_fixture' })) },
    drug2_targets: { targets: targetSets[1][1].map(target_id => ({ target_id, evidence_status: 'synthetic_fixture' })) },
    evidence_graphs: { explanation_evidence_graph: { nodes, edges, semantics: 'synthetic concept fixture; not biological evidence' } },
    local_evidence: {
      cross_drug_paths: [
        { drug: 'Demo-A', destination_drug: 'Demo-B', path: ['SYN_A', 'SYN_BRIDGE_1', 'SYN_SHARED', 'SYN_B'], length: 3, path_kind: 'cross_drug_target' },
        { drug: 'Demo-A', destination_drug: 'Demo-B', path: ['SYN_SHARED'], length: 0, shared_target: true },
      ],
      top_paths: [
        { drug: 'Demo-A', path: ['SYN_A', 'SYN_BRIDGE_1', 'SYN_DISEASE'], length: 2, path_kind: 'target_to_disease' },
        { drug: 'Demo-B', path: ['SYN_B', 'SYN_BRIDGE_2', 'SYN_DISEASE'], length: 2, path_kind: 'target_to_disease' },
      ],
    },
    network_connectivity: {
      status: 'synthetic_example_not_computed',
      drugs: [
        { drug: 'Demo-A', targets: [{ target: 'SYN_A', status: 'shown' }, { target: 'SYN_SHARED', status: 'shown' }, { target: 'SYN_UNMAPPED', status: 'target_not_in_ppi' }] },
        { drug: 'Demo-B', targets: [{ target: 'SYN_B', status: 'shown' }, { target: 'SYN_SHARED', status: 'display_limit' }, { target: 'SYN_ISOLATED', status: 'disconnected_from_disease' }] },
      ],
      note: '这些诊断为界面测试占位，不是解释器运行结果。',
    },
    network_statistics: { valid: false, quality_status: 'not_evaluated', reason: 'conceptual_example_no_calculation' },
  };
});
