/* Prepare a local run for the offline viewer. Does not query services or change results. */
'use strict';
const fs = require('node:fs');
const path = require('node:path');
const crypto = require('node:crypto');
const G = require('./graph-model.js');

function prepareRun(file, requestedPaths = []) {
  const bytes = fs.readFileSync(file);
  if (bytes.length > G.LIMITS.fileBytes) throw new Error('Local run exceeds the 20 MiB viewer limit');
  const records = G.records(JSON.parse(bytes.toString('utf8').replace(/^\uFEFF/, '')));
  if (records.length !== 1) throw new Error('Choose one completed run for the embedded example');
  const payload = records[0], model = G.build(payload);
  if (payload.error || G.isSynthetic(payload) || !['test', 'research_hypothesis'].includes(payload.run_mode)) {
    throw new Error('A run example requires a completed non-synthetic test or research-hypothesis result');
  }
  const graph = payload.evidence_graphs && payload.evidence_graphs.explanation_evidence_graph;
  if (!graph) throw new Error('The run has no recorded evidence graph');
  const drugNames = [payload.drug1_name, payload.drug2_name];
  const nodes = new Map(graph.nodes.map(node => [node.id, node]));
  const drugNodes = drugNames.map(name => graph.nodes.find(node => node.type === 'drug' && node.name === name));
  if (drugNodes.some(node => !node)) throw new Error('Both named drugs must exist in the recorded graph');
  function linked(type, a, b) {
    return graph.edges.some(edge => G.relationship(edge.type) === type &&
      ((edge.source === a && edge.target === b) || (edge.source === b && edge.target === a)));
  }
  function contextFor(row) {
    const last = 'gene:' + row.genes[row.genes.length - 1];
    return graph.nodes.find(node => node.type === 'disease_module' && linked('disease', last, node.id));
  }
  const eligible = model.paths.filter(row => row.valid && row.kind === 'disease' &&
    drugNames.includes(row.raw.drug) && contextFor(row) &&
    linked('target', drugNodes[drugNames.indexOf(row.raw.drug)].id, 'gene:' + row.genes[0]));
  const paths = requestedPaths.length ? requestedPaths.map(id => eligible.find(row => row.id === id))
    : drugNames.map(name => eligible.find(row => row.raw.drug === name && !row.zeroHop) || eligible.find(row => row.raw.drug === name));
  if (paths.length !== 2 || paths.some(row => !row) || new Set(paths.map(row => row.raw.drug)).size !== 2) {
    throw new Error('Choose one recorded disease path for each drug; missing edges are never supplied');
  }
  const contexts = paths.map(contextFor);
  if (contexts[0].id !== contexts[1].id) throw new Error('The selected paths must lead to the same recorded disease module');
  const ids = new Set(drugNodes.map(node => node.id));
  const positions = {}, compactPositions = {}, targets = [], members = [];
  const maxGenes = Math.max(...paths.map(row => row.genes.length));
  for (const [i, drug] of drugNodes.entries()) {
    const row = paths.find(item => item.raw.drug === drugNames[i]);
    const genes = row.genes.map(gene => 'gene:' + gene), y = i ? 110 : -110, x = i ? 100 : -100;
    positions[drug.id] = [-430, y]; compactPositions[drug.id] = [x, -240];
    targets.push(genes[0]); members.push(genes[genes.length - 1]);
    genes.forEach((id, index) => {
      if (!nodes.has(id)) throw new Error('A focus node is missing from the evidence graph');
      ids.add(id);
      // Position is editorial only; all nodes and edges remain the original recorded data.
      if (!positions[id]) {
        positions[id] = [genes.length === 1 ? 160 : -130 + index * 290, y];
        compactPositions[id] = [x, genes.length === 1 ? 70 : -80 + index * 150];
      }
    });
  }
  const endpoint = contexts[0].id; ids.add(endpoint);
  positions[endpoint] = [-130 + Math.max(1, maxGenes - 1) * 290 + 360, 0]; compactPositions[endpoint] = [0, -80 + Math.max(1, maxGenes - 1) * 150 + 170];
  const labels = Array.isArray(payload.disease_context_labels) ? payload.disease_context_labels.join(' / ') : '';
  const source = payload.disease_context_source || '输入记录中的疾病集合';
  const contextLabel = labels === 'pan_cancer_tier1' ? '泛癌 Tier 1' : labels || contexts[0].name;
  const caseId = payload.audit && payload.audit.engineering_case && payload.audit.engineering_case.case_id || 'LOCAL-RUN';
  return {
    schema: 'interpreterplato-local-run/v1',
    source: { case_id: caseId, filename: path.basename(file), sha256: crypto.createHash('sha256').update(bytes).digest('hex'), timestamp: payload.timestamp || null },
    presentation: { nodes: [...ids], pathIds: paths.map(row => row.id), endpoint, targets, members, positions, compactPositions,
      contextLabel, contextSource: source, contextSize: payload.disease_context_size || null,
      contextDisplayLabel: '疾病模块 / ' + contextLabel + '\n' + (source === 'COSMIC Cancer Gene Census' ? 'COSMIC CGC' : source) + ' · ' + (payload.disease_context_size || '?') + ' 个基因' },
    payload,
  };
}
module.exports = { prepareRun };
