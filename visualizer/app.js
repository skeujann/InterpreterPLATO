/* Browser-only viewer. No fetch, upload endpoint, telemetry, or model execution. */
(function () {
  'use strict';
  const G = window.InterpreterGraph;
  const $ = id => document.getElementById(id);
  const state = { cy: null, rows: [], model: null, source: '', selectedPath: null };
  const GROUP_LABELS = { gene: '基因／蛋白节点', drug: '药物', pathway: '通路集合', disease_module: '疾病基因集合', prediction_evidence: '预测证据记录', mechanism: '报告摘要结构', counterfactual: '敏感性分析记录' };
  const EDGE_LABELS = { ppi: '文件记录的无向 PPI', target: '候选靶点关联', disease: '疾病集合归属', pathway: '通路集合归属', report: '报告结构／证据支持', hypothesis: '规则生成假说', other: '未分类关系' };
  function element(tag, className, value) {
    const result = document.createElement(tag);
    if (className) result.className = className;
    if (value !== undefined) result.textContent = String(value);
    return result;
  }
  function error(message) { $('error-banner').textContent = message; $('error-banner').hidden = !message; }
  function download(blob, name) {
    const url = URL.createObjectURL(blob), link = element('a');
    link.href = url; link.download = name; document.body.appendChild(link); link.click(); link.remove();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
  }
  function pretty(value, limit) {
    const result = JSON.stringify(value, null, 2);
    return result.length > limit ? result.slice(0, limit) + '\n…（界面截断；完整属性保留在输入文件）' : result;
  }
  function tab(name) {
    for (const target of ['details', 'diagnostics']) {
      $('tab-' + target).setAttribute('aria-selected', String(name === target));
      $(target + '-panel').hidden = name !== target;
    }
  }
  function inspect(item) {
    const data = item.data();
    $('detail-title').textContent = item.isNode() ? data.label : EDGE_LABELS[data.category];
    $('detail-description').textContent = item.isNode()
      ? (GROUP_LABELS[data.type] || data.type) + ' · ' + data.originalId
      : data.type + (data.category === 'ppi' ? ' · 无方向；来源未在此独立核验' : ' · 关联或结构关系，不作因果箭头解读');
    $('detail-json').textContent = pretty(data.raw, 6500);
    $('detail-json').hidden = false;
  }
  function clearHighlight() {
    state.selectedPath = null;
    if (state.cy) state.cy.elements().removeClass('dimmed highlighted');
    document.querySelectorAll('.path-button').forEach(button => button.classList.remove('selected'));
  }
  function renderCounts() {
    if (!state.cy) return;
    const visible = state.cy.nodes().filter(node => !node.hasClass('filtered'));
    const ppi = state.cy.edges().filter(edge => edge.data('category') === 'ppi' && !edge.hasClass('filtered'));
    $('graph-counts').textContent = visible.length + ' / ' + state.cy.nodes().length + ' 节点 · ' + ppi.length + ' 条 PPI';
    $('graph-empty').hidden = visible.length > 0;
    $('graph-empty').querySelector('strong').textContent = state.cy.nodes().length ? '当前筛选未显示任何节点' : '该记录没有可展示的证据图';
    $('export-png').disabled = !visible.length;
  }
  function applyFilters() {
    if (!state.cy) return;
    const nodeGroups = new Set([...document.querySelectorAll('[data-node-group]:checked')].map(input => input.dataset.nodeGroup));
    const edgeGroups = new Set([...document.querySelectorAll('[data-edge-group]:checked')].map(input => input.dataset.edgeGroup));
    state.cy.batch(() => {
      state.cy.nodes().forEach(node => node.toggleClass('filtered', !nodeGroups.has(node.data('group'))));
      state.cy.edges().forEach(edge => {
        const category = edge.data('category');
        const group = ['disease', 'pathway'].includes(category) ? 'membership' : ['report', 'hypothesis', 'other'].includes(category) ? 'report' : category;
        edge.toggleClass('filtered', !edgeGroups.has(group) || edge.source().hasClass('filtered') || edge.target().hasClass('filtered'));
      });
    });
    renderCounts();
  }
  function layout() {
    if (!state.cy || !state.model.nodeCount) return;
    const name = $('layout-select').value;
    const options = { name, animate: false, fit: true, padding: 40, nodeDimensionsIncludeLabels: true };
    const visibleNodes = state.cy.nodes().filter(node => !node.hasClass('filtered')).length;
    if (name === 'cose' && visibleNodes > 400) {
      options.name = 'grid'; $('layout-select').value = 'grid';
      $('search-status').textContent = '节点较多，使用轻量网格布局。可筛选后重新布局。';
    }
    if (options.name === 'cose') Object.assign(options, { randomize: false, nodeRepulsion: 16000, idealEdgeLength: 110, numIter: 800 });
    if (name === 'breadthfirst') Object.assign(options, { directed: false, roots: state.cy.nodes().filter(node => node.data('type') === 'drug'), spacingFactor: 1.4 });
    if (name === 'concentric') Object.assign(options, { minNodeSpacing: 35, concentric: node => node.data('type') === 'drug' ? 10 : node.degree(), levelWidth: () => 2 });
    state.cy.elements().filter(item => !item.hasClass('filtered')).layout(options).run();
  }
  function selectPath(path, button) {
    clearHighlight(); $('node-search').value = '';
    if (!path.valid) return;
    state.selectedPath = path.id;
    const ids = new Set([...path.nodeIds, ...path.edgeIds]);
    state.cy.elements().addClass('dimmed');
    const selected = state.cy.elements().filter(item => ids.has(item.id()));
    selected.removeClass('dimmed').addClass('highlighted');
    button.classList.add('selected');
    const hidden = selected.filter(item => item.hasClass('filtered')).length;
    $('search-status').textContent = path.zeroHop ? '零跳记录：没有 PPI 边，也不会补画自环。' : '高亮 ' + path.edgeIds.length + ' 条文件中已有的 PPI 边。';
    if (hidden) $('search-status').textContent += ' 当前筛选隐藏了部分路径元素。';
    $('detail-title').textContent = path.kind === 'cross' ? '双药靶点连接' : '靶点—疾病路径';
    $('detail-description').textContent = path.zeroHop ? '零跳仅表示共享节点或集合成员关系。' : '路径仅用于高亮已记录的边，不增添网络关系。';
    $('detail-json').textContent = pretty(path.raw, 6500); $('detail-json').hidden = false; tab('details');
  }
  function renderPaths() {
    const list = $('path-list'); list.replaceChildren();
    if (!state.model.paths.length) { list.appendChild(element('p', 'small muted', '该记录未包含可展示的路径。')); return; }
    for (const [kind, title] of [['cross', '双药靶点连接'], ['disease', '靶点—疾病路径']]) {
      const paths = state.model.paths.filter(path => path.kind === kind);
      if (!paths.length) continue;
      list.appendChild(element('div', 'path-group-label', title + ' · ' + paths.length));
      for (const path of paths) {
        const button = element('button', 'path-button'); button.type = 'button'; button.disabled = !path.valid;
        button.appendChild(element('strong', '', path.raw.drug || title));
        button.appendChild(element('span', 'path-route', path.genes.join(' — ')));
        if (path.zeroHop) button.appendChild(element('span', 'path-zero', kind === 'cross' ? '共享节点零跳 · 0 条 PPI' : '疾病集合成员零跳 · 0 条 PPI'));
        else button.appendChild(element('span', 'path-zero', (path.genes.length - 1) + ' 跳 · 仅高亮文件中已有的边'));
        if (!path.valid) button.appendChild(element('span', 'path-issue', '无法高亮：' + path.missing.join('；')));
        button.addEventListener('click', () => selectPath(path, button)); list.appendChild(button);
      }
    }
  }
  function diagnostics() {
    const container = $('diagnostic-summary'); container.replaceChildren();
    for (const warning of state.model.warnings) container.appendChild(element('p', 'diagnostic-item', warning));
    const invalid = state.model.paths.filter(path => !path.valid).length;
    if (invalid) container.appendChild(element('p', 'diagnostic-item', invalid + ' 条路径与证据图不一致，已禁用高亮；没有自动补边。'));
    const payload = state.model.payload;
    if (payload.error) container.appendChild(element('p', 'diagnostic-item', '上游错误：' + String(payload.error)));
    $('diagnostic-json').textContent = pretty({
      network_connectivity: payload.network_connectivity || (payload.local_evidence && payload.local_evidence.connectivity) || { status: 'not_recorded' },
      warning_flags: payload.warning_flags || [],
      network_statistics_status: payload.network_statistics ? { valid: payload.network_statistics.valid, quality_status: payload.network_statistics.quality_status } : 'not_recorded',
      viewer_warnings: state.model.warnings,
    }, 16000);
  }
  function display(index, suppliedModel) {
    const payload = state.rows[index], model = suppliedModel || G.build(payload);
    state.model = model; state.selectedPath = null;
    if (state.cy) state.cy.destroy();
    $('node-search').value = ''; $('search-status').textContent = '输入名称可高亮匹配节点';
    $('case-title').textContent = (payload.drug1_name || payload.drug1 || '未标注药物 A') + ' + ' + (payload.drug2_name || payload.drug2 || '未标注药物 B');
    $('case-meta').textContent = (payload.cell_line ? payload.cell_line + ' · ' : '') + '记录 ' + (index + 1) + ' / ' + state.rows.length;
    $('source-name').textContent = state.source;
    const concept = payload.metadata && payload.metadata.conceptual_example === true;
    $('mode-banner').className = 'mode-banner ' + (model.synthetic || payload.run_mode === 'test' || payload.run_mode === 'demo' ? 'synthetic' : 'research');
    $('mode-title').textContent = concept ? '概念示意（非真实计算结果）' : model.mode;
    $('mode-description').textContent = concept ? '全部节点与关系为合成测试数据，不用于科研证据。' : '只展示输入文件记录的关系，不独立核验原始 PPI，不推导因果或机制真实性。';
    state.cy = window.cytoscape({ container: $('cy'), elements: model.elements, style: G.styles, layout: { name: 'grid', padding: 40 }, minZoom: 0.12, maxZoom: 3.5 });
    state.cy.on('tap', 'node, edge', event => { inspect(event.target); tab('details'); });
    state.cy.on('mouseover', 'node, edge', event => {
      const item = event.target; inspect(item);
      const position = event.renderedPosition || (item.isNode() ? item.renderedPosition() : item.renderedMidpoint());
      $('graph-tooltip').textContent = item.isNode() ? item.data('label') + '\n' + (GROUP_LABELS[item.data('type')] || item.data('type')) : EDGE_LABELS[item.data('category')] + '\n' + item.data('type');
      $('graph-tooltip').style.left = Math.max(8, Math.min(position.x + 16, $('cy').clientWidth - 220)) + 'px';
      $('graph-tooltip').style.top = Math.max(60, Math.min(position.y + 65, $('cy').clientHeight + 15)) + 'px';
      $('graph-tooltip').hidden = false;
    });
    state.cy.on('mouseout', 'node, edge', () => { $('graph-tooltip').hidden = true; });
    $('detail-title').textContent = '探索证据图'; $('detail-description').textContent = '悬停或点击节点与边，查看 JSON 中记录的类型、来源及属性。';
    $('detail-json').hidden = true; $('graph-tooltip').hidden = true;
    applyFilters(); renderPaths(); diagnostics(); layout();
    if (!model.nodeCount) tab('diagnostics');
  }
  function load(input, name) {
    const rows = G.records(input); G.build(rows[0]);
    state.rows = rows; state.source = name;
    $('result-select').replaceChildren();
    rows.forEach((row, index) => $('result-select').appendChild(new Option((index + 1) + '. ' + (row.drug1_name || row.drug1 || '药物 A') + ' + ' + (row.drug2_name || row.drug2 || '药物 B') + (row.error ? '（未完成）' : ''), String(index))));
    $('result-select').hidden = rows.length < 2; $('result-label').hidden = rows.length < 2;
    error(''); display(0);
  }
  async function openFile(file) {
    if (!file) return;
    try {
      if (file.size > G.LIMITS.fileBytes) throw new Error('文件超过 20 MiB。本查看器用于轻量解释结果，不读取原始数据集。');
      load(JSON.parse((await file.text()).replace(/^\uFEFF/, '')), file.name);
    } catch (exception) { error('文件未载入：' + exception.message + ' 当前视图仍为之前的记录。'); }
    $('file-input').value = '';
  }
  async function exportPng() {
    if (!state.cy || !state.model.nodeCount) return;
    const dataUrl = state.cy.png({ full: true, scale: 2, bg: '#ffffff', maxWidth: 2400, maxHeight: 1700 });
    const image = new Image();
    image.onload = () => {
      const canvas = document.createElement('canvas'); canvas.width = Math.max(image.width, 720); canvas.height = image.height + 76;
      const context = canvas.getContext('2d'); context.fillStyle = '#ffffff'; context.fillRect(0, 0, canvas.width, canvas.height); context.drawImage(image, (canvas.width - image.width) / 2, 0);
      context.fillStyle = '#233952'; context.font = '14px "Microsoft YaHei", sans-serif';
      context.fillText('InterpreterPLATO · ' + state.model.mode, 15, image.height + 27, canvas.width - 30);
      context.fillStyle = '#7d8ca0'; context.font = '11px "Microsoft YaHei", sans-serif';
      context.fillText('当前筛选视图 · 文件关系未独立核验 · PPI 无方向，连通不构成因果证据', 15, image.height + 51, canvas.width - 30);
      canvas.toBlob(blob => { if (blob) download(blob, 'InterpreterPLATO-' + (state.model.synthetic ? 'SYNTHETIC-' : '') + 'network.png'); else error('浏览器未能导出图像。'); }, 'image/png');
    };
    image.onerror = () => error('浏览器未能导出图像。'); image.src = dataUrl;
  }
  function boot() {
    if (!G || typeof window.cytoscape !== 'function') { error('查看器资源未加载完整。请保留目录内的资源文件，或使用打包后的单文件 HTML。'); return; }
    $('file-input').addEventListener('change', event => openFile(event.target.files[0]));
    $('load-example').addEventListener('click', () => load(window.InterpreterExample, '内置合成概念示例 · 无真实计算'));
    $('download-example').addEventListener('click', () => download(new Blob([JSON.stringify(window.InterpreterExample, null, 2)], { type: 'application/json' }), 'synthetic-explanation.json'));
    $('result-select').addEventListener('change', event => {
      const index = Number(event.target.value);
      try { error(''); display(index); } catch (exception) {
        const payload = state.rows[index];
        const empty = G.build({ ...payload, evidence_graphs: null });
        empty.payload = payload; empty.warnings = ['证据图未展示：' + exception.message];
        display(index, empty);
        error('该条记录无法展示：' + exception.message);
      }
    });
    document.querySelectorAll('[data-node-group], [data-edge-group]').forEach(input => input.addEventListener('change', applyFilters));
    $('layout-select').addEventListener('change', layout);
    $('fit-graph').addEventListener('click', () => state.cy && state.cy.fit(state.cy.elements().filter(item => !item.hasClass('filtered')), 45));
    $('clear-selection').addEventListener('click', () => { clearHighlight(); $('node-search').value = ''; $('search-status').textContent = '已取消高亮'; });
    $('node-search').addEventListener('input', () => {
      clearHighlight(); const query = $('node-search').value.trim().toLocaleLowerCase();
      if (!query || !state.cy) { $('search-status').textContent = '输入名称可高亮匹配节点'; return; }
      const matches = state.cy.nodes().filter(node => !node.hasClass('filtered') && (node.data('label') + ' ' + node.data('originalId')).toLocaleLowerCase().includes(query));
      state.cy.elements().addClass('dimmed'); matches.removeClass('dimmed').addClass('highlighted');
      $('search-status').textContent = '匹配 ' + matches.length + ' 个当前可见节点。';
    });
    $('node-search').addEventListener('keydown', event => {
      if (event.key !== 'Enter' || !state.cy) return;
      const matches = state.cy.nodes('.highlighted').filter(node => !node.hasClass('filtered'));
      if (matches.length) { inspect(matches[0]); tab('details'); state.cy.fit(matches, 80); }
    });
    $('export-png').addEventListener('click', () => exportPng().catch(exception => error(exception.message)));
    document.addEventListener('dragover', event => { event.preventDefault(); });
    document.addEventListener('drop', event => { event.preventDefault(); if (event.dataTransfer.files.length !== 1) error('请每次拖入一个解释 JSON 文件。'); else openFile(event.dataTransfer.files[0]); });
    document.querySelectorAll('label.button[for="file-input"]').forEach(label => { label.tabIndex = 0; label.setAttribute('role', 'button'); label.addEventListener('keydown', event => { if (event.key === 'Enter' || event.key === ' ') { event.preventDefault(); $('file-input').click(); } }); });
    $('tab-details').addEventListener('click', () => tab('details')); $('tab-diagnostics').addEventListener('click', () => tab('diagnostics'));
    const resize = () => { if (state.cy) { state.cy.resize(); } };
    if (typeof ResizeObserver === 'function') new ResizeObserver(resize).observe($('cy'));
    window.addEventListener('resize', resize);
    load(window.InterpreterExample, '内置合成概念示例 · 无真实计算');
  }
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', boot); else boot();
})();
