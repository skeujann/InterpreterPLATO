/* Optional distribution helper. The built HTML needs no Node or Python runtime. */
'use strict';
const fs = require('node:fs');
const path = require('node:path');
const crypto = require('node:crypto');

function build(output, options = {}) {
  const cases = options.cases || (options.caseFile ? [{ file: options.caseFile, pathIds: options.pathIds || [] }] : []);
  if (cases.length > 100) throw new Error('A local viewer supports at most 100 run examples');
  if (cases.some(item => path.resolve(output) === path.resolve(item.file))) throw new Error('The output must not overwrite the source run JSON');
  const root = __dirname;
  const dependency = JSON.parse(fs.readFileSync(path.join(root, 'vendor', 'dependency.json'), 'utf8'));
  for (const [file, expected] of Object.entries(dependency.files_sha256)) {
    const actual = crypto.createHash('sha256').update(fs.readFileSync(path.join(root, 'vendor', file))).digest('hex');
    if (actual !== expected) throw new Error('Vendored dependency integrity mismatch: ' + file);
  }
  let html = fs.readFileSync(path.join(root, 'index.html'), 'utf8');
  const localRuns = cases.map(item => require('./prepare-run.cjs').prepareRun(item.file, item.pathIds || []));
  const css = fs.readFileSync(path.join(root, 'styles.css'), 'utf8');
  html = html.replace('<link rel="stylesheet" href="styles.css">', '<style>\n' + css + '\n</style>');
  for (const source of ['vendor/cytoscape.min.js', 'example.js', 'public-demo.js', 'reference.js', 'graph-model.js', 'app.js']) {
    const embedded = source === 'app.js' && localRuns.length ? 'window.InterpreterRunExamples = ' + JSON.stringify(localRuns).replace(/</g, '\\u003c') + ';\nwindow.InterpreterRunExample = window.InterpreterRunExamples[0];\n' : '';
    const code = (embedded + fs.readFileSync(path.join(root, source), 'utf8')).replace(/<\/script/gi, '<\\/script');
    // A replacement callback preserves literal $&, $$ and similar sequences in JS.
    html = html.replace('<script src="' + source + '"></script>', () => '<script>\n' + code + '\n</script>');
  }
  const license = fs.readFileSync(path.join(root, 'vendor', 'CYTOSCAPE-LICENSE.txt'), 'utf8').replace(/--/g, '—');
  html = html.replace('</body>', '<!-- Cytoscape.js ' + dependency.version + ' (MIT)\n' + license + '\n-->\n</body>');
  fs.mkdirSync(path.dirname(output), { recursive: true });
  fs.writeFileSync(output, html);
  return output;
}
if (require.main === module) {
  const args = process.argv.slice(2);
  const options = { cases: [] }; let output = path.join(__dirname, 'dist', 'InterpreterPLATO-Visualizer.html');
  for (let i = 0; i < args.length; i += 2) {
    if (!args[i + 1]) throw new Error('Every build option needs a value');
    if (args[i] === '--output') output = path.resolve(args[i + 1]);
    else if (args[i] === '--case') options.cases.push({ file: path.resolve(args[i + 1]), pathIds: [] });
    else if (args[i] === '--focus-path') {
      if (!options.cases.length) throw new Error('--focus-path requires a preceding --case');
      options.cases[options.cases.length - 1].pathIds.push(args[i + 1]);
    }
    else throw new Error('Usage: node visualizer/build.cjs [--output PATH] [--case PATH] [--focus-path PATH_ID]');
  }
  console.log(build(output, options));
}
module.exports = { build };
