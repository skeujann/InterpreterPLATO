/* Optional distribution helper. The built HTML needs no Node or Python runtime. */
'use strict';
const fs = require('node:fs');
const path = require('node:path');
const crypto = require('node:crypto');

function build(output) {
  const root = __dirname;
  const dependency = JSON.parse(fs.readFileSync(path.join(root, 'vendor', 'dependency.json'), 'utf8'));
  for (const [file, expected] of Object.entries(dependency.files_sha256)) {
    const actual = crypto.createHash('sha256').update(fs.readFileSync(path.join(root, 'vendor', file))).digest('hex');
    if (actual !== expected) throw new Error('Vendored dependency integrity mismatch: ' + file);
  }
  let html = fs.readFileSync(path.join(root, 'index.html'), 'utf8');
  const css = fs.readFileSync(path.join(root, 'styles.css'), 'utf8');
  html = html.replace('<link rel="stylesheet" href="styles.css">', '<style>\n' + css + '\n</style>');
  for (const source of ['vendor/cytoscape.min.js', 'example.js', 'graph-model.js', 'app.js']) {
    const code = fs.readFileSync(path.join(root, source), 'utf8').replace(/<\/script/gi, '<\\/script');
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
  if (args.length && (args.length !== 2 || args[0] !== '--output')) throw new Error('Usage: node visualizer/build.cjs [--output PATH]');
  const output = args.length ? path.resolve(args[1]) : path.join(__dirname, 'dist', 'InterpreterPLATO-Visualizer.html');
  console.log(build(output));
}
module.exports = { build };
