/* Generate the offline browser module from the frozen, downloadable JSON. */
'use strict';
const fs = require('node:fs');
const path = require('node:path');
const crypto = require('node:crypto');
const input = path.join(__dirname, 'examples/dabrafenib-trametinib.json');
const bytes = fs.readFileSync(input);
const expected = fs.readFileSync(input + '.sha256', 'utf8').trim().split(/\s+/)[0];
if (crypto.createHash('sha256').update(bytes).digest('hex') !== expected) throw new Error('Reference JSON hash mismatch');
const data = JSON.parse(bytes.toString('utf8'));
const source = '/* Public reference annotations. See REFERENCE_CASE.md for attribution and scope. */\n' +
  '(function (root, factory) {\n' +
  "  if (typeof module === 'object' && module.exports) module.exports = factory();\n" +
  '  else root.InterpreterReference = factory();\n' +
  "})(typeof globalThis !== 'undefined' ? globalThis : this, function () {\n" +
  "  'use strict';\n  return " + JSON.stringify(data, null, 2) + ';\n});\n';
fs.writeFileSync(path.join(__dirname, 'reference.js'), source);
console.log('Generated reference.js from verified reference JSON');
