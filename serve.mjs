// serve.mjs — static file server for local work.
//
// Sends Cache-Control: no-store on everything. The default Python http.server
// lets the browser keep an ES module across a reload, so an edited shader can
// look like it changed nothing at all; several minutes of a debugging session
// went into that once and it is not going in again.
//
//   node serve.mjs [port]

import { createServer } from 'node:http';
import { readFile, stat, writeFile, mkdir } from 'node:fs/promises';
import { extname, join, normalize, resolve, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';

const port = Number(process.argv[2] ?? 8811);
// fileURLToPath, not url.pathname: on Windows the latter keeps the leading
// slash and the percent-encoding, and every request 404s.
const root = resolve(dirname(fileURLToPath(import.meta.url)));

const TYPES = {
  '.html': 'text/html; charset=utf-8',
  '.js': 'text/javascript; charset=utf-8',
  '.mjs': 'text/javascript; charset=utf-8',
  '.css': 'text/css; charset=utf-8',
  '.json': 'application/json; charset=utf-8',
  '.md': 'text/plain; charset=utf-8',
  '.svg': 'image/svg+xml',
  '.png': 'image/png',
  '.ico': 'image/x-icon',
};

createServer(async (req, res) => {
  try {
    // POST /__capture — save a still of the canvas so a render can be LOOKED at.
    // Development only; nothing the page ships uses it. Writes go to shots/,
    // which is gitignored.
    if (req.method === 'POST' && req.url.startsWith('/__capture')) {
      const name = (new URL(req.url, 'http://x').searchParams.get('name') || 'shot')
        .replace(/[^A-Za-z0-9_-]/g, '');
      const chunks = [];
      for await (const c of req) chunks.push(c);
      const body = Buffer.concat(chunks).toString();
      const b64 = body.includes(',') ? body.split(',').pop() : body;
      await mkdir(join(root, 'shots'), { recursive: true });
      const out = join(root, 'shots', name + '.jpg');
      await writeFile(out, Buffer.from(b64, 'base64'));
      res.writeHead(200, { 'Content-Type': 'text/plain' }).end(out);
      return;
    }
    const url = new URL(req.url, 'http://x');
    let p = decodeURIComponent(url.pathname);
    if (p.endsWith('/')) p += 'index.html';
    const file = join(root, normalize(p).replace(/^(\.\.[/\\])+/, ''));
    if (!file.startsWith(root)) { res.writeHead(403).end('forbidden'); return; }
    const s = await stat(file);
    if (s.isDirectory()) { res.writeHead(302, { Location: p + '/' }).end(); return; }
    const body = await readFile(file);
    res.writeHead(200, {
      'Content-Type': TYPES[extname(file).toLowerCase()] ?? 'application/octet-stream',
      'Cache-Control': 'no-store, must-revalidate',
      'Content-Length': body.length,
    });
    res.end(body);
  } catch (e) {
    res.writeHead(e.code === 'ENOENT' ? 404 : 500).end(String(e.code ?? e));
  }
}).listen(port, () => console.log(`lava-lamp-lab on http://localhost:${port}`));
