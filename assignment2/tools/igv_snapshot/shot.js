const { chromium } = require('playwright');
(async () => {
  const b = await chromium.launch(); const p = await b.newPage({viewport: {width: 1400, height: 1000}});
  p.on('pageerror', e => console.log('err', e.message)); p.on('requestfailed', r => console.log('fail', r.url())); p.on('console', m => console.log('console:', m.text()));
  await p.goto('http://127.0.0.1:8765/index.html');
  await p.waitForTimeout(20000); console.log('ready', await p.evaluate(() => window.igvReady));
  await p.waitForTimeout(5000);
  await p.screenshot({path: process.argv[2], clip: {x: 0, y: 0, width: 1400, height: await p.evaluate(() => document.getElementById('igv').getBoundingClientRect().bottom + 10)}}); await b.close();
})();
