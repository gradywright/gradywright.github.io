/* Grady B. Wright — site behaviour. No dependencies.
   1. Mobile nav toggle   2. Header shadow on scroll
   3. Scroll reveal       4. Hero node-sphere (canvas) */
(function () {
  'use strict';
  var reduce = window.matchMedia && matchMedia('(prefers-reduced-motion: reduce)').matches;

  /* 1. Mobile nav */
  var toggle = document.querySelector('[data-nav-toggle]');
  var nav = document.querySelector('[data-nav]');
  if (toggle && nav) {
    toggle.addEventListener('click', function () {
      var open = toggle.getAttribute('aria-expanded') !== 'true';
      toggle.setAttribute('aria-expanded', String(open));
      nav.classList.toggle('is-open', open);
    });
    document.addEventListener('keydown', function (e) {
      if (e.key === 'Escape' && toggle.getAttribute('aria-expanded') === 'true') {
        toggle.setAttribute('aria-expanded', 'false');
        nav.classList.remove('is-open');
        toggle.focus();
      }
    });
  }

  /* 2. Header shadow */
  var header = document.querySelector('[data-header]');
  if (header) {
    var onScroll = function () { header.classList.toggle('is-scrolled', window.scrollY > 12); };
    onScroll();
    window.addEventListener('scroll', onScroll, { passive: true });
  }

  /* 3. Scroll reveal */
  var items = document.querySelectorAll('[data-reveal]');
  if (items.length) {
    if (!('IntersectionObserver' in window) || reduce) {
      items.forEach(function (el) { el.classList.add('is-in'); });
    } else {
      var io = new IntersectionObserver(function (entries) {
        entries.forEach(function (en) {
          if (en.isIntersecting) { en.target.classList.add('is-in'); io.unobserve(en.target); }
        });
      }, { rootMargin: '0px 0px -8% 0px', threshold: 0.08 });
      items.forEach(function (el, i) {
        if (!el.style.getPropertyValue('--d')) el.style.setProperty('--d', (i % 6) * 70 + 'ms');
        io.observe(el);
      });
    }
  }

  /* 3b. Research page: mark the year currently in view in the year nav */
  var years = document.querySelectorAll('.year[id]');
  var yearLinks = document.querySelectorAll('.year-nav a');
  if (years.length && yearLinks.length && 'IntersectionObserver' in window) {
    var mark = function (id) {
      yearLinks.forEach(function (a) {
        var on = a.getAttribute('href') === '#' + id;
        if (on) { a.setAttribute('aria-current', 'true'); } else { a.removeAttribute('aria-current'); }
        if (on) { var bar = a.parentNode.parentNode; bar.scrollLeft = a.offsetLeft - bar.clientWidth / 2; }
      });
    };
    var yearObserver = new IntersectionObserver(function (entries) {
      entries.forEach(function (en) { if (en.isIntersecting) mark(en.target.id); });
    }, { rootMargin: '-25% 0px -65% 0px' });
    years.forEach(function (y) { yearObserver.observe(y); });
  }

  /* 4. Node sphere: scattered "mesh-free" nodes on S^2, slowly rotating. */
  var canvas = document.querySelector('canvas[data-sphere]');
  if (!canvas || !canvas.getContext) return;
  var ctx = canvas.getContext('2d');
  var N = 560, golden = Math.PI * (3 - Math.sqrt(5)), pts = [];
  for (var i = 0; i < N; i++) {
    var y = 1 - 2 * (i + 0.5) / N, r = Math.sqrt(1 - y * y), t = golden * i;
    pts.push({ x: Math.cos(t) * r, y: y, z: Math.sin(t) * r, hot: i % 41 === 0 });
  }
  var w = 0, h = 0, angle = 0.8, tilt = -0.42, tiltTarget = tilt, last = 0, raf = 0, visible = true;

  function resize() {
    var rect = canvas.getBoundingClientRect(), dpr = Math.min(window.devicePixelRatio || 1, 2);
    w = rect.width; h = rect.height;
    canvas.width = Math.round(w * dpr); canvas.height = Math.round(h * dpr);
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    draw();
  }

  function draw() {
    ctx.clearRect(0, 0, w, h);
    var R = Math.min(w, h) * 0.43, cx = w / 2, cy = h / 2;
    ctx.globalAlpha = 1;
    ctx.strokeStyle = 'rgba(180, 193, 214, 0.14)'; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.arc(cx, cy, R, 0, 6.2832); ctx.stroke();

    var ca = Math.cos(angle), sa = Math.sin(angle), ct = Math.cos(tilt), st = Math.sin(tilt);
    for (var k = 0; k < N; k++) {
      var p = pts[k];
      var x = p.x * ca + p.z * sa, z = -p.x * sa + p.z * ca;      // spin about the vertical axis
      var y2 = p.y * ct - z * st, z2 = p.y * st + z * ct;          // tip toward the viewer
      var d = (z2 + 1) / 2;                                        // 0 = far side, 1 = near side
      ctx.globalAlpha = 0.1 + 0.9 * Math.pow(d, 1.7);
      ctx.fillStyle = p.hot ? '#F4A168' : '#BED0EB';
      ctx.beginPath();
      ctx.arc(cx + x * R, cy + y2 * R, (0.7 + 1.7 * d) * (p.hot ? 1.45 : 1), 0, 6.2832);
      ctx.fill();
    }
    ctx.globalAlpha = 1;
  }

  function frame(now) {
    raf = 0;
    if (!visible) return;
    var dt = last ? Math.min(now - last, 50) : 16; last = now;
    angle += dt * 0.00009;
    tilt += (tiltTarget - tilt) * 0.05;
    draw();
    raf = requestAnimationFrame(frame);
  }
  function start() { if (!raf && !reduce) { last = 0; raf = requestAnimationFrame(frame); } }

  resize();
  window.addEventListener('resize', resize);
  if (reduce) return;                                              // static frame only

  window.addEventListener('pointermove', function (e) {
    var rect = canvas.getBoundingClientRect();
    var ny = (e.clientY - (rect.top + rect.height / 2)) / Math.max(window.innerHeight, 1);
    tiltTarget = -0.42 + Math.max(-0.5, Math.min(0.5, ny)) * 0.25;
  }, { passive: true });
  if ('IntersectionObserver' in window) {
    new IntersectionObserver(function (en) { visible = en[0].isIntersecting; if (visible) start(); }).observe(canvas);
  }
  start();
})();
