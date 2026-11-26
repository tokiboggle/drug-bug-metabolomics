/**
 * 27200 Data-driven Bioengineering — Dataset ME2 Oral Exam Slide Deck Builder (QA-Polished v3)
 * Generates: slides/27200_ME2_Omeprazole_Metabolomics.pptx (15 slides, LAYOUT_WIDE 13.333" x 7.5")
 *
 * Run with:
 *   NODE_PATH=/tmp/pptx-tools/node_modules node src/build_slides.js
 */

const path = require("path");
const pptxgen = require("pptxgenjs");
const React = require("react");
const ReactDOMServer = require("react-dom/server");
const sharp = require("sharp");
const {
  FaFlask,
  FaMicroscope,
  FaExclamationTriangle,
  FaRobot,
  FaShieldAlt,
  FaClock,
  FaFilter,
  FaBalanceScale,
  FaSearch,
  FaLayerGroup,
  FaArrowRight,
  FaBug,
  FaBullseye,
} = require("react-icons/fa");

const ROOT = path.resolve(__dirname, "..");
const FIG_DIR = path.join(ROOT, "results", "figures");
const OUT_PATH = path.join(ROOT, "slides", "27200_ME2_Omeprazole_Metabolomics.pptx");

// ---------------------------------------------------------------------------
// Palette — Bio-Slate + Analytical Teal + Omeprazole Vermilion
// ---------------------------------------------------------------------------
const C = {
  darkBg: "0F2430",       // Title & closing slide deep slate-navy
  darkCard: "173445",     // Subtle dark card surface
  lightBg: "F7F9FB",      // Content slide off-white
  white: "FFFFFF",
  primary: "143A4A",      // Dominant deep analytical teal-slate
  omep: "C85200",         // Sharp vermilion accent (Omeprazole)
  dmso: "0072B2",         // Deep blue (DMSO control)
  green: "008B68",        // Emerald (Medium / Confirmed bio)
  violet: "6A3D9A",       // Violet (ESI suppression / Regioisomers)
  crimson: "B2182B",      // Crimson (QC failure / False hit)
  amber: "9A6108",        // Rich high-contrast amber
  text: "1E293B",         // Primary body charcoal
  textSoft: "334155",     // High-contrast secondary text
  muted: "52606D",        // Captions & footers
  border: "D5DFE8",       // Subtle card border
  tintTeal: "EAF2F5",     // Soft teal tint
  lightOnDark: "EAF2F7",  // High-contrast light text on dark bg
  mutedOnDark: "BFD3DF",  // High-contrast secondary text on dark bg
};

const FONT_HEAD = "Georgia";
const FONT_BODY = "Calibri";

// Fresh shadow factory (NEVER reuse PptxGenJS option objects!)
const makeShadow = (opacity = 0.09) => ({
  type: "outer",
  color: "0F2430",
  blur: 6,
  offset: 2,
  angle: 135,
  opacity,
});

// Parse markup _{sub} and ^{sup} into native PptxGenJS rich-text runs
function chemRuns(str, baseOpts = {}) {
  const runs = [];
  const regex = /(_\{[^}]+\}|\^\{[^}]+\})/g;
  let last = 0;
  let m;
  while ((m = regex.exec(str)) !== null) {
    if (m.index > last) {
      runs.push({ text: str.slice(last, m.index), options: { ...baseOpts } });
    }
    const tok = m[0];
    const inner = tok.slice(2, -1);
    if (tok.startsWith("_")) {
      runs.push({ text: inner, options: { ...baseOpts, subscript: true } });
    } else {
      runs.push({ text: inner, options: { ...baseOpts, superscript: true } });
    }
    last = regex.lastIndex;
  }
  if (last < str.length) {
    runs.push({ text: str.slice(last), options: { ...baseOpts } });
  }
  return runs;
}

// Build clean multi-paragraph bullet runs with paraSpaceAfter on FIRST run of each paragraph!
function makeChemBullets(items, fontSize = 10.2, paraGap = 8, textColor = C.text) {
  const out = [];
  items.forEach((txt, idx) => {
    const runs = chemRuns(txt, { fontFace: FONT_BODY, fontSize, color: textColor });
    runs[0].options.bullet = true;
    if (idx < items.length - 1) {
      runs[0].options.paraSpaceAfter = paraGap;
      runs[runs.length - 1].options.breakLine = true;
    }
    out.push(...runs);
  });
  return out;
}

// Rasterize react-icons SVG to high-DPI base64 PNG
async function iconPng(IconComp, hexColor, size = 256) {
  const svg = ReactDOMServer.renderToStaticMarkup(
    React.createElement(IconComp, { color: `#${hexColor}`, size: String(size) })
  );
  const buf = await sharp(Buffer.from(svg)).png().toBuffer();
  return "image/png;base64," + buf.toString("base64");
}

// Exact aspect-ratio-preserving image placement
async function fitImageBox(imgPath, bx, by, bw, bh) {
  const meta = await sharp(imgPath).metadata();
  const ratio = meta.width / meta.height;
  let w = bw;
  let h = w / ratio;
  if (h > bh) {
    h = bh;
    w = h * ratio;
  }
  const x = bx + (bw - w) / 2;
  const y = by + (bh - h) / 2;
  return { x, y, w, h };
}

// Standard content-slide chrome: >=0.5" L/R margins, single-line title, NO underline accent
function addChrome(pres, slide, { section, title, speakerTag, speakerColor, slideNum, totalSlides = 15 }) {
  slide.background = { color: C.lightBg };

  slide.addText(section.toUpperCase(), {
    x: 0.55,
    y: 0.30,
    w: 9.6,
    h: 0.22,
    fontFace: FONT_BODY,
    fontSize: 9.0,
    bold: true,
    color: C.omep,
    charSpacing: 1.4,
    margin: 0,
  });

  slide.addText(title, {
    x: 0.55,
    y: 0.52,
    w: 9.75,
    h: 0.44,
    fontFace: FONT_HEAD,
    fontSize: 20.5,
    bold: true,
    color: C.primary,
    margin: 0,
    valign: "middle",
  });

  if (speakerTag) {
    slide.addShape(pres.shapes.ROUNDED_RECTANGLE, {
      x: 10.50,
      y: 0.40,
      w: 2.28,
      h: 0.38,
      rectRadius: 0.19,
      fill: { color: speakerColor || C.primary },
      line: { color: speakerColor || C.primary, width: 0 },
    });
    slide.addText(speakerTag, {
      x: 10.50,
      y: 0.40,
      w: 2.28,
      h: 0.38,
      fontFace: FONT_BODY,
      fontSize: 9.2,
      bold: true,
      color: C.white,
      align: "center",
      valign: "middle",
      margin: 0,
    });
  }

  slide.addText("27200 Data-Driven Bioengineering  |  Dataset ME2: Omeprazole & Com20 Gut Community", {
    x: 0.55,
    y: 7.08,
    w: 9.5,
    h: 0.22,
    fontFace: FONT_BODY,
    fontSize: 8.5,
    color: C.muted,
    margin: 0,
  });
  slide.addText(`${slideNum} / ${totalSlides}`, {
    x: 11.78,
    y: 7.08,
    w: 1.0,
    h: 0.22,
    fontFace: FONT_BODY,
    fontSize: 8.5,
    bold: true,
    color: C.primary,
    align: "right",
    margin: 0,
  });
}

// Clean white card with thick single-side left accent bar (RECTANGLE)
function addAccentCard(pres, slide, { x, y, w, h, accentColor, fillColor = C.white, title, body, titleSize = 11.2, bodySize = 9.9, lineSpacing = 1.18 }) {
  slide.addShape(pres.shapes.RECTANGLE, {
    x,
    y,
    w,
    h,
    fill: { color: fillColor },
    line: { color: C.border, width: 0.75 },
    shadow: makeShadow(0.08),
  });
  slide.addShape(pres.shapes.RECTANGLE, {
    x,
    y,
    w: 0.07,
    h,
    fill: { color: accentColor },
    line: { color: accentColor, width: 0 },
  });
  const padL = 0.18;
  const padR = 0.12;
  const titleRuns = typeof title === "string" ? chemRuns(title, { fontFace: FONT_BODY, fontSize: titleSize, bold: true, color: accentColor }) : title;
  slide.addText(titleRuns, {
    x: x + padL,
    y: y + 0.08,
    w: w - padL - padR,
    h: 0.27,
    margin: 0,
    valign: "top",
  });
  const bodyRuns = typeof body === "string" ? chemRuns(body, { fontFace: FONT_BODY, fontSize: bodySize, color: C.text }) : body;
  slide.addText(bodyRuns, {
    x: x + padL,
    y: y + 0.36,
    w: w - padL - padR,
    h: h - 0.42,
    fontFace: FONT_BODY,
    fontSize: bodySize,
    color: C.text,
    margin: 0,
    valign: "top",
    lineSpacingMultiple: lineSpacing,
  });
}

// Tight-framed figure card (no awkward empty white letterboxing!)
async function addFigureCard(pres, slide, imgFile, { x, y, w, h, tightFrame = true }) {
  const imgPath = path.join(FIG_DIR, imgFile);
  const pad = 0.07;
  const box = await fitImageBox(imgPath, x + pad, y + pad, w - 2 * pad, h - 2 * pad);
  const fx = tightFrame ? box.x - pad : x;
  const fy = tightFrame ? box.y - pad : y;
  const fw = tightFrame ? box.w + 2 * pad : w;
  const fh = tightFrame ? box.h + 2 * pad : h;

  slide.addShape(pres.shapes.RECTANGLE, {
    x: fx,
    y: fy,
    w: fw,
    h: fh,
    fill: { color: C.white },
    line: { color: C.border, width: 0.8 },
    shadow: makeShadow(0.09),
  });
  slide.addImage({ path: imgPath, ...box });
  return { fx, fy, fw, fh };
}

async function buildDeck() {
  const pres = new pptxgen();
  pres.layout = "LAYOUT_WIDE"; // 13.333" x 7.5"
  pres.author = "Tóki, Magda, Caroline & Irene — Group ME2 (27200)";
  pres.title = "Omeprazole and a Bacterial Community Over Time — Untargeted LC-MS/MS Metabolomics";

  const ic = {
    flask: await iconPng(FaFlask, C.white),
    micro: await iconPng(FaMicroscope, C.white),
    warn: await iconPng(FaExclamationTriangle, C.white),
    robot: await iconPng(FaRobot, C.white),
    shield: await iconPng(FaShieldAlt, C.white),
    clock: await iconPng(FaClock, C.white),
    filter: await iconPng(FaFilter, C.white),
    scale: await iconPng(FaBalanceScale, C.white),
    search: await iconPng(FaSearch, C.white),
    layers: await iconPng(FaLayerGroup, C.white),
    bug: await iconPng(FaBug, C.white),
    target: await iconPng(FaBullseye, C.white),
    arrowTeal: await iconPng(FaArrowRight, C.primary),
  };

  // =========================================================================
  // SLIDE 1 — DARK TITLE SLIDE (Zero decorative collisions + crisp chips)
  // =========================================================================
  {
    const s = pres.addSlide();
    s.background = { color: C.darkBg };

    // Clean geometric bokeh motif shifted strictly into top-right/far-right margin (zero overlap!)
    s.addShape(pres.shapes.OVAL, { x: 10.35, y: -0.95, w: 4.1, h: 4.1, fill: { color: "17384A", transparency: 30 }, line: { color: "265973", width: 1 } });
    s.addShape(pres.shapes.OVAL, { x: 11.15, y: 2.35, w: 2.9, h: 2.9, fill: { color: "C85200", transparency: 80 }, line: { color: "C85200", width: 1 } });
    s.addShape(pres.shapes.OVAL, { x: 9.75, y: 3.75, w: 1.55, h: 1.55, fill: { color: "008B68", transparency: 76 }, line: { color: "008B68", width: 1 } });

    // Left vertical vermilion motif bar
    s.addShape(pres.shapes.RECTANGLE, { x: 0.72, y: 1.15, w: 0.09, h: 4.30, fill: { color: C.omep }, line: { color: C.omep, width: 0 } });

    // Course pill
    s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: 1.08, y: 1.15, w: 4.65, h: 0.40, rectRadius: 0.20, fill: { color: "1A3D52" }, line: { color: "387190", width: 1 } });
    s.addText("27200 DATA-DRIVEN BIOENGINEERING · DATASET ME2", {
      x: 1.08, y: 1.15, w: 4.65, h: 0.40, fontFace: FONT_BODY, fontSize: 9.6, bold: true, color: "BCE4F7", align: "center", valign: "middle", charSpacing: 0.8, margin: 0,
    });

    // Main Title
    s.addText("Omeprazole & a Synthetic Gut\nBacterial Community Over Time", {
      x: 1.08, y: 1.78, w: 8.8, h: 1.62, fontFace: FONT_HEAD, fontSize: 38, bold: true, color: C.white, lineSpacingMultiple: 1.12, margin: 0,
    });

    // Clean balanced 2-line subtitle
    s.addText("Untargeted LC-MS/MS metabolomics of Com20 (omeprazole vs. DMSO, 0–8 h):\nseparating active microbial drug biotransformation from LC-MS matrix artifacts", {
      x: 1.08, y: 3.56, w: 8.4, h: 0.85, fontFace: FONT_BODY, fontSize: 15, color: C.lightOnDark, lineSpacingMultiple: 1.28, margin: 0,
    });

    // Compact Group/Source metadata bar — cleanly clears decorative circles by >1.4"!
    s.addShape(pres.shapes.RECTANGLE, { x: 1.08, y: 4.78, w: 7.45, h: 0.64, fill: { color: C.darkCard }, line: { color: "30617B", width: 0.9 } });
    s.addText([
      { text: "Group: ", options: { bold: true, color: "F4A261" } },
      { text: "Tóki · Magda · Caroline · Irene   |   ", options: { color: C.white } },
      { text: "Source: ", options: { bold: true, color: "F4A261" } },
      { text: "Pakkir Shah et al. (2026), bioRxiv", options: { color: C.lightOnDark } },
    ], {
      x: 1.24, y: 4.78, w: 7.15, h: 0.64, fontFace: FONT_BODY, fontSize: 11.5, valign: "middle", margin: 0,
    });

    // High-contrast Course-Toolkit chips with generous spacing
    const weeks = ["W2 · Structures", "W3 · Networks", "W4 · Design & Stats", "W6–7 · Enrichment", "W8 · LC-MS/MS", "W9 · ML & PCA", "W13 · Graphs"];
    const chipW = 1.62;
    const chipGap = 0.14;
    const startX = 0.72;
    weeks.forEach((wk, i) => {
      const wx = startX + i * (chipW + chipGap);
      s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: wx, y: 6.18, w: chipW, h: 0.40, rectRadius: 0.18, fill: { color: "1C4359" }, line: { color: "4A86A6", width: 1.0 } });
      s.addText(wk, { x: wx, y: 6.18, w: chipW, h: 0.40, fontFace: FONT_BODY, fontSize: 9.2, bold: true, color: C.white, align: "center", valign: "middle", margin: 0 });
    });

    s.addNotes(
      "Welcome examiners and peer group. Today we present our analysis of Dataset ME2 from Pakkir Shah et al. (2026): untargeted LC-MS/MS metabolomics of the Com20 synthetic gut bacterial community incubated with omeprazole versus DMSO over 0 to 8 hours. Our analysis integrates six course modules—from experimental design and LC-MS QC to ML, chemical-class enrichment, multi-layer networks, and structural zoom-in."
    );
  }

  // =========================================================================
  // SLIDE 2 — BIOLOGICAL QUESTION & EXECUTIVE ANSWER
  // =========================================================================
  {
    const s = pres.addSlide();
    addChrome(pres, s, {
      section: "01 · Biological Question & Main Answer",
      title: "How Does Com20 Respond to Omeprazole Over 0–8 Hours?",
      speakerTag: "OVERVIEW · ALL",
      speakerColor: C.primary,
      slideNum: 2,
    });

    // Hero question banner
    s.addShape(pres.shapes.RECTANGLE, { x: 0.55, y: 1.18, w: 12.23, h: 0.72, fill: { color: C.primary }, shadow: makeShadow(0.10) });
    s.addImage({ data: ic.target, x: 0.75, y: 1.33, w: 0.42, h: 0.42 });
    s.addText([
      { text: "Research Question:  ", options: { bold: true, color: "F4A261", fontSize: 13 } },
      { text: "How do omeprazole-treated Com20 bacterial cultures change over time (T0–T8) compared with DMSO controls?", options: { color: C.white, fontSize: 13 } },
    ], { x: 1.32, y: 1.18, w: 11.2, h: 0.72, fontFace: FONT_BODY, valign: "middle", margin: 0 });

    // Three stat callouts — darker crisp descriptions, balanced height, zero split formulas!
    const statsCards = [
      {
        num: "50.8%",
        sub: "Global Variance from Treatment",
        desc: "PERMANOVA across 3,414 filtered features (p = 0.001) — driven primarily by microbial drug transformation.",
        col: C.omep,
      },
      {
        num: "59.1%",
        sub: "Clean Endogenous Variance from Time",
        desc: "Outside the RT 1.80–2.25 min suppression zone (n = 1,915), Com20 preserves core 0–8 h growth in both arms.",
        col: C.dmso,
      },
      {
        num: "9× Surge",
        sub: "Ufiprazole (Omeprazole Sulfide)",
        desc: "Anaerobic bacterial sulfoxide reduction (ID 2962, 330.13 Da) reaches ~47% of late TIC and feeds a cascade.",
        col: C.green,
      },
    ];
    statsCards.forEach((sc, idx) => {
      const cx = 0.55 + idx * 4.19;
      s.addShape(pres.shapes.RECTANGLE, { x: cx, y: 2.20, w: 3.85, h: 1.98, fill: { color: C.white }, line: { color: C.border, width: 0.8 }, shadow: makeShadow(0.08) });
      s.addShape(pres.shapes.RECTANGLE, { x: cx, y: 2.20, w: 3.85, h: 0.08, fill: { color: sc.col }, line: { color: sc.col, width: 0 } });
      s.addText(sc.num, { x: cx + 0.20, y: 2.38, w: 3.45, h: 0.62, fontFace: FONT_HEAD, fontSize: 34, bold: true, color: sc.col, margin: 0 });
      s.addText(sc.sub, { x: cx + 0.20, y: 3.04, w: 3.45, h: 0.28, fontFace: FONT_BODY, fontSize: 12, bold: true, color: C.primary, margin: 0 });
      s.addText(sc.desc, { x: cx + 0.20, y: 3.35, w: 3.45, h: 0.72, fontFace: FONT_BODY, fontSize: 10.6, color: C.textSoft, lineSpacingMultiple: 1.18, margin: 0 });
    });

    // Two bottom cards with generous 0.32" top gap, clean paragraph spacing & balanced fill!
    addAccentCard(pres, s, {
      x: 0.55, y: 4.50, w: 5.95, h: 2.28, accentColor: C.omep,
      title: "Biological Story: Two Intertwined Responses",
      titleSize: 13.0,
      bodySize: 11.2,
      lineSpacing: 1.20,
      body: [
        { text: "1. What bacteria do to omeprazole: ", options: { bold: true, color: C.omep, paraSpaceAfter: 10 } },
        { text: "Anaerobic sulfoxide reduction to ufiprazole, O-demethylation (ID 2465), C–S cleavage (5-MeO-2-MBI, ID 2480), regioselective methylation (ID 3605 / ID 2500), and 2S conjugates (ID 2368, +95×).", options: { color: C.text, breakLine: true } },
        { text: "2. What omeprazole does to bacteria: ", options: { bold: true, color: C.dmso } },
        { text: "Core growth metabolism (bile acid deconjugation, N-acyl amino acids, lipid uptake) stays conserved, with selective shifts (2× attenuated indole-3-propionic acid, delayed PE 13:0).", options: { color: C.text } },
      ],
    });

    addAccentCard(pres, s, {
      x: 6.83, y: 4.50, w: 5.95, h: 2.28, accentColor: C.violet,
      title: "Why Rigorous QC Was Essential Before Biology",
      titleSize: 13.0,
      bodySize: 11.0,
      lineSpacing: 1.18,
      body: [
        { text: "Four hidden analytical traps would invert conclusions if missed:", options: { bold: true, color: C.primary, paraSpaceAfter: 8, breakLine: true } },
        ...makeChemBullets([
          "4 failed blank-like T4/T5 injections (TIC < solvent blanks)",
          "5.3-day LC-MS queue decay of acid-labile parent omeprazole (T1/T2 = 0)",
          "TIC normalization artifact (83% fake down-regulated endogenous hits)",
          "Co-elution ESI ion suppression under ufiprazole at RT 1.80–2.25 min (OR = 46.5)",
        ], 10.8, 6),
      ],
    });

    s.addNotes(
      "Before diving into details, here is our central takeaway: omeprazole does NOT collapse Com20 bacterial metabolism over 8 hours—59.1% of clean endogenous variance is conserved biological growth time. Instead, the 50.8% global treatment variance comes from active multi-step microbial biotransformation of omeprazole itself, led by a 9-fold surge in ufiprazole. Equally importantly, four major LC-MS artifacts had to be caught and corrected so we wouldn't mistake matrix effects for biology."
    );
  }

  // =========================================================================
  // SLIDE 3 — EXPERIMENTAL DESIGN & CONTROL ARCHITECTURE (Tóki)
  // =========================================================================
  {
    const s = pres.addSlide();
    addChrome(pres, s, {
      section: "02 · Experimental Design & Data Audit (Week 4 & Week 8)",
      title: "Com20 Longitudinal Design & Control Architecture",
      speakerTag: "TÓKI · DESIGN & QC",
      speakerColor: C.dmso,
      slideNum: 3,
    });

    // Top 4-step process flow — ALL 4 headers fit cleanly on 1 single line!
    const steps = [
      { icon: ic.flask, col: C.primary, h: "Com20 Gut SynCom", t: "20 commensal gut strains\n(6 phyla, ~61% gut pathways)\nAnaerobic mGAM, 37 °C" },
      { icon: ic.clock, col: C.omep, h: "2 Treatments × 9 Hours", t: "Omeprazole (8 μM) vs DMSO\nHourly T0–T8 (0–8 h)\n3 biological replicates (n = 54)" },
      { icon: ic.micro, col: C.dmso, h: "Untargeted LC-MS/MS", t: "EtOAc extract → 50% MeOH\nC18 RP (5 min), ESI+ DDA\nQ Exactive HF Orbitrap" },
      { icon: ic.filter, col: C.green, h: "MZmine / IIMN / GNPS", t: "4,525 features × 69 injections\n335 annotated (7.4% raw)\n12-drug study subset table" },
    ];
    steps.forEach((st, i) => {
      const sx = 0.55 + i * 3.16;
      s.addShape(pres.shapes.RECTANGLE, { x: sx, y: 1.18, w: 2.76, h: 1.62, fill: { color: C.white }, line: { color: C.border, width: 0.8 }, shadow: makeShadow(0.08) });
      s.addShape(pres.shapes.OVAL, { x: sx + 0.14, y: 1.32, w: 0.40, h: 0.40, fill: { color: st.col }, line: { color: st.col, width: 0 } });
      s.addImage({ data: st.icon, x: sx + 0.22, y: 1.40, w: 0.24, h: 0.24 });
      s.addText(st.h, { x: sx + 0.62, y: 1.32, w: 2.06, h: 0.40, fontFace: FONT_BODY, fontSize: 11, bold: true, color: st.col, valign: "middle", margin: 0 });
      s.addText(st.t, { x: sx + 0.16, y: 1.84, w: 2.46, h: 0.88, fontFace: FONT_BODY, fontSize: 9.8, color: C.text, lineSpacingMultiple: 1.22, margin: 0 });
      if (i < 3) {
        s.addImage({ data: ic.arrowTeal, x: sx + 2.86, y: 1.89, w: 0.20, h: 0.20 });
      }
    });

    // Bottom section heading aligned cleanly across BOTH columns
    s.addText("Control Architecture & Filtering Roles (69 Injections)", {
      x: 0.55, y: 3.12, w: 8.15, h: 0.28, fontFace: FONT_BODY, fontSize: 12.2, bold: true, color: C.primary, margin: 0,
    });
    s.addText("Key Dataset Audit Findings", {
      x: 9.00, y: 3.12, w: 3.78, h: 0.28, fontFace: FONT_BODY, fontSize: 12.2, bold: true, color: C.primary, margin: 0,
    });

    const headerOpts = { fill: { color: C.primary }, color: C.white, bold: true, fontSize: 10, fontFace: FONT_BODY, valign: "middle" };
    const cellOpts = (bg) => ({ fill: { color: bg }, color: C.text, fontSize: 9.6, fontFace: FONT_BODY, valign: "middle" });
    const rows = [
      [
        { text: "Sample / Control Group", options: headerOpts },
        { text: "n", options: { ...headerOpts, align: "center" } },
        { text: "What It Contains", options: headerOpts },
        { text: "Analytical Role in Our Pipeline", options: headerOpts },
      ],
      [
        { text: "Bacterial (Omepr / DMSO)", options: { ...cellOpts(C.white), bold: true, color: C.omep } },
        { text: "54 → 50", options: { ...cellOpts(C.white), align: "center", bold: true } },
        { text: "Com20 + drug/vehicle (T0–T8 × 3 reps)", options: cellOpts(C.white) },
        { text: "Primary modeling cohort (4 failed T4/T5 removed)", options: cellOpts(C.white) },
      ],
      [
        { text: "Medium Control (T0 only)", options: { ...cellOpts(C.tintTeal), bold: true, color: C.green } },
        { text: "3", options: { ...cellOpts(C.tintTeal), align: "center", bold: true } },
        { text: "Sterile mGAM medium, no bacteria, no drug", options: cellOpts(C.tintTeal) },
        { text: "Separates mGAM background (DKPs) from bacterial ions", options: cellOpts(C.tintTeal) },
      ],
      [
        { text: "Solvent Blanks", options: { ...cellOpts(C.white), bold: true, color: C.muted } },
        { text: "4", options: { ...cellOpts(C.white), align: "center", bold: true } },
        { text: "Solvent / LC-MS background (Jul 2–3)", options: cellOpts(C.white) },
        { text: "Blank filter: retain max group mean > 3× mean Blank", options: cellOpts(C.white) },
      ],
      [
        { text: "Pooled QCs (Polled_Qcs)", options: { ...cellOpts(C.tintTeal), bold: true, color: C.violet } },
        { text: "4", options: { ...cellOpts(C.tintTeal), align: "center", bold: true } },
        { text: "Pooled across ALL 12 study drugs (incl. erythromycin!)", options: cellOpts(C.tintTeal) },
        { text: "Diagnostic only — NOT presence filter (~12× dilution!)", options: cellOpts(C.tintTeal) },
      ],
      [
        { text: "QC-Mix (6 standards)", options: { ...cellOpts(C.white), bold: true, color: C.amber } },
        { text: "4", options: { ...cellOpts(C.white), align: "center", bold: true } },
        { text: "Synthetic sulfonamides + amitriptyline (332 QC-only)", options: cellOpts(C.white) },
        { text: "RT & mass-accuracy check; excluded from bio statistics", options: cellOpts(C.white) },
      ],
    ];
    // Total table height = 0.40 + 5*0.56 = 3.20" -> y=3.50 to y=6.70!
    s.addTable(rows, {
      x: 0.55,
      y: 3.50,
      w: 8.15,
      colW: [1.92, 0.90, 2.48, 2.85],
      rowH: [0.40, 0.56, 0.56, 0.56, 0.56, 0.56],
      border: { pt: 0.6, color: C.border },
    });

    // Right cards aligned top (3.50) and bottom (6.70) identically to the table!
    addAccentCard(pres, s, {
      x: 9.00, y: 3.50, w: 3.78, h: 1.66, accentColor: C.omep,
      title: "Why 4,525 vs 6,993 Features?",
      titleSize: 11.4,
      bodySize: 9.9,
      body: "The handout cites 6,993 features, but row ID spans 16–11,816: this CSV is a subset of the 12-drug ChemProp2 table. Pooled QCs still carry 170 features absent from our samples (e.g., erythromycin ID 4135) — so 629/1,060 omeprazole-only ions are missing in Pooled QCs.",
    });

    addAccentCard(pres, s, {
      x: 9.00, y: 5.42, w: 3.78, h: 1.28, accentColor: C.dmso,
      title: "Metadata Quirks Handled Cleanly",
      titleSize: 11.4,
      bodySize: 9.9,
      body: "Column typo ATRRIBUTE_SampleType (3 R's); lowercase t7_*.mzML filenames; Polled_Qcs naming; missing seconds in T3_Omepr_rep2 run_date (parsed via ISO8601); trailing empty column dropped.",
    });

    s.addNotes(
      "Tóki speaking. This slide resolves two major early design questions. First, why does the quant table have 4,525 rows instead of 6,993? Row IDs span 16 to 11,816 because this is a subset of Pakkir Shah et al.'s 12-drug study. That immediately answers how to use our controls: Pooled QCs were pooled across all 12 drugs—we literally detected erythromycin in them—which dilutes omeprazole-specific features 12-fold so 629 of 1,060 omeprazole-only features are absent from Pooled QCs. Thus we filter on replicate prevalence plus >3x Blank, use T0 Medium to separate mGAM background, and treat QC-Mix strictly as synthetic standards."
    );
  }

  // =========================================================================
  // SLIDE 4 — QC AUDIT: FAILED INJECTIONS & RUN-ORDER CONFOUNDING (Tóki)
  // =========================================================================
  {
    const s = pres.addSlide();
    addChrome(pres, s, {
      section: "03 · Quality Control & Confounder Audit (Week 4 & Week 8)",
      title: "Four Failed Injections & Queue Decay of Parent Omeprazole",
      speakerTag: "TÓKI · DESIGN & QC",
      speakerColor: C.dmso,
      slideNum: 4,
    });

    await addFigureCard(pres, s, "fig01_qc_failed_injections_and_run_order.png", {
      x: 0.55, y: 1.14, w: 12.23, h: 4.30, tightFrame: true,
    });

    addAccentCard(pres, s, {
      x: 0.55, y: 5.72, w: 3.88, h: 1.12, accentColor: C.crimson,
      title: "1. Four Blank-Like Failed Injections",
      titleSize: 10.8,
      bodySize: 9.7,
      body: "T4/T5 rep3 & T5_Omepr_rep2 have TIC 4.4–4.7×10^{7} (< all Blanks!), >81% zeros, r = 0.83–0.89 with Blanks (3 back-to-back Jul 4). Removed → n = 50.",
    });

    addAccentCard(pres, s, {
      x: 4.72, y: 5.72, w: 3.88, h: 1.12, accentColor: C.amber,
      title: "2. Run-Order Confounding (Jul 2–7)",
      titleSize: 10.8,
      bodySize: 9.7,
      body: "Timepoint correlates negatively with 5.3-day queue time (r = −0.47 DMSO, −0.43 Omepr): T4/T5/T8 ran Days 1–2 while T1/T2/T6 ran Days 5–6.",
    });

    addAccentCard(pres, s, {
      x: 8.90, y: 5.72, w: 3.88, h: 1.12, accentColor: C.dmso,
      title: "3. Acid-Labile Parent vs. Stable Product",
      titleSize: 10.8,
      bodySize: 9.7,
      body: "Parent omeprazole (ID 2828) decays in queue (r = −0.62, p = 9.4×10^{-4}) hitting 0 at T1 & T2! Microbial ufiprazole (ID 2962) tracks bio time (r = +0.98).",
    });

    s.addNotes(
      "Tóki speaking. Auditing sample-level QC revealed two critical confounders. Panel A & B: four bacterial injections at T4 and T5 failed completely—their TIC is lower than solvent blanks and they correlate 0.83 to 0.89 with blanks versus 0.21 to 0.40 with their own biological replicates; three occurred back-to-back within 30 minutes on July 4. Removing them leaves 50 valid samples. Panel C: over the 5.3-day acquisition queue, timepoint is confounded with run order (r = -0.47), and acid-labile parent omeprazole decays with autosampler wait time (r = -0.62), dropping to zero in T1 and T2 which ran last on days 5–6. Meanwhile ufiprazole tracks true biological incubation time cleanly (r = +0.98)."
    );
  }

  // =========================================================================
  // SLIDE 5 — NORMALIZATION TRAP (PQN vs TIC) & ESI ION SUPPRESSION (Tóki)
  // =========================================================================
  {
    const s = pres.addSlide();
    addChrome(pres, s, {
      section: "04 · Filtering, Normalization & Matrix Effects (Week 4 & Week 8)",
      title: "PQN Beats TIC — Plus Co-Elution ESI Ion Suppression",
      speakerTag: "TÓKI · PROCESSING",
      speakerColor: C.dmso,
      slideNum: 5,
    });

    await addFigureCard(pres, s, "fig02_normalization_pqn_vs_tic_and_ion_suppression.png", {
      x: 0.55, y: 1.14, w: 12.23, h: 4.30, tightFrame: true,
    });

    addAccentCard(pres, s, {
      x: 0.55, y: 5.72, w: 5.95, h: 1.12, accentColor: C.omep,
      title: "Compositional Trap: TIC Fabricates ~1,300 False Down Hits",
      titleSize: 11.0,
      bodySize: 9.8,
      body: "Ufiprazole + its 9 co-eluting IIMN ions (cg=220) surge 9× to ~47% of total TIC, doubling late Omeprazole TIC (1.05× → 1.99×). Blind TIC pushes 83% (1,780/2,146) of shared features down; PQN stays flat (~1.0×).",
    });

    addAccentCard(pres, s, {
      x: 6.83, y: 5.72, w: 5.95, h: 1.12, accentColor: C.violet,
      title: "Co-Elution ESI Ion Suppression Notch (RT 1.80–2.25 min, OR = 46.5)",
      titleSize: 11.0,
      bodySize: 9.8,
      body: "Even after PQN, 114/172 (66.3%) of negative shared Trt×Time hits cluster under Ufiprazole's RT 1.80–2.25 min peak (only 9% of features; Fisher p = 2.1×10^{-91}). Flagging this window prevents false depletion calls.",
    });

    s.addNotes(
      "Tóki speaking. After filtering to 3,414 features (>=2 replicates in >=1 group, >3x Blank), we discovered two major LC-MS measurement pitfalls. First (Panels A-B): because ufiprazole and its in-source fragments double late omeprazole TIC, standard TIC normalization falsely depresses 83% of shared endogenous features. Switching to Probabilistic Quotient Normalization (PQN) anchored on 1,902 ubiquitous DMSO features eliminates this compositional distortion. Second (Panel C): plotting interaction slopes against retention time exposed a razor-sharp notch right under ufiprazole's massive 4.9e8 peak at RT 1.90 min—66.3% of all negative shared interactions occur in just 9% of RT space (OR = 46.5, p = 2.1e-91). This is classic ESI charge-competition ion suppression, so we explicitly flag RT 1.80–2.25 min before interpreting endogenous biology."
    );
  }

  // =========================================================================
  // SLIDE 6 — PCA & MULTIVARIATE PERMANOVA (Magda)
  // =========================================================================
  {
    const s = pres.addSlide();
    addChrome(pres, s, {
      section: "05 · Unsupervised Ordination & Variance Partitioning (Week 4 & Week 9)",
      title: "Global Treatment Split vs. Conserved Endogenous Growth",
      speakerTag: "MAGDA · STATS & ML",
      speakerColor: C.omep,
      slideNum: 6,
    });

    await addFigureCard(pres, s, "fig03_pca_and_variance_partitioning.png", {
      x: 0.55, y: 1.14, w: 12.23, h: 4.30, tightFrame: true,
    });

    addAccentCard(pres, s, {
      x: 0.55, y: 5.72, w: 3.88, h: 1.12, accentColor: C.omep,
      title: "Global PCA (3,414 Features)",
      titleSize: 10.8,
      bodySize: 9.7,
      body: "PC1 (32.6%) separates Treatment (50.8% R^{2}); PC2 (18.7%) traces T0→T8 Time (28.4%). Interaction 11.7% (all p = 0.001; residual 9.0%).",
    });

    addAccentCard(pres, s, {
      x: 4.72, y: 5.72, w: 3.88, h: 1.12, accentColor: C.dmso,
      title: "Clean Endogenous PCA (n = 1,915)",
      titleSize: 10.8,
      bodySize: 9.7,
      body: "Outside RT 1.80–2.25 suppression, Time dominates at 59.1% (Trt 16.1%, Int 7.1%, p = 0.022). Both arms move in parallel away from T0 Medium.",
    });

    addAccentCard(pres, s, {
      x: 8.90, y: 5.72, w: 3.88, h: 1.12, accentColor: C.green,
      title: "Run-Order Confounder Verified",
      titleSize: 10.8,
      bodySize: 9.7,
      body: "Conditional on Treatment × Time, residual acquisition run order explains ≤ 1.1% (p > 0.06, n.s.) — confirming biological variance dominates.",
    });

    s.addNotes(
      "Magda speaking. Turning to multivariate analysis on PQN-normalized log2 intensities: Panel A shows global PCA across all 3,414 features—PC1 separates treatment (50.8% PERMANOVA R²) while PC2 tracks biological time from T0 to T8 (28.4%), with 91% of total variance explained. Crucially, Panel B restricts PCA to the 1,915 clean shared endogenous features outside the ion-suppression window plus T0 Medium controls: now biological time dominates at 59.1% while treatment drops to 16.1%, and both DMSO and omeprazole cultures progress along parallel growth trajectories departing from T0 medium. Conditional run-order variance is under 1.1% and non-significant."
    );
  }

  // =========================================================================
  // SLIDE 7 — DIFFERENTIAL STATISTICS & NON-LINEAR MODELS (Magda)
  // =========================================================================
  {
    const s = pres.addSlide();
    addChrome(pres, s, {
      section: "06 · Univariate Statistical Inference (Week 4)",
      title: "Linear & Categorical ANOVA Capture Slopes and Transients",
      speakerTag: "MAGDA · STATS & ML",
      speakerColor: C.omep,
      slideNum: 7,
    });

    await addFigureCard(pres, s, "fig04_statistical_volcano_and_model_comparison.png", {
      x: 0.55, y: 1.14, w: 12.23, h: 4.30, tightFrame: true,
    });

    addAccentCard(pres, s, {
      x: 0.55, y: 5.72, w: 3.88, h: 1.12, accentColor: C.dmso,
      title: "Run-Adjusted Linear Model (BH-FDR)",
      titleSize: 10.8,
      bodySize: 9.7,
      body: "y ~ Trt + Time + Trt:Time + RunOrder (q < 0.05): 2,064 Treatment, 807 Time, 904 Trt×Time, and 526 queue-drift hits cleanly isolated.",
    });

    addAccentCard(pres, s, {
      x: 4.72, y: 5.72, w: 3.88, h: 1.12, accentColor: C.omep,
      title: "Categorical ANOVA (+774 Transients)",
      titleSize: 10.8,
      bodySize: 9.7,
      body: "y ~ Trt * factor(Time) catches 1,595 Trt×Time hits — rescuing +774 non-linear cascade intermediates peaking at T3–T5 (e.g., ID 4470, 2719).",
    });

    addAccentCard(pres, s, {
      x: 8.90, y: 5.72, w: 3.88, h: 1.12, accentColor: C.violet,
      title: "Model Awareness: Log_{2} Carryover Trap",
      titleSize: 10.8,
      bodySize: 9.7,
      body: "For Ufiprazole (ID 2962), 0.01% DMSO carryover fluctuates 7k→80k (+3.5 log_{2}), masking its log_{2} interaction despite a 6,100× ratio (within-Omepr r = +0.98).",
    });

    s.addNotes(
      "Magda speaking. For univariate inference across all 3,414 features with Benjamini-Hochberg FDR control (q < 0.05), we paired two complementary models. First, a run-order-adjusted linear model estimates directional treatment and slope differences while controlling for queue drift. Second, a two-way categorical ANOVA captures +774 non-linear treatment-by-time hits—transient biotransformation intermediates like ID 4470 and ID 2719 that rise to T3 and then deplete by T8, which a straight-line slope misses. We also caught a subtle log-scale artifact on ufiprazole itself: tiny 0.01% column carryover in DMSO fluctuates on a log2 scale almost as much as real accumulation in omeprazole, masking its log2 interaction term—which is why we also run within-treatment temporal tests."
    );
  }

  // =========================================================================
  // SLIDE 8 — TEMPORAL HEATMAP & 4 KINETIC MODULES (Magda) — FULL WIDTH FOR LEGIBILITY!
  // =========================================================================
  {
    const s = pres.addSlide();
    addChrome(pres, s, {
      section: "07 · Temporal Clustering & Cascade Dynamics (Week 9)",
      title: "606 Unique Dynamic Entities Split Into 4 Kinetic Waves",
      speakerTag: "MAGDA · STATS & ML",
      speakerColor: C.omep,
      slideNum: 8,
    });

    await addFigureCard(pres, s, "fig05_temporal_heatmap_and_kinetic_clusters.png", {
      x: 0.55, y: 1.12, w: 12.23, h: 4.52, tightFrame: true,
    });

    addAccentCard(pres, s, {
      x: 0.55, y: 5.84, w: 3.88, h: 1.04, accentColor: C.primary,
      title: "IIMN Ion-Identity Collapse",
      titleSize: 10.6,
      bodySize: 9.5,
      body: "3,414 filtered ions collapse to 1,814 unique entities → 606 omeprazole-dynamic entities clustered by k-means (prevents Ufiprazole's 10 ions dominating).",
    });

    addAccentCard(pres, s, {
      x: 4.72, y: 5.84, w: 3.88, h: 1.04, accentColor: C.violet,
      title: "Sequential Cascade Waves (M4 → M3 → M1)",
      titleSize: 10.6,
      bodySize: 9.5,
      body: "Early depleted precursors (M4, n=84) feed T2–T3 intermediates (M3, n=204: C_{17}H_{17}N_{3}O_{4} ID 4470) and monotonic end-products (M1, n=215).",
    });

    addAccentCard(pres, s, {
      x: 8.90, y: 5.84, w: 3.88, h: 1.04, accentColor: C.amber,
      title: "M2 Captures Queue-Labile Parent (2828*)",
      titleSize: 10.6,
      bodySize: 9.5,
      body: "M2 (n=103, peak T4–T5) includes mid-late transients plus acid-labile parent omeprazole (ID 2828), whose dip at T1–T2 reflects Day 5–6 queue decay.",
    });

    s.addNotes(
      "Magda speaking. Before clustering, we collapsed co-eluting IIMN correlation groups from 3,414 features down to 1,814 unique molecular entities so that a single molecule like ufiprazole (10 features) wouldn't distort cluster sizes. Giving this figure full slide width makes all 19 landmark heatmap rows crisp and legible. K-means clustering of the 606 dynamic omeprazole-enriched entities reveals four sequential kinetic waves: 84 early depleted precursors, 204 early-to-mid cascade intermediates peaking at T2–T3 and getting consumed by T8, 103 mid-late transients peaking at T4–T5 (including queue-labile parent omeprazole), and 215 monotonic end-products—hallmarks of multi-step microbial biotransformation."
    );
  }

  // =========================================================================
  // SLIDE 9 — SUPERVISED ML & DARK METABOLOME PRIORITIZATION (Magda)
  // =========================================================================
  {
    const s = pres.addSlide();
    addChrome(pres, s, {
      section: "08 · Supervised Machine Learning & Unknown Ranking (Week 9)",
      title: "ML Classifiers Confirm Endogenous Signal & Rank Unknowns",
      speakerTag: "MAGDA · STATS & ML",
      speakerColor: C.omep,
      slideNum: 9,
    });

    await addFigureCard(pres, s, "fig06_ml_prioritization_known_and_unknown.png", {
      x: 0.55, y: 1.14, w: 12.23, h: 4.30, tightFrame: true,
    });

    addAccentCard(pres, s, {
      x: 0.55, y: 5.72, w: 3.88, h: 1.12, accentColor: C.dmso,
      title: "Perfect Treatment & Stage Separation",
      titleSize: 10.8,
      bodySize: 9.7,
      body: "5-fold stratified CV achieves AUC = 1.000 (RF & L1-LogReg) on all unique features AND on clean endogenous features alone—confirming real community shifts.",
    });

    addAccentCard(pres, s, {
      x: 4.72, y: 5.72, w: 3.88, h: 1.12, accentColor: C.omep,
      title: "#1 Ranked Unknown Matches Source Paper",
      titleSize: 10.8,
      bodySize: 9.7,
      body: "Composite ranking (effect size + FDR + RF importance + consistency) puts unannotated ID 2618 (m/z 360.14, Me-omeprazole from Pakkir Shah Fig. 1!) at #1.",
    });

    addAccentCard(pres, s, {
      x: 8.90, y: 5.72, w: 3.88, h: 1.12, accentColor: C.green,
      title: "Dark Metabolome Dominates Top Hits",
      titleSize: 10.8,
      bodySize: 9.7,
      body: "14 of the top 15 discriminant features lack MS2 library names (ID 2678 +CH_{2} of 2368, ID 1553 r = +0.99, ID 2368 2S)—showing why unknowns must be retained!",
    });

    s.addNotes(
      "Magda speaking. Following the dataset handout's recommendation to retain and rank unknown features alongside known ones, we trained stratified 5-fold cross-validated Random Forest and regularized logistic regression models. Both separate omeprazole from DMSO with AUC 1.000—even when restricted strictly to clean endogenous features outside the suppression zone—and classify bacterial growth stages with 96–100% balanced accuracy. Ranking IIMN-collapsed features by a composite score of effect size, FDR, RF importance, and replicate consistency surfaces 14 unannotated dark-metabolome candidates among the top 15 hits—including #1 ID 2618 (m/z 360.14, the exact methylated omeprazole highlighted in Pakkir Shah et al. Figure 1)."
    );
  }

  // =========================================================================
  // SLIDE 10 — ANNOTATION AUDIT, ORA & ENDOGENOUS BIOLOGY (Caroline) — FULL WIDTH FIG!
  // =========================================================================
  {
    const s = pres.addSlide();
    addChrome(pres, s, {
      section: "09 · Interpretation, Enrichment & Endogenous Pathways (Weeks 6–7 & 13)",
      title: "Chemical-Class Enrichment & Selective Endogenous Shifts",
      speakerTag: "CAROLINE · ENRICHMENT",
      speakerColor: C.green,
      slideNum: 10,
    });

    // Full-width tight-fitted fig07 so all 6 subplots & ORA labels are crystal-clear!
    await addFigureCard(pres, s, "fig07_chemical_class_enrichment_and_endogenous_trajectories.png", {
      x: 0.55, y: 1.10, w: 12.23, h: 4.56, tightFrame: true,
    });

    addAccentCard(pres, s, {
      x: 0.55, y: 5.82, w: 3.88, h: 1.08, accentColor: C.crimson,
      title: "1. Annotation Audit (259 / 3,414 = 7.6%)",
      titleSize: 10.6,
      bodySize: 9.5,
      body: "RDKit curation removed trace carryover from 10 other study drugs + 6 QC-Mix standards (n = 228 bio background) & flagged false hit 'Sertindole'.",
    });

    addAccentCard(pres, s, {
      x: 4.72, y: 5.82, w: 3.88, h: 1.08, accentColor: C.green,
      title: "2. Fisher's Exact ORA (BH-FDR q < 0.05)",
      titleSize: 10.6,
      bodySize: 9.5,
      body: "Omepr: Benzimidazoles 10.1× (q=2.2e−11). Produced: Bile Acids 5.5× & N-Acyl AAs 2.9× (q=0.012). Stable mGAM: DKPs 3.7× (q=1.1e−18).",
    });

    addAccentCard(pres, s, {
      x: 8.90, y: 5.82, w: 3.88, h: 1.08, accentColor: C.dmso,
      title: "3. Clean Endogenous vs. Artifact",
      titleSize: 10.6,
      bodySize: 9.5,
      body: "Conserved: N-Ac-histamine, Gly-C17:1, oleic acid. Genuine shifts: IPA (ID 4553, RT 2.60) 2× attenuated & PE 13:0 delayed — unlike suppressed ID 2972.",
    });

    s.addNotes(
      "Caroline speaking. To move from a list of hits to biological pathways, we first audited all 259 kept MS2 annotations using RDKit SMARTS and chemical curation, separating out trace carryover of 10 other drugs from the 12-drug batch and 6 QC-Mix standards so our Fisher's exact background reflects the 228 measured biological compounds. Over-representation analysis confirms three crisp biological themes: omeprazole-enriched features are 10.1-fold enriched for benzimidazole/prazole derivatives (q = 2.2e-11); bacterial growth over T0–T8 actively produces bile acids (5.5x, bile salt hydrolase deconjugation) and microbial N-acyl amino acids (2.9x); and abundant cyclic dipeptides (3.7x, q = 1.1e-18) represent stable mGAM medium background. On the right, genuine trajectories like N-acetylhistamine, Gly-C17:1, and oleic acid are conserved, while indole-3-propionic acid (RT 2.60) is selectively reduced ~2-fold and PE 13:0 is delayed—contrasted directly with N-phenethylacetamide (RT 1.87), which drops to zero purely from ufiprazole ion suppression."
    );
  }

  // =========================================================================
  // SLIDE 11 — MULTI-LAYER NETWORK ANALYSIS (Irene)
  // =========================================================================
  {
    const s = pres.addSlide();
    addChrome(pres, s, {
      section: "10 · Network Analysis: Hubs, Modules & Cascades (Week 3 & Week 13)",
      title: "Biotransformation & Endogenous Co-Abundance Networks",
      speakerTag: "IRENE · NETWORKS",
      speakerColor: C.violet,
      slideNum: 11,
    });

    await addFigureCard(pres, s, "fig08_biotransformation_and_coabundance_networks.png", {
      x: 0.55, y: 1.16, w: 8.65, h: 5.06, tightFrame: false,
    });

    s.addShape(pres.shapes.RECTANGLE, { x: 0.55, y: 6.42, w: 8.65, h: 0.44, fill: { color: C.tintTeal }, line: { color: C.border, width: 0.6 } });
    s.addText("Cytoscape GraphML exports: results/networks/omeprazole_biotransformation_network.graphml & endogenous_coabundance_network.graphml", {
      x: 0.68, y: 6.42, w: 8.40, h: 0.44, fontFace: FONT_BODY, fontSize: 8.8, italic: true, color: C.primary, valign: "middle", margin: 0,
    });

    const bioCardRuns = [
      ...chemRuns("Full graph: 57 IIMN-collapsed nodes, 82 exact-mass (< 4.5 mDa) edges (Panel A displays the 34-node, 76-edge connected core).", { fontFace: FONT_BODY, fontSize: 10.0, color: C.text }),
      { text: "Sub-0.1 mDa enzymatic shifts: ", options: { fontFace: FONT_BODY, fontSize: 10.0, bold: true, color: C.primary } },
      ...chemRuns("±O sulfoxide reduction (Omeprazole ↔ Ufiprazole; Me-Omeprazole 2618 ↔ 3605/2500), ±CH_{2} methylation (2368 ↔ 2678, r = 0.996; 2480 ↔ 1677), ±H_{2} redox, and C–S cleavage.", { fontFace: FONT_BODY, fontSize: 10.0, color: C.text }),
    ];
    bioCardRuns[0].options.paraSpaceAfter = 8;
    bioCardRuns[0].options.breakLine = true;

    addAccentCard(pres, s, {
      x: 9.50, y: 1.16, w: 3.28, h: 2.70, accentColor: C.omep,
      title: "A. Xenobiotic Transformation Graph",
      titleSize: 11.2,
      bodySize: 10.0,
      lineSpacing: 1.20,
      body: bioCardRuns,
    });

    addAccentCard(pres, s, {
      x: 9.50, y: 4.16, w: 3.28, h: 2.70, accentColor: C.dmso,
      title: "B. Endogenous Com20 Core Network",
      titleSize: 11.2,
      bodySize: 10.0,
      lineSpacing: 1.20,
      body: [
        { text: "91 annotated metabolites, 834 edges (|r| ≥ 0.72, q < 0.001), 8 communities (Q = 0.224; Panel B shows 58-node giant core at |r| ≥ 0.78).", options: { color: C.text, paraSpaceAfter: 7, breakLine: true } },
        { text: "Growth Hubs: ", options: { bold: true, color: C.omep, paraSpaceAfter: 6 } },
        { text: "15-HEDE (deg 47), LPE 13:0 & Phe-C4:0 (44), IPA, N-Ac-histamine & Oleic acid (42).", options: { color: C.text, breakLine: true } },
        { text: "Bridge Hub: ", options: { bold: true, color: C.violet } },
        { text: "Guanine (highest betweenness = 0.285), linking uptake & synthesis modules.", options: { color: C.text } },
      ],
    });

    s.addNotes(
      "Irene speaking. Applying Week 3 and Week 13 network analysis, we built two complementary graphs—both exported as Cytoscape GraphML files. Panel A shows the 34-node, 76-edge connected core of our 57-node ChemProp-inspired omeprazole biotransformation mass-difference network: parent omeprazole and ufiprazole sit at the center of dense ±O, ±CH2, ±H2, and C-S cleavage cascades linking both annotated products and high-priority unknowns like ID 2368 and ID 2678 (0.01 mDa error, r = 0.996). Panel B shows the giant connected core of the endogenous Com20 co-abundance network (91 total nodes, 834 edges, 8 modules), cleanly separating nutrient-uptake blue hubs like oleic acid and Phe-Pro from bacterially produced red hubs like LPE 13:0, Phe-C4:0, N-acetylhistamine, and indole-3-propionic acid, bridged by guanine."
    );
  }

  // =========================================================================
  // SLIDE 12 — MOLECULAR ZOOM-IN: STRUCTURES, REGIOISOMERS & ISOTOPES (Irene)
  // =========================================================================
  {
    const s = pres.addSlide();
    addChrome(pres, s, {
      section: "11 · Molecular Zoom-In & Chemical Verification (Week 2 & Week 8)",
      title: "Zoom-In: Regioisomers, 2S Conjugate & False-Hit Audit",
      speakerTag: "IRENE · HIT ZOOM-IN",
      speakerColor: C.violet,
      slideNum: 12,
    });

    await addFigureCard(pres, s, "fig09_hit_zoomin_structures_and_regioisomers.png", {
      x: 0.55, y: 1.16, w: 8.45, h: 5.70, tightFrame: false,
    });

    addAccentCard(pres, s, {
      x: 9.28, y: 1.16, w: 3.50, h: 1.78, accentColor: C.omep,
      title: "1. Ufiprazole Cascade & Multiplicity",
      titleSize: 10.8,
      bodySize: 9.3,
      lineSpacing: 1.24,
      body: "Bacterial sulfoxide reduction (Watanabe et al., 1995) converts Omeprazole (346.12) → Ufiprazole (ID 2962, 330.13, +9×). One molecule spawns 10 IIMN ions in cg=220 ([M+H]^{+}, ^{13}C, ^{34}S, [M+Na]^{+}, dimer + 180.04 & 151.10 frags), branching to 4-desmethyl (2465) & 5-MeO-2-MBI (2480, ^{34}S 4.2%).",
    });

    const regIntro = chemRuns("Two methyl-ufiprazole isomers (344.1427 Da, 0.0 ppm; C_{18}H_{21}N_{3}O_{2}S) localized via IIMN fragments:", { fontFace: FONT_BODY, fontSize: 9.2, color: C.text });
    regIntro[0].options.paraSpaceAfter = 5;
    regIntro[regIntro.length - 1].options.breakLine = true;

    addAccentCard(pres, s, {
      x: 9.28, y: 3.16, w: 3.50, h: 1.80, accentColor: C.violet,
      title: "2. Regioisomer Proof via IIMN Frags",
      titleSize: 10.8,
      bodySize: 9.2,
      lineSpacing: 1.24,
      body: [
        ...regIntro,
        ...makeChemBullets([
          "ID 3605 (RT 2.17, cg=264): 194.051 (+CH_{2} benzimid.) + 151.099 (intact pyr.) → Benzimidazole-Me*",
          "ID 2500 (RT 1.55, cg=175): 180.035 (benzimid.) + 165.115 (+CH_{2} pyr.) → Pyridine-Me*",
        ], 9.2, 5),
      ],
    });

    addAccentCard(pres, s, {
      x: 9.28, y: 5.18, w: 3.50, h: 1.68, accentColor: C.crimson,
      title: "3. Dark Hit ID 2368 & 'Sertindole' Audit",
      titleSize: 10.8,
      bodySize: 9.3,
      lineSpacing: 1.24,
      body: makeChemBullets([
        "ID 2368 (274.0565, +95×): 8.2% ^{34}S M+2 proves 2S + pyridine 151.0992 → unique C_{11}H_{15}NO_{3}S_{2} (−0.46 ppm; putative thioglycolate disulfide*).",
        "ID 5664 ('Sertindole'): 0.0% M+2 vs 32% required for ^{37}Cl → refutes false library hit!",
      ], 9.3, 6),
    });

    s.addNotes(
      "Irene speaking. Zooming in structurally with RDKit (Panel A & C) and high-resolution isotope physics (Panel D) gives us mechanistic confidence. First, ufiprazole (ID 2962) undergoes both O-demethylation (ID 2465) and C-S cleavage to 5-methoxy-2-mercaptobenzimidazole (ID 2480, 34S = 4.2% confirming 1 sulfur). Second, we regio-localized two abundant methylated ufiprazole isomers (m/z 344.1427) purely from their co-eluting IIMN fragments: ID 3605 carries the +14 Da methyl on the benzimidazole fragment (194.051 + 151.099), whereas ID 2500 carries it on the pyridine fragment (180.035 + 165.115). Third, our top unannotated hit ID 2368 (+95-fold) has an 8.2% 34S M+2 peak confirming TWO sulfurs plus the intact pyridine fragment 151.0992, uniquely solving to C11H15NO3S2 (-0.46 ppm)—consistent with a mixed disulfide with mGAM thioglycolate (a structural hypothesis awaiting MS2 standard verification). Finally, Panel D proves why we never trust library names blindly: ID 5664 matched chlorinated Sertindole by MS2 cosine, but has 0% M+2 instead of 32% 37Cl."
    );
  }

  // =========================================================================
  // SLIDE 13 — MEASUREMENT CONSTRAINTS, LIMITATIONS & COURSE TOOLKIT
  // =========================================================================
  {
    const s = pres.addSlide();
    addChrome(pres, s, {
      section: "12 · Model Awareness & Measurement Constraints",
      title: "LC-MS/MS Measurement Limits & Bounded Conclusions",
      speakerTag: "LIMITATIONS · ALL",
      speakerColor: C.primary,
      slideNum: 13,
    });

    // ALL 6 titles fit crisply on 1 single line! Uniform 0.30" gaps!
    const limits = [
      {
        icon: ic.micro, col: C.dmso,
        t: "1. ESI+ & Ethyl Acetate Bias",
        b: "Favors basic/neutral N-heterocycles (prazoles, indoles, DKPs) and lipids; blind to anionic short-chain fatty acids (acetate, propionate, butyrate) and polar sugar-phosphates.",
      },
      {
        icon: ic.scale, col: C.omep,
        t: "2. Relative Intensity ≠ Concentration",
        b: "Peak areas reflect ionization efficiency, not molarity. High-ionizing ufiprazole accounts for ~47% of late TIC — requiring PQN normalization instead of TIC.",
      },
      {
        icon: ic.warn, col: C.violet,
        t: "3. Co-Elution ESI Ion Suppression",
        b: "Overloaded ufiprazole peak suppresses 66.3% of negative shared interactions inside RT 1.80–2.25 min (OR = 46.5), masking endogenous trends in that 0.45-min slice.",
      },
      {
        icon: ic.clock, col: C.amber,
        t: "4. Queue Decay & Abiotic Control",
        b: "Acid-labile parent omeprazole degraded over the 5.3-day queue (T1/T2 = 0). No abiotic Omeprazole + mGAM time series is in this CSV (source paper confirms microbial enrichment).",
      },
      {
        icon: ic.search, col: C.crimson,
        t: "5. 92.4% Dark Metabolome",
        b: "Only 259/3,414 features match MS2 libraries, and cosine hits can misassign novel derivatives (e.g., Sertindole lacking ^{37}Cl). Putative structures await authentic standards.",
      },
      {
        icon: ic.layers, col: C.green,
        t: "6. Power (n = 2–3) & No 16S Here",
        b: "Removing 4 failed injections reduces T4 & T5 to n = 2 per arm, limiting power for subtle endogenous shifts; metabolomics alone cannot assign reactions to single Com20 strains.",
      },
    ];

    limits.forEach((lm, idx) => {
      const colIdx = idx % 3;
      const rowIdx = Math.floor(idx / 3);
      const lx = 0.55 + colIdx * 4.19;
      const ly = 1.24 + rowIdx * 2.28;
      s.addShape(pres.shapes.RECTANGLE, { x: lx, y: ly, w: 3.85, h: 1.98, fill: { color: C.white }, line: { color: C.border, width: 0.8 }, shadow: makeShadow(0.08) });
      s.addShape(pres.shapes.RECTANGLE, { x: lx, y: ly, w: 0.07, h: 1.98, fill: { color: lm.col }, line: { color: lm.col, width: 0 } });
      s.addShape(pres.shapes.OVAL, { x: lx + 0.18, y: ly + 0.16, w: 0.40, h: 0.40, fill: { color: lm.col }, line: { color: lm.col, width: 0 } });
      s.addImage({ data: lm.icon, x: lx + 0.26, y: ly + 0.24, w: 0.24, h: 0.24 });
      s.addText(lm.t, { x: lx + 0.68, y: ly + 0.16, w: 3.06, h: 0.40, fontFace: FONT_BODY, fontSize: 11.2, bold: true, color: lm.col, valign: "middle", margin: 0 });
      s.addText(chemRuns(lm.b, { fontFace: FONT_BODY, fontSize: 10.2, color: C.text }), {
        x: lx + 0.18, y: ly + 0.66, w: 3.52, h: 1.22, valign: "top", lineSpacingMultiple: 1.22, margin: 0,
      });
    });

    // Clean 2-line Course Toolkit bar with zero dangling middots or double spaces!
    s.addShape(pres.shapes.RECTANGLE, { x: 0.55, y: 5.86, w: 12.23, h: 0.94, fill: { color: C.primary }, shadow: makeShadow(0.08) });
    s.addText([
      { text: "Course Representations Applied (6 Modules; ≥3 Required): ", options: { bold: true, color: "F4A261", fontSize: 10.6 } },
      { text: "Week 2 (RDKit structures & Morgan similarity) · Week 3 & 13 (Biotransformation & co-abundance graphs, hubs, modules)", options: { color: C.white, fontSize: 10.2, breakLine: true, paraSpaceAfter: 4 } },
      { text: "Week 4 (Design, confounders, linear/ANOVA, PERMANOVA, BH-FDR) · Weeks 6–7 (Chemical-class Fisher ORA) · Week 8 (LC-MS QC, IIMN, PQN, isotopes) · Week 9 (PCA, k-means, CV Random Forest & LogReg)", options: { color: C.lightOnDark, fontSize: 10.0 } },
    ], { x: 0.75, y: 5.92, w: 11.83, h: 0.82, fontFace: FONT_BODY, valign: "middle", lineSpacingMultiple: 1.20, margin: 0 });

    s.addNotes(
      "Model awareness and measurement constraints are central to our conclusions. Positive-mode ESI after ethyl acetate extraction is blind to acidic short-chain fatty acids and polar central carbon metabolites, so we cannot claim omeprazole has zero effect on SCFA fermentation. Peak areas are relative intensities, not molar concentrations. We explicitly bounded our endogenous conclusions outside the RT 1.80–2.25 ion-suppression window, noted that parent omeprazole queue decay and the absence of an abiotic omeprazole time series in this CSV require relying on Pakkir Shah et al. to separate spontaneous vs enzymatic hydrolysis, and treated all unannotated structures as hypotheses pending authentic standards. Across the project we integrated six distinct course modules."
    );
  }

  // =========================================================================
  // SLIDE 14 — USE OF AI & VERIFICATION DISCIPLINE
  // =========================================================================
  {
    const s = pres.addSlide();
    addChrome(pres, s, {
      section: "13 · Responsible AI Disclosure (Week 1 & Week 12)",
      title: "AI Accelerated Coding — Chemistry & Audits Verified Every Claim",
      speakerTag: "AI REFLECTION · ALL",
      speakerColor: C.amber,
      slideNum: 14,
    });

    const pillars = [
      {
        icon: ic.robot,
        col: C.dmso,
        head: "1. What AI Was Used For",
        items: [
          "Conceptual design diagrams (bottles & timeline drafts prior to coding).",
          "Drafting and debugging reproducible Python scripts (src/pipeline.py, src/make_figures.py) and executed Jupyter notebooks.",
          "Literature retrieval for experimental context (Pakkir Shah et al., 2026; Watanabe et al., 1995).",
          "Structuring clear, colorblind-safe matplotlib & PptxGenJS layouts.",
        ],
      },
      {
        icon: ic.shield,
        col: C.green,
        head: "2. How We Checked Every Output",
        items: [
          "Re-ran the full pipeline deterministically from raw CSVs (seed = 42); every p-value, R^{2}, OR, and count computed directly from data/raw/.",
          "Checked every chemical claim via exact monoisotopic masses (< 1 ppm), ^{34}S / ^{37}Cl M+2 isotope physics, and co-eluting IIMN fragments.",
          "Cross-checked all paper methods (Com20, mGAM, hourly T0–T8, Q Exactive HF) directly against PMC12889426.",
          "Visually inspected every rendered figure and slide via independent QA.",
        ],
      },
      {
        icon: ic.bug,
        col: C.crimson,
        head: "3. Real Errors Caught & Fixed",
        items: [
          "Caught a variable-shadowing bug in summary JSON where ion-suppression stats printed 2/54 instead of 114/172 (66.3%) — fixed & re-verified!",
          "Caught k-means cluster label inversion by inspecting T0–T8 centroid arrays before naming modules.",
          "Refuted GNPS/AI-surfaced 'Sertindole' hit (ID 5664) by proving 0% M+2 instead of 32% ^{37}Cl.",
          "Clearly separated library-confirmed hits from putative structural hypotheses (*).",
        ],
      },
    ];

    pillars.forEach((p, idx) => {
      const px = 0.55 + idx * 4.19;
      s.addShape(pres.shapes.RECTANGLE, { x: px, y: 1.22, w: 3.85, h: 4.42, fill: { color: C.white }, line: { color: C.border, width: 0.8 }, shadow: makeShadow(0.09) });
      s.addShape(pres.shapes.RECTANGLE, { x: px, y: 1.22, w: 3.85, h: 0.74, fill: { color: p.col }, line: { color: p.col, width: 0 } });
      s.addShape(pres.shapes.OVAL, { x: px + 0.18, y: 1.35, w: 0.48, h: 0.48, fill: { color: "FFFFFF", transparency: 82 }, line: { color: "FFFFFF", transparency: 40, width: 0.8 } });
      s.addImage({ data: p.icon, x: px + 0.28, y: 1.45, w: 0.28, h: 0.28 });
      s.addText(p.head, { x: px + 0.76, y: 1.22, w: 2.98, h: 0.74, fontFace: FONT_BODY, fontSize: 12.6, bold: true, color: C.white, valign: "middle", margin: 0 });

      s.addText(makeChemBullets(p.items, 10.9, 11, C.text), {
        x: px + 0.20,
        y: 2.12,
        w: 3.45,
        h: 3.38,
        valign: "top",
        lineSpacingMultiple: 1.24,
        margin: 0,
      });
    });

    // Crisp bottom principle banner with exact vertical centering & >=0.30" top gap
    s.addShape(pres.shapes.RECTANGLE, { x: 0.55, y: 5.96, w: 12.23, h: 0.84, fill: { color: C.primary }, shadow: makeShadow(0.08) });
    s.addText([
      { text: "Scientific Ownership (Full Log in notes/ai-use.md): ", options: { bold: true, color: "F4A261", fontSize: 11.0 } },
      { text: "AI accelerated engineering and exploration, but every statistic, formula, isotope ratio, and biological claim was independently verified from raw data — we defend our choices and evidence, not the chatbot.", options: { color: C.white, fontSize: 10.6 } },
    ], { x: 0.75, y: 5.96, w: 11.83, h: 0.84, fontFace: FONT_BODY, valign: "middle", lineSpacingMultiple: 1.18, margin: 0 });

    s.addNotes(
      "Following the Week 1 and Week 12 AI guidelines, we logged all AI tool usage in notes/ai-use.md. AI accelerated Python coding, visualization, and literature lookup, but every scientific result was verified independently against the raw data and physical chemistry. Concrete examples of catches during QA: visual inspection caught a Python variable-shadowing bug that printed 2/54 instead of 114/172 in Figure 2C; inspecting k-means centroids fixed inverted cluster labels; and checking M+2 isotopes refuted the false Sertindole library hit. We own and can defend every equation, cutoff, and biological interpretation."
    );
  }

  // =========================================================================
  // SLIDE 15 — REQUIRED LAST SLIDE: CONCLUSIONS & GROUP CONTRIBUTIONS (Dark)
  // =========================================================================
  {
    const s = pres.addSlide();
    s.background = { color: C.darkBg };

    s.addText("14 · CONCLUSION & GROUP CONTRIBUTIONS (REQUIRED LAST SLIDE)", {
      x: 0.60, y: 0.30, w: 9.6, h: 0.22, fontFace: FONT_BODY, fontSize: 9.2, bold: true, color: "F4A261", charSpacing: 1.4, margin: 0,
    });
    s.addText("Take-Home Message & Team Contributions", {
      x: 0.60, y: 0.52, w: 9.6, h: 0.48, fontFace: FONT_HEAD, fontSize: 25, bold: true, color: C.white, margin: 0, valign: "middle",
    });

    // Top-right badge matching rest of deck!
    s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: 10.45, y: 0.40, w: 2.28, h: 0.38, rectRadius: 0.19, fill: { color: C.omep }, line: { color: C.omep, width: 0 } });
    s.addText("CONCLUSION · ALL", {
      x: 10.45, y: 0.40, w: 2.28, h: 0.38, fontFace: FONT_BODY, fontSize: 9.2, bold: true, color: C.white, align: "center", valign: "middle", margin: 0,
    });

    // Crisp Conclusion Banner
    s.addShape(pres.shapes.RECTANGLE, { x: 0.60, y: 1.24, w: 12.13, h: 1.30, fill: { color: C.darkCard }, line: { color: "2E5B73", width: 1 } });
    s.addShape(pres.shapes.RECTANGLE, { x: 0.60, y: 1.24, w: 0.09, h: 1.30, fill: { color: C.omep }, line: { color: C.omep, width: 0 } });
    s.addText([
      { text: "Answer to the Biological Question: ", options: { bold: true, color: "F4A261", fontSize: 12 } },
      { text: "Over 0–8 h, Com20 gut bacteria maintain their core endogenous metabolic growth trajectory (59.1% time variance; bile acid & N-acyl amino acid production) with selective modest shifts (2× attenuated IPA, delayed PE 13:0). The dominant omeprazole-specific divergence (50.8% global variance) is active anaerobic microbial drug metabolism — led by 9× sulfoxide reduction to ufiprazole (ID 2962), O-demethylation, C–S cleavage, regioselective methylation (ID 3605 / ID 2500), and 2S dark-metabolome conjugates (ID 2368) — isolated only after correcting failed injections, queue decay, TIC bias, and RT 1.80–2.25 ESI suppression.", options: { color: C.lightOnDark, fontSize: 11 } },
    ], { x: 0.86, y: 1.30, w: 11.70, h: 1.18, fontFace: FONT_BODY, valign: "middle", lineSpacingMultiple: 1.18, margin: 0 });

    // Clean breathing room before table heading & table border
    s.addText("Group Member Contributions (Replace [sXXXXXX] with DTU student numbers before exam)", {
      x: 0.60, y: 2.86, w: 11.0, h: 0.28, fontFace: FONT_BODY, fontSize: 12, bold: true, color: "9FD4EB", margin: 0,
    });

    const th = { fill: { color: "1E4A5E" }, color: C.white, bold: true, fontFace: FONT_BODY, fontSize: 10.6, valign: "middle" };
    const td = (bg, extra = {}) => ({ fill: { color: bg }, color: C.lightOnDark, fontFace: FONT_BODY, fontSize: 10.2, valign: "middle", ...extra });

    const contribRows = [
      [
        { text: "Group Member", options: th },
        { text: "Student Number", options: { ...th, align: "center" } },
        { text: "Assigned Project Sections & Deliverables", options: th },
        { text: "Slides Presented", options: { ...th, align: "center" } },
      ],
      [
        { text: "Tóki", options: td("132E3C", { bold: true, color: C.white }) },
        { text: "[sXXXXXX]", options: td("132E3C", { align: "center", color: "F4A261", bold: true }) },
        { text: "Design & control audit, failed-injection QC, run-order confounder analysis, blank/prevalence filtering, PQN vs TIC normalization, ESI ion-suppression diagnosis (notebooks/02, fig01–02)", options: td("132E3C") },
        { text: "Slides 3–5", options: td("132E3C", { align: "center", bold: true }) },
      ],
      [
        { text: "Magda", options: td("18384A", { bold: true, color: C.white }) },
        { text: "[sXXXXXX]", options: td("18384A", { align: "center", color: "F4A261", bold: true }) },
        { text: "Global & clean-endogenous PCA, 999-perm PERMANOVA, run-adjusted linear & categorical ANOVA (BH-FDR), IIMN-collapsed k-means kinetic clustering, 5-fold CV ML feature ranking (notebooks/02, fig03–06)", options: td("18384A") },
        { text: "Slides 6–9", options: td("18384A", { align: "center", bold: true }) },
      ],
      [
        { text: "Caroline", options: td("132E3C", { bold: true, color: C.white }) },
        { text: "[sXXXXXX]", options: td("132E3C", { align: "center", color: "F4A261", bold: true }) },
        { text: "MS2 annotation curation & carryover audit, RDKit SMARTS chemical-class ontology, Fisher's exact over-representation enrichment (BH-FDR), endogenous Com20 pathway interpretation (notebooks/03, fig07)", options: td("132E3C") },
        { text: "Slide 10", options: td("132E3C", { align: "center", bold: true }) },
      ],
      [
        { text: "Irene", options: td("18384A", { bold: true, color: C.white }) },
        { text: "[sXXXXXX]", options: td("18384A", { align: "center", color: "F4A261", bold: true }) },
        {
          text: chemRuns(
            "Biotransformation mass-difference network & endogenous co-abundance hubs/modules (Cytoscape GraphML), RDKit 2D hit zoom-in, regioisomer & ^{34}S/^{37}Cl isotope verification, limitations synthesis (notebooks/04, fig08–09)",
            { fontFace: FONT_BODY, fontSize: 10.2, color: C.lightOnDark }
          ),
          options: td("18384A"),
        },
        { text: "Slides 11–13", options: td("18384A", { align: "center", bold: true }) },
      ],
    ];

    // Compact balanced row heights: 0.42 + 4*0.70 = 3.22" -> y=3.32 to y=6.54
    s.addTable(contribRows, {
      x: 0.60,
      y: 3.32,
      w: 12.13,
      colW: [1.45, 1.55, 7.68, 1.45],
      rowH: [0.42, 0.70, 0.70, 0.70, 0.70],
      border: { pt: 0.7, color: "2E5B73" },
    });

    // High-contrast footer + matching slide number 15 / 15
    s.addText("All group members contributed jointly to framing (Slides 1–2), AI verification & disclosure (Slide 14), and oral exam preparation. | Reproducible code: src/pipeline.py & notebooks/02–04", {
      x: 0.60, y: 6.96, w: 10.8, h: 0.26, fontFace: FONT_BODY, fontSize: 9.0, color: C.mutedOnDark, italic: true, margin: 0,
    });
    s.addText("15 / 15", {
      x: 11.73, y: 6.96, w: 1.0, h: 0.26, fontFace: FONT_BODY, fontSize: 9.0, bold: true, color: "F4A261", align: "right", margin: 0,
    });

    s.addNotes(
      "Closing slide summarizing our answer to the biological question and each member's contributions. NOTE FOR GROUP BEFORE EXAM: replace the four [sXXXXXX] placeholders with your DTU student numbers (and remove the parenthetical reminder in the heading)!"
    );
  }

  await pres.writeFile({ fileName: OUT_PATH });
  console.log(`Successfully wrote ${OUT_PATH}`);
}

buildDeck().catch((err) => {
  console.error(err);
  process.exit(1);
});
